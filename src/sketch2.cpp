#include "sketch2.hpp"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fcntl.h>
#include <numeric>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_set>

#include "enc.hpp"
#include "exthash.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "rqseq.hpp"
#include "tpool.hpp"

extern uint32_t num_threads;

struct Sketch2::MappedFile
{
  int fd = -1;
  void* ptr = nullptr;
  size_t len = 0;

  MappedFile() = default;
  MappedFile(const MappedFile&) = delete;
  MappedFile& operator=(const MappedFile&) = delete;

  ~MappedFile()
  {
    if (ptr && ptr != MAP_FAILED && len) munmap(ptr, len);
    if (fd >= 0) close(fd);
  }

  const char* begin() const { return static_cast<const char*>(ptr); }
  const char* end() const { return begin() + len; }

  static std::shared_ptr<MappedFile> open(const std::filesystem::path& path)
  {
    auto m = std::make_shared<MappedFile>();
    m->fd = ::open(path.c_str(), O_RDONLY);
    if (m->fd < 0) error_exit("Cannot open sketch2 file: " + path.string());
    struct stat st
    {};
    if (fstat(m->fd, &st) != 0) error_exit("Cannot fstat sketch2 file: " + path.string());
    m->len = static_cast<size_t>(st.st_size);
    if (m->len == 0) error_exit("Empty sketch2 file: " + path.string());
    m->ptr = mmap(nullptr, m->len, PROT_READ, MAP_PRIVATE, m->fd, 0);
    if (m->ptr == MAP_FAILED) error_exit("mmap failed for sketch2 file: " + path.string());
#ifdef MADV_SEQUENTIAL
    madvise(m->ptr, m->len, MADV_SEQUENTIAL);
#endif
    return m;
  }
};

namespace {

  vec<uint64_t> sample_coords(uint64_t npos, uint64_t nsamples, std::mt19937& rng)
  {
    assert(npos >= 1);
    const uint64_t n = std::min(npos, nsamples);
    if (n < npos / 2) {
      vec<uint64_t> out;
      out.reserve(n);
      std::unordered_set<uint64_t> seen;
      seen.reserve(static_cast<size_t>(n) * 2);
      std::uniform_int_distribution<uint64_t> pick(0, npos - 1);
      while (out.size() < n) {
        const uint64_t x = pick(rng);
        if (seen.insert(x).second) out.push_back(x);
      }
      return out;
    }
    vec<uint64_t> out(npos);
    std::iota(out.begin(), out.end(), uint64_t(0));
    std::shuffle(out.begin(), out.end(), rng);
    out.resize(n);
    return out;
  }

  // Emit packed (bix, enc) hashes for every frac-retained k-mer in [j0, j1).
  template<bool StrandAware>
  void collect_range_hashes(const char* cseq,
                            uint64_t len,
                            uint64_t j0,
                            uint64_t j1,
                            const LSHF& lshf,
                            uint32_t nrows,
                            bool canonical,
                            uint64_t mask_bp,
                            uint64_t mask_lr,
                            vec<uint64_t>& out_fw,
                            vec<uint64_t>& out_rc,
                            uint32_t& nvalid_fw,
                            uint32_t& nvalid_rc)
  {
    const uint32_t k = lshf.get_k();
    out_fw.clear();
    if constexpr (StrandAware) out_rc.clear();
    nvalid_fw = 0;
    nvalid_rc = 0;
    if (j1 <= j0 || len < k) return;
    out_fw.reserve(static_cast<size_t>(j1 - j0));
    if constexpr (StrandAware) out_rc.reserve(static_cast<size_t>(j1 - j0));

    uint64_t enc_lr = 0, enc_bp = 0;
    uint64_t l = 0;
    const uint64_t i1 = std::min(j1 + k - 1, len);
    for (uint64_t i = j0; i < i1; ++i) {
      if (__builtin_expect(SEQ_NT4_TABLE[static_cast<uint8_t>(cseq[i])] >= 4, 0)) {
        l = 0;
        continue;
      }
      ++l;
      if (l < k) continue;
      const uint64_t j = i - k + 1;
      if (l == k) {
        compute_encoding(cseq + j, cseq + i + 1, enc_lr, enc_bp);
      } else {
        update_encoding(cseq + i, enc_lr, enc_bp);
      }
      enc_bp &= mask_bp;
      enc_lr &= mask_lr;
      if (__builtin_expect(j >= j1, 0)) break;

      // Valid k-mer reached. The window's nvalid counts every valid k-mer
      // (observed or frac-unobserved); dist2 infers unobserved from it and the
      // pool size per window. Only observed k-mers (bix < nrows) enter the pool.
      ++nvalid_fw;
      const uint64_t rc_bp = revcomp_bp64(enc_bp, static_cast<uint8_t>(k));
      if constexpr (StrandAware) {
        ++nvalid_rc;
        const uint32_t bix_fw = lshf.compute_hash_bp(enc_bp);
        if (bix_fw < nrows) out_fw.push_back(pack_win_hash(bix_fw, lshf.drop_ppos_lr(enc_lr)));
        const uint32_t bix_rc = lshf.compute_hash_bp(rc_bp);
        if (bix_rc < nrows) out_rc.push_back(pack_win_hash(bix_rc, lshf.drop_ppos_lr(bp64_to_lr64(rc_bp))));
      } else if (canonical) {
        uint32_t bix;
        enc_t enc;
        if (rc_bp < enc_bp) {
          bix = lshf.compute_hash_bp(enc_bp);
          enc = lshf.drop_ppos_lr(enc_lr);
        } else {
          bix = lshf.compute_hash_bp(rc_bp);
          enc = lshf.drop_ppos_lr(bp64_to_lr64(rc_bp));
        }
        if (bix < nrows) out_fw.push_back(pack_win_hash(bix, enc));
      } else {
        const uint32_t bix = lshf.compute_hash_bp(enc_bp);
        if (bix < nrows) out_fw.push_back(pack_win_hash(bix, lshf.drop_ppos_lr(enc_lr)));
      }
    }
  }

} // namespace

void read_path_list(const std::filesystem::path& list_path, std::vector<str>& paths, std::vector<str>& names)
{
  std::ifstream in(list_path);
  check_fstream(in, "Cannot open input list", list_path.string());
  std::string line;
  std::string name, path;
  while (std::getline(in, line)) {
    const size_t first_nonspace = line.find_first_not_of(" \t\r");
    if (first_nonspace == std::string::npos) continue;
    if (line[first_nonspace] == '#') continue;
    line.erase(0, first_nonspace); // trim leading whitespace
    while (!line.empty() && (line.back() == ' ' || line.back() == '\t' || line.back() == '\r'))
      line.pop_back();
    const size_t tab = line.find('\t');
    if (tab != std::string::npos) {
      name = line.substr(0, tab);
      path = line.substr(tab + 1);
      while (!name.empty() && (name.back() == ' ' || name.back() == '\t'))
        name.pop_back();
      while (!path.empty() && (path.front() == ' ' || path.front() == '\t' || path.front() == '\r'))
        path.erase(0, 1);
    } else {
      name.clear();
      path = line;
    }
    if (path.empty()) continue;
    paths.push_back(path);
    names.push_back(name);
  }
}

namespace {

  void sort_pool_by_lsh(win_hash_pool_t& pool)
  {
    const size_t n = pool.hashes.size();
    if (n <= 1) return;
    vec<size_t> order(n);
    std::iota(order.begin(), order.end(), size_t(0));
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return pool.hashes[a] < pool.hashes[b]; });
    vec<uint64_t> hashes(n);
    vec<uint32_t> win_ix(n);
    for (size_t i = 0; i < n; ++i) {
      hashes[i] = pool.hashes[order[i]];
      win_ix[i] = pool.win_ix[order[i]];
    }
    pool.hashes = std::move(hashes);
    pool.win_ix = std::move(win_ix);
  }

} // namespace

Sketch2::Sketch2(std::filesystem::path sketch_path)
  : sketch_path(std::move(sketch_path))
{
}

const char* Sketch2::need(const char* p, const char* end, size_t n, const str& what)
{
  if (p == nullptr || end == nullptr || static_cast<size_t>(end - p) < n) {
    error_exit("Truncated sketch2 while reading " + what);
  }
  return p;
}

void Sketch2::build_nonempty_bitmap()
{
  nonempty_bits.assign((static_cast<size_t>(nrows) + 63) / 64, 0);
  if (!sfhm || nrows == 0) return;
  sfhm->fill_nonempty_bitmap(nonempty_bits.data(), nrows);
}

const char* Sketch2::skip_sfhm_bytes(const char* p, const char* end)
{
  p = need(p, end, sizeof(uint64_t), "SFHM nkmers");
  uint64_t nkmers = 0;
  std::memcpy(&nkmers, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  const size_t enc_bytes = static_cast<size_t>(nkmers) * sizeof(enc_t);
  p = need(p, end, enc_bytes, "SFHM encodings");
  p += enc_bytes;
  p = need(p, end, sizeof(uint32_t), "SFHM nrows");
  uint32_t sfhm_nrows = 0;
  std::memcpy(&sfhm_nrows, p, sizeof(uint32_t));
  p += sizeof(uint32_t);
  const size_t inc_bytes = static_cast<size_t>(sfhm_nrows) * sizeof(inc_t);
  p = need(p, end, inc_bytes, "SFHM increments");
  return p + inc_bytes;
}

const char* Sketch2::skip_windows_bytes(const char* p, const char* end, uint64_t& tau_out, uint64_t& bin_shift_out)
{
  auto skip_pool = [&](const char* q) {
    q = need(q, end, sizeof(uint64_t), "pool size");
    uint64_t n = 0;
    std::memcpy(&n, q, sizeof(uint64_t));
    q += sizeof(uint64_t);
    const size_t hash_bytes = static_cast<size_t>(n) * sizeof(uint64_t);
    const size_t ix_bytes = static_cast<size_t>(n) * sizeof(uint32_t);
    q = need(q, end, hash_bytes + ix_bytes, "pool payloads");
    return q + hash_bytes + ix_bytes;
  };
  p = need(p, end, 2 * sizeof(uint64_t), "tau/bin_shift");
  std::memcpy(&tau_out, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  std::memcpy(&bin_shift_out, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  p = need(p, end, sizeof(uint64_t), "nwins");
  uint64_t nwins = 0;
  std::memcpy(&nwins, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  for (uint64_t wi = 0; wi < nwins; ++wi) {
    p = need(p, end, sizeof(uint64_t), "qid_len");
    uint64_t qid_len = 0;
    std::memcpy(&qid_len, p, sizeof(uint64_t));
    p += sizeof(uint64_t);
    p = need(p, end, static_cast<size_t>(qid_len) + 2 * sizeof(uint64_t) + 2 * sizeof(uint32_t), "window meta");
    p += static_cast<size_t>(qid_len) + 2 * sizeof(uint64_t) + 2 * sizeof(uint32_t);
  }
  p = skip_pool(p);
  return skip_pool(p);
}

// Zero-copy pool view into a live external buffer (mmap). The caller must keep
// the underlying memory alive for the returned pool's lifetime.
const char* Sketch2::view_pool_mem(const char* p, const char* end, win_hash_pool_t& pool)
{
  p = need(p, end, sizeof(uint64_t), "pool size");
  uint64_t n = 0;
  std::memcpy(&n, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  pool.clear();
  pool.hashes_v = reinterpret_cast<const uint64_t*>(p);
  pool.win_ix_v = reinterpret_cast<const uint32_t*>(p + static_cast<size_t>(n) * sizeof(uint64_t));
  pool.n_ = n;
  const size_t payload = static_cast<size_t>(n) * (sizeof(uint64_t) + sizeof(uint32_t));
  p = need(p, end, payload, "pool payloads");
  return p + payload;
}

const char* Sketch2::read_windows_mem(const char* p, const char* end)
{
  p = need(p, end, 2 * sizeof(uint64_t), "tau/bin_shift");
  std::memcpy(&tau, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  std::memcpy(&bin_shift, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  const uint64_t tau_bin = std::max<uint64_t>(1, (tau + bin_size - 1) >> bin_shift);
  nwinmers = tau_bin << bin_shift;

  p = need(p, end, sizeof(uint64_t), "nwins");
  uint64_t nwins = 0;
  std::memcpy(&nwins, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  wins_v.clear();
  wins_v.reserve(static_cast<size_t>(nwins));
  for (uint64_t wi = 0; wi < nwins; ++wi) {
    win2_t win;
    p = need(p, end, sizeof(uint64_t), "qid_len");
    uint64_t qid_len = 0;
    std::memcpy(&qid_len, p, sizeof(uint64_t));
    p += sizeof(uint64_t);
    p = need(p, end, static_cast<size_t>(qid_len) + 2 * sizeof(uint64_t) + 2 * sizeof(uint32_t), "window meta");
    win.qid.assign(p, p + qid_len);
    p += static_cast<size_t>(qid_len);
    std::memcpy(&win.start, p, sizeof(uint64_t));
    p += sizeof(uint64_t);
    std::memcpy(&win.end, p, sizeof(uint64_t));
    p += sizeof(uint64_t);
    std::memcpy(&win.nvalid_fw, p, sizeof(uint32_t));
    p += sizeof(uint32_t);
    std::memcpy(&win.nvalid_rc, p, sizeof(uint32_t));
    p += sizeof(uint32_t);
    wins_v.push_back(std::move(win));
  }
  p = view_pool_mem(p, end, pool_fw);
  p = view_pool_mem(p, end, pool_rc);
  loaded_windows = true;
  return p;
}

void Sketch2::write_pool(std::ostream& stream, const win_hash_pool_t& pool)
{
  uint64_t n = pool.hashes.size();
  stream.write(reinterpret_cast<const char*>(&n), sizeof(uint64_t));
  if (n) {
    stream.write(reinterpret_cast<const char*>(pool.hashes.data()), static_cast<std::streamsize>(n * sizeof(uint64_t)));
    stream.write(reinterpret_cast<const char*>(pool.win_ix.data()), static_cast<std::streamsize>(n * sizeof(uint32_t)));
  }
}

void Sketch2::load_from_offset(const sketch2_entry_t& e, const sketch2_index_t& idx, Sketch2Part part)
{
  mapped = MappedFile::open(sketch_path);
  const char* base = mapped->begin();
  const char* end = mapped->end();
  if (e.offset + e.len > mapped->len) error_exit("sketch2 record past EOF: " + sketch_path.string());
  const char* p = base + static_cast<ptrdiff_t>(e.offset);

  p = need(p, end, sizeof(uint64_t), "rname length");
  uint64_t len_rname = 0;
  std::memcpy(&len_rname, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  p = need(p, end, static_cast<size_t>(len_rname) + 3 * sizeof(uint64_t) + sizeof(double) + sizeof(uint8_t), "record meta");
  rname.assign(p, p + len_rname);
  p += static_cast<ptrdiff_t>(len_rname);
  std::memcpy(&timestamp, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  std::memcpy(&genome_length_bp, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  std::memcpy(&nvalid_bases, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  std::memcpy(&rho, p, sizeof(double));
  p += sizeof(double);
  uint8_t has_sfhm = 0;
  std::memcpy(&has_sfhm, p, sizeof(uint8_t));
  p += sizeof(uint8_t);

  // LSH geometry + window parameters are file-level (shared params block).
  k = idx.k;
  w = idx.w;
  h = idx.h;
  nrows = idx.nrows;
  canonical = idx.canonical;
  tau = idx.tau;
  bin_shift = idx.bin_shift;
  lshf = std::make_shared<LSHF>(idx.ppos, idx.npos);
  const uint64_t bin_size2 = uint64_t(1) << bin_shift;
  const uint64_t tau_bin2 = std::max<uint64_t>(1, (tau + bin_size2 - 1) >> bin_shift);
  nwinmers = tau_bin2 << bin_shift;

  const bool want_buckets = (static_cast<uint8_t>(part) & static_cast<uint8_t>(Sketch2Part::Buckets)) != 0;
  const bool want_windows = (static_cast<uint8_t>(part) & static_cast<uint8_t>(Sketch2Part::Windows)) != 0;

  if (has_sfhm) {
    if (want_buckets) {
      if (e.buckets_len == 0) error_exit("Corrupt sketch2 record (SFHM without bucket bytes): " + sketch_path.string());
      const char* sp = base + static_cast<ptrdiff_t>(e.buckets_off);
      const char* sp_end = sp + e.buckets_len;
      sfhm = std::make_shared<SFHM>();
      sfhm->view_mem(sp, sp_end);
      build_nonempty_bitmap();
    } else {
      sfhm = nullptr;
      nonempty_bits.clear();
    }
  } else if (want_buckets) {
    error_exit("sketch2 has no SFHM buckets (windows-only record cannot be a reference): " + sketch_path.string());
  } else {
    sfhm = nullptr;
    nonempty_bits.clear();
  }

  if (want_windows) {
    if (e.windows_len == 0) {
      wins_v.clear();
      pool_fw.clear();
      pool_rc.clear();
      loaded_windows = false;
    } else {
      const char* wp = base + static_cast<ptrdiff_t>(e.windows_off);
      const char* wp_end = wp + e.windows_len;
      (void)read_windows_mem(wp, wp_end);
      loaded_windows = true;
    }
  } else {
    uint64_t tau_tmp = 0, bin_tmp = 0;
    if (e.windows_len)
      (void)skip_windows_bytes(base + static_cast<ptrdiff_t>(e.windows_off),
                               base + static_cast<ptrdiff_t>(e.windows_off + e.windows_len),
                               tau_tmp,
                               bin_tmp);
    wins_v.clear();
    pool_fw.clear();
    pool_rc.clear();
    loaded_windows = false;
  }
}

void Sketch2::advise_range(uint64_t off, uint64_t len, int advice) const noexcept
{
  if (!mapped) return;
  if (off >= mapped->len) return;
  const size_t l = static_cast<size_t>(std::min<uint64_t>(len, mapped->len - off));
  if (l == 0) return;
  void* addr = const_cast<char*>(mapped->begin()) + off;
  ::madvise(addr, l, advice);
}

bool Sketch2::scan_bucket(uint32_t bix, enc_t enc_lr, uint32_t& hdist_min) const noexcept
{
  if (bix == INVALID_BIX || bix >= nrows) return false;
  const enc_t* ix1 = sfhm->bucket_ptr_start(bix);
  const enc_t* ix2 = sfhm->bucket_ptr_next(bix);
  uint32_t hdist_curr = std::numeric_limits<uint32_t>::max();
  for (; ix1 < ix2; ++ix1) {
    const uint32_t hdist = popcount_lr32((*ix1) ^ enc_lr);
    hdist_curr = hdist < hdist_curr ? hdist : hdist_curr;
  }
  hdist_min = hdist_curr;
  return true;
}

bool Sketch2::compatible_with(const Sketch2& other) const
{
  if (k != other.k || h != other.h || w != other.w) return false;
  if (nrows != other.nrows || canonical != other.canonical) return false;
  if (tau != other.tau || bin_shift != other.bin_shift) return false;
  if (lshf->get_ppos_v() != other.lshf->get_ppos_v()) return false;
  if (lshf->get_npos_v() != other.lshf->get_npos_v()) return false;
  return true;
}

sketch2_index_t read_sketch2_index(const std::filesystem::path& sketch_path)
{
  auto mapped = Sketch2::MappedFile::open(sketch_path);
  const char* p = mapped->begin();
  const char* end = mapped->end();
  sketch2_index_t idx;
  p = Sketch2::need(p, end, 2 * sizeof(uint32_t) + sizeof(uint64_t), "magic/version/nsketches");
  uint32_t magic = 0, version = 0;
  std::memcpy(&magic, p, sizeof(uint32_t));
  p += sizeof(uint32_t);
  std::memcpy(&version, p, sizeof(uint32_t));
  p += sizeof(uint32_t);
  if (magic != SKETCH2_MAGIC) {
    error_exit("Not a sketch2 file (bad magic / rebuild with sketch2): " + sketch_path.string());
  }
  if (version != SKETCH2_VERSION) {
    error_exit("Unsupported sketch2 format version " + std::to_string(version) + " (expected " +
               std::to_string(SKETCH2_VERSION) + "): " + sketch_path.string());
  }
  std::memcpy(&idx.nsketches, p, sizeof(uint64_t));
  p += sizeof(uint64_t);

  auto need8 = [&](const char* what) {
    p = Sketch2::need(p, end, sizeof(uint64_t), what);
    uint64_t v = 0;
    std::memcpy(&v, p, sizeof(uint64_t));
    p += sizeof(uint64_t);
    return v;
  };
  idx.seed = need8("seed");
  idx.timestamp = need8("timestamp");
  p = Sketch2::need(
    p, end, 3 * sizeof(uint8_t) + sizeof(bool) + sizeof(uint32_t) + sizeof(double) + 3 * sizeof(uint64_t), "shared params");
  std::memcpy(&idx.k, p, sizeof(uint8_t));
  p += sizeof(uint8_t);
  std::memcpy(&idx.w, p, sizeof(uint8_t));
  p += sizeof(uint8_t);
  std::memcpy(&idx.h, p, sizeof(uint8_t));
  p += sizeof(uint8_t);
  std::memcpy(&idx.canonical, p, sizeof(bool));
  p += sizeof(bool);
  std::memcpy(&idx.nrows, p, sizeof(uint32_t));
  p += sizeof(uint32_t);
  std::memcpy(&idx.frac, p, sizeof(double));
  p += sizeof(double);
  std::memcpy(&idx.tau, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  std::memcpy(&idx.bin_shift, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  std::memcpy(&idx.sample_size, p, sizeof(uint64_t));
  p += sizeof(uint64_t);

  if (idx.h > idx.k) error_exit("Corrupt sketch2 LSH config in " + sketch_path.string());
  p = Sketch2::need(p, end, static_cast<size_t>(idx.k), "ppos/npos");
  idx.ppos.assign(p, p + idx.h);
  p += idx.h;
  idx.npos.assign(p, p + (idx.k - idx.h));
  p += (idx.k - idx.h);

  idx.records.resize(static_cast<size_t>(idx.nsketches));
  for (auto& r : idx.records) {
    p = Sketch2::need(p, end, sizeof(sketch2_entry_t), "index entry");
    std::memcpy(&r, p, sizeof(sketch2_entry_t));
    p += sizeof(sketch2_entry_t);
  }
  return idx;
}

bool Sketch2SC::validate_configuration()
{
  bool is_invalid = false;
  if (frac <= 0.0 || frac > 1.0) {
    is_invalid = true;
    cerr_msg("--frac must be in (0, 1]; got ", frac);
  }
  if (w < k) {
    is_invalid = true;
    cerr_msg("The minimum minimizer window size (-w) is k (-k)!");
  }
  if (h < 3) {
    is_invalid = true;
    cerr_msg("The minimum number of LSH positions (-h) is 3!");
  }
  if (h > 15) {
    is_invalid = true;
    cerr_msg("The maximum number of LSH positions (-h) is 15!");
  }
  if (k > 31) {
    is_invalid = true;
    cerr_msg("The maximum allowed k-mer length (-k) is 31!");
  }
  if (k < 19) {
    is_invalid = true;
    cerr_msg("The minimum allowed k-mer length (-k) is 19!");
  }
  if ((k - h) > 16) {
    is_invalid = true;
    cerr_msg("For compact k-mer encodings, h must be >= k-16!");
  }
  if (tau == 0) {
    is_invalid = true;
    cerr_msg("-l (window length in k-mers) is required and must be positive");
  }
  if (bin_shift > 16) {
    is_invalid = true;
    cerr_msg("--bin-shift must be <= 16; got ", bin_shift);
  }
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  if (bin_size > tau) {
    is_invalid = true;
    cerr_msg("--bin-shift gives bin_size=", bin_size, ", which exceeds -l=", tau);
  }
  return !is_invalid;
}

win_sample_t Sketch2SC::sample_windows(const str& input_path, uint64_t& genome_length_bp, uint64_t& nvalid_bases)
{
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  const uint64_t tau_bin = std::max<uint64_t>(1, (tau + bin_size - 1) >> bin_shift);
  const uint64_t nwinmers = tau_bin << bin_shift;
  const uint64_t xtau = nwinmers + k - 1;
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  const uint64_t mask_bp = u64m >> ((32 - k) * 2);
  const uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));

  qseq_sptr_t qs = std::make_shared<QSeq>(input_path);
  while (qs->read_next_batch()) {
  }
  const vec<qseq_t>& batch_v = qs->get_batch_v();

  // Per-genome metadata: total sequence length (bp) and the count of ACGT bases
  // (N-ambiguous positions excluded), gathered while the batch is already in
  // memory. Written into the record preamble and surfaced by dist2.
  genome_length_bp = 0;
  nvalid_bases = 0;
  for (const qseq_t& bq : batch_v) {
    genome_length_bp += bq.seq.size();
    for (char bc : bq.seq) {
      if (SEQ_NT4_TABLE[static_cast<uint8_t>(bc)] < 4) ++nvalid_bases;
    }
  }

  uint64_t total_npos = 0;
  vec<uint64_t> lenc_v(batch_v.size());
  for (size_t bix = 0; bix < batch_v.size(); ++bix) {
    const uint64_t L = batch_v[bix].seq.size();
    if (L >= xtau) {
      const uint64_t enmers = L - k + 1;
      total_npos += (enmers - nwinmers) / bin_size + 1;
    }
    lenc_v[bix] = total_npos;
  }

  win_sample_t sample;
  if (total_npos == 0) {
    warn_msg(concat_msg("No eligible windows in ", input_path, " for -l=", tau));
    return sample;
  }

  vec<uint64_t> positions_v = sample_coords(total_npos, sample_size, gen);
  std::sort(positions_v.begin(), positions_v.end());
  sample.wins_v.reserve(positions_v.size());

  size_t pidx = 0;
  vec<uint64_t> tmp_fw, tmp_rc;
  for (size_t bix = 0; bix < batch_v.size() && pidx < positions_v.size(); ++bix) {
    const uint64_t L = batch_v[bix].seq.size();
    if (L < xtau) continue;
    const uint64_t rend = lenc_v[bix];
    const uint64_t roff = bix == 0 ? 0 : lenc_v[bix - 1];
    const uint64_t enmers = L - k + 1;
    const char* cseq = batch_v[bix].seq.data();

    while (pidx < positions_v.size() && positions_v[pidx] < rend) {
      const uint64_t start_bin = positions_v[pidx] - roff;
      const uint64_t jx = start_bin << bin_shift;
      const uint64_t jy = std::min(jx + nwinmers, enmers);
      const uint32_t wix = static_cast<uint32_t>(sample.wins_v.size());

      win2_t win;
      win.qid = batch_v[bix].qid;
      win.start = jx;
      win.end = jy;

      if (canonical) {
        uint32_t nvalid_rc = 0;
        collect_range_hashes<false>(
          cseq, L, jx, jy, *lshf, nrows, true, mask_bp, mask_lr, tmp_fw, tmp_rc, win.nvalid_fw, nvalid_rc);
      } else {
        collect_range_hashes<true>(
          cseq, L, jx, jy, *lshf, nrows, false, mask_bp, mask_lr, tmp_fw, tmp_rc, win.nvalid_fw, win.nvalid_rc);
      }
      for (uint64_t h : tmp_fw) {
        sample.pool_fw.hashes.push_back(h);
        sample.pool_fw.win_ix.push_back(wix);
      }
      for (uint64_t h : tmp_rc) {
        sample.pool_rc.hashes.push_back(h);
        sample.pool_rc.win_ix.push_back(wix);
      }

      sample.wins_v.push_back(std::move(win));
      ++pidx;
    }
  }

  // Global LSH order so dist2 walks reference buckets consecutively.
  sort_pool_by_lsh(sample.pool_fw);
  sort_pool_by_lsh(sample.pool_rc);
  return sample;
}

void Sketch2SC::write_file_header(std::ostream& sout, uint64_t nsketches)
{
  // GSK5 file header: magic, version, nsketches, then the shared params block.
  // All records in one bundle share one LSH geometry + window sampling config,
  // so readers parse it once instead of per record.
  const uint32_t magic = SKETCH2_MAGIC;
  const uint32_t version = SKETCH2_VERSION;
  sout.write(reinterpret_cast<const char*>(&magic), sizeof(uint32_t));
  sout.write(reinterpret_cast<const char*>(&version), sizeof(uint32_t));
  sout.write(reinterpret_cast<const char*>(&nsketches), sizeof(uint64_t));
  const uint64_t seed64 = seed;
  const uint64_t timestamp =
    std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  sout.write(reinterpret_cast<const char*>(&seed64), sizeof(uint64_t));
  sout.write(reinterpret_cast<const char*>(&timestamp), sizeof(uint64_t));
  sout.write(reinterpret_cast<const char*>(&k), sizeof(uint8_t));
  sout.write(reinterpret_cast<const char*>(&w), sizeof(uint8_t));
  sout.write(reinterpret_cast<const char*>(&h), sizeof(uint8_t));
  sout.write(reinterpret_cast<char*>(&canonical), sizeof(bool));
  sout.write(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));
  sout.write(reinterpret_cast<const char*>(&frac), sizeof(double));
  sout.write(reinterpret_cast<const char*>(&tau), sizeof(uint64_t));
  sout.write(reinterpret_cast<const char*>(&bin_shift), sizeof(uint64_t));
  sout.write(reinterpret_cast<const char*>(&sample_size), sizeof(uint64_t));
  sout.write(reinterpret_cast<const char*>(lshf->ppos_data()), h * sizeof(uint8_t));
  sout.write(reinterpret_cast<const char*>(lshf->npos_data()), (k - h) * sizeof(uint8_t));
}

sketch2_entry_t Sketch2SC::write_record(std::string& bytes,
                                        const str& rname,
                                        uint64_t timestamp,
                                        uint64_t genome_len,
                                        uint64_t nvalid,
                                        double rho,
                                        const sfhm_sptr_t& sketch_sfhm,
                                        const win_sample_t& sample)
{
  // Serialize one record into `bytes` (produce once per genome in parallel)
  // with the entry's section offsets recorded relative to the byte start. The
  // index still carries absolute file offsets, resolved at flush time.
  std::ostringstream os(std::ios::out | std::ios::binary);
  sketch2_entry_t e;
  e.offset = 0;
  e.len = 0;
  const uint64_t len_rname = rname.length();
  os.write(reinterpret_cast<const char*>(&len_rname), sizeof(uint64_t));
  os.write(rname.c_str(), static_cast<std::streamsize>(len_rname));
  os.write(reinterpret_cast<const char*>(&timestamp), sizeof(uint64_t));
  os.write(reinterpret_cast<const char*>(&genome_len), sizeof(uint64_t));
  os.write(reinterpret_cast<const char*>(&nvalid), sizeof(uint64_t));
  os.write(reinterpret_cast<const char*>(&rho), sizeof(double));
  const uint8_t has_sfhm = 1;
  os.write(reinterpret_cast<const char*>(&has_sfhm), sizeof(uint8_t));

  e.buckets_off = static_cast<uint64_t>(os.tellp());
  sketch_sfhm->save(os);
  e.buckets_len = static_cast<uint64_t>(os.tellp()) - e.buckets_off;
  e.windows_off = static_cast<uint64_t>(os.tellp());
  write_windows(os, sample);
  e.windows_len = static_cast<uint64_t>(os.tellp()) - e.windows_off;
  e.len = static_cast<uint64_t>(os.tellp());
  bytes = os.str();
  return e;
}

void Sketch2SC::write_windows(std::ostream& sout, const win_sample_t& sample)
{
  sout.write(reinterpret_cast<const char*>(&tau), sizeof(uint64_t));
  sout.write(reinterpret_cast<const char*>(&bin_shift), sizeof(uint64_t));
  uint64_t nwins = sample.wins_v.size();
  sout.write(reinterpret_cast<const char*>(&nwins), sizeof(uint64_t));
  for (const win2_t& win : sample.wins_v) {
    uint64_t qid_len = win.qid.size();
    sout.write(reinterpret_cast<const char*>(&qid_len), sizeof(uint64_t));
    sout.write(win.qid.data(), static_cast<std::streamsize>(qid_len));
    sout.write(reinterpret_cast<const char*>(&win.start), sizeof(uint64_t));
    sout.write(reinterpret_cast<const char*>(&win.end), sizeof(uint64_t));
    sout.write(reinterpret_cast<const char*>(&win.nvalid_fw), sizeof(uint32_t));
    sout.write(reinterpret_cast<const char*>(&win.nvalid_rc), sizeof(uint32_t));
  }
  Sketch2::write_pool(sout, sample.pool_fw);
  Sketch2::write_pool(sout, sample.pool_rc);
}

void Sketch2SC::process()
{
  // Fold --input-list entries into the input set (kept before the -i paths so
  // record order stays deterministic: list entries first, then CLI paths).
  if (!input_list_path.empty()) {
    std::vector<str> list_paths, list_names;
    read_path_list(input_list_path, list_paths, list_names);
    paths_v.insert(paths_v.begin(), list_paths.begin(), list_paths.end());
    rnames_v.insert(rnames_v.begin(), list_names.begin(), list_names.end());
  }
  if (paths_v.empty()) error_exit("No input files provided! Use -i <path> and/or --input-list <file>.");

  const uint64_t nsketches = paths_v.size();
  const uint32_t nthreads = std::max(1u, num_threads);
  cerr_msg("Preparing to sketch2 ", nsketches, " file(s) (", nthreads, " thread(s), global LSH-sorted windows", ")");

  std::ofstream sketch_stream(sketch_path, std::ofstream::binary);
  write_file_header(sketch_stream, nsketches);

  // Index block placeholder at a fixed location; patched in place after all
  // records are laid down so readers never walk per-record headers.
  const std::streampos index_pos = sketch_stream.tellp();
  vec<sketch2_entry_t> entries(nsketches);
  sketch_stream.write(reinterpret_cast<const char*>(entries.data()),
                      static_cast<std::streamsize>(sizeof(sketch2_entry_t) * nsketches));

  struct pending_t
  {
    std::string bytes;
    sketch2_entry_t rel; // section offsets relative to this record's start
    uint64_t windows = 0;
    uint64_t hashes = 0;
  };
  vec<pending_t> pending(nsketches);

  rho_v.assign(nsketches, 0.0);

  // Produce per-genome records in parallel; each worker seeds the thread-local
  // RNG with (1 + genome index) so stream 0 is identical to the single-genome
  // legacy build and the bundle is fully deterministic across thread counts.
  {
    ThreadPool pool(nthreads);
    pool.parallel_for(nsketches, 1, [&](uint64_t i) {
      init_thread_rng(static_cast<uint32_t>(1 + i));
      const str& input_path = paths_v[i];
      rseq_sptr_t rs = std::make_shared<RSeq>(input_path, lshf, w, nrows, canonical);
      sdhm_sptr_t sdhm = std::make_shared<SDHM>();
      sdhm->fill_table(nrows, rs);
      sfhm_sptr_t sketch_sfhm = std::make_shared<SFHM>(sdhm);
      rho_v[i] = rs->get_rho();

      uint64_t genome_len = 0, nvalid = 0;
      win_sample_t sample = sample_windows(input_path, genome_len, nvalid);

      str rname;
      if (i < rnames_v.size() && !rnames_v[i].empty()) {
        rname = rnames_v[i];
      } else {
        rname = std::filesystem::path(input_path).filename().string();
      }
      const uint64_t timestamp =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      pending[i].rel = write_record(pending[i].bytes, rname, timestamp, genome_len, nvalid, rho_v[i], sketch_sfhm, sample);
      pending[i].windows = sample.wins_v.size();
      pending[i].hashes = sample.pool_fw.size();
    });
  }

  // Commit records in index order (deterministic); resolve relative section
  // offsets to absolute file offsets as the stream position advances.
  uint64_t pos = static_cast<uint64_t>(sketch_stream.tellp());
  for (uint64_t i = 0; i < nsketches; ++i) {
    pending_t& rec = pending[i];
    entries[i] = rec.rel;
    entries[i].offset = pos;
    if (entries[i].buckets_len) entries[i].buckets_off += pos;
    if (entries[i].windows_len) entries[i].windows_off += pos;
    pos += rec.bytes.size();
    sketch_stream.write(rec.bytes.data(), static_cast<std::streamsize>(rec.bytes.size()));

    std::cerr << "\rCreated sketch2 [" << (i + 1) << "/" << nsketches << "] "
              << "(windows=" << rec.windows << ", hashes=" << rec.hashes << ")..." << std::flush;
  }
  std::cerr << std::endl;

  // Patch the index with the laid-down record offsets.
  sketch_stream.seekp(index_pos);
  sketch_stream.write(reinterpret_cast<const char*>(entries.data()),
                      static_cast<std::streamsize>(sizeof(sketch2_entry_t) * nsketches));
  sketch_stream.seekp(0, std::ios::end);

  check_fstream(sketch_stream, std::string("Failed to write the sketch2!"), sketch_path.string());
  sketch_stream.close();
  cerr_msg("Sketch2 file saved to ", sketch_path.string(), " with ", nsketches, " sketch(es)");
}

Sketch2SC::Sketch2SC(CLI::App& sc)
{
  set_sketch_defaults();
  sc.add_option("-i,--input-path", paths_v, "Input FASTA/FASTQ file(s) <path> (or URL) (gzip compatible)")
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("--input-list",
                input_list_path,
                "Read input paths from a file, one per line "
                "(optional record name + TAB + path); combines with -i");
  sc.add_option("-o,--output-path", sketch_path, "Path to store the resulting binary sketch2 file")->required();
  sc.add_option("-k,--mer-len", k, "Length of k-mers [27]")->check(CLI::Range(19, 31))->check(CLI::PositiveNumber);
  sc.add_option("-w,--win-len", w, "Length of the minimizer window (w>=k) [k+6]")->check(CLI::PositiveNumber);
  sc.add_option("-h,--num-positions", h, "Number of positions for the LSH [k-16]")->check(CLI::PositiveNumber);
  sc.add_option("--frac", frac, "Keep a k-mer if LSH(x) < frac * 2^(2h); i.e., subsampling ratio [1.0]")
    ->check(CLI::Range(std::numeric_limits<double>::min(), 1.0));
  sc.add_flag(
    "--strand-agnostic,!--strand-aware", canonical, "A (canonical) strand-agnostic (default) or strand-aware sketch");
  sc.add_option("-l", tau, "Length of sampled windows in k-mers")->required()->check(CLI::PositiveNumber);
  sc.add_option("-b,--bin-shift", bin_shift, "Group consecutive k-mers into bins of size 2^b [0]")->check(CLI::Range(0, 16));
  sc.add_option("--sample-size", sample_size, "Windows sampled across the genome [200]")->check(CLI::PositiveNumber);
  sc.callback([&]() {
    if (!(sc.count("-w") + sc.count("--win-len"))) {
      w = k + 6;
      h = k - 16;
    }
    if (!validate_configuration()) {
      error_exit("Invalid configuration!");
    }
  });
}
