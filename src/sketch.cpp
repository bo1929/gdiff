#include "sketch.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <fcntl.h>
#include <map>
#include <numeric>
#include <sstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_set>

#include "common.hpp"
#include "enc.hpp"
#include "random.hpp"
#include "rqseq.hpp"
#include "tpool.hpp"

extern uint32_t num_threads;

class MappedFile
{
public:
  MappedFile() = default;
  MappedFile(const MappedFile&) = delete;
  MappedFile& operator=(const MappedFile&) = delete;

  ~MappedFile()
  {
    if (ptr && ptr != MAP_FAILED && len) munmap(ptr, len);
    if (fd >= 0) close(fd);
  }

  [[nodiscard]] const char* begin() const noexcept { return static_cast<const char*>(ptr); }
  [[nodiscard]] const char* end() const noexcept { return begin() + len; }
  [[nodiscard]] size_t size() const noexcept { return len; }

  void advise(uint64_t off, uint64_t nbytes, int advice) const noexcept
  {
    if (off >= len) return;
    const size_t n = static_cast<size_t>(std::min<uint64_t>(nbytes, len - off));
    if (n == 0) return;
    ::madvise(const_cast<char*>(begin()) + off, n, advice);
  }

  static std::shared_ptr<const MappedFile> open(const std::filesystem::path& path)
  {
    auto m = std::make_shared<MappedFile>();
    m->fd = ::open(path.c_str(), O_RDONLY);
    if (m->fd < 0) error_exit("Cannot open sketch file: " + path.string());
    struct stat st
    {};
    if (fstat(m->fd, &st) != 0) error_exit("Cannot fstat sketch file: " + path.string());
    m->len = static_cast<size_t>(st.st_size);
    if (m->len == 0) error_exit("Empty sketch file: " + path.string());
    m->ptr = mmap(nullptr, m->len, PROT_READ, MAP_PRIVATE, m->fd, 0);
    if (m->ptr == MAP_FAILED) error_exit("mmap failed for sketch file: " + path.string());
    return m;
  }

private:
  int fd = -1;
  void* ptr = nullptr;
  size_t len = 0;
};

namespace {

  constexpr uint64_t round_up_8(uint64_t x) noexcept { return (x + 7) & ~uint64_t(7); }

  const char* need(const char* p, const char* end, size_t n, const char* what)
  {
    if (p == nullptr || end == nullptr || static_cast<size_t>(end - p) < n) {
      error_exit(concat_msg("Truncated sketch while reading ", what));
    }
    return p;
  }

  template<typename T>
  T read_pod(const char*& p, const char* end, const char* what)
  {
    p = need(p, end, sizeof(T), what);
    T v{};
    std::memcpy(&v, p, sizeof(T));
    p += sizeof(T);
    return v;
  }

  template<typename T>
  void write_pod(std::ostream& os, const T& v)
  {
    os.write(reinterpret_cast<const char*>(&v), sizeof(T));
  }

  // View a length-prefixed array of T; sections are padded to 8 bytes.
  template<typename T>
  const T* view_array(const char*& p, const char* end, uint64_t n, const char* what)
  {
    const uint64_t bytes = round_up_8(n * sizeof(T));
    p = need(p, end, bytes, what);
    const T* out = n ? reinterpret_cast<const T*>(p) : nullptr;
    p += bytes;
    return out;
  }

  template<typename T>
  void write_array(std::ostream& os, const T* data, uint64_t n)
  {
    if (n) os.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(n * sizeof(T)));
    pad_to_8(os);
  }

  vec<uint64_t> sample_coords(uint64_t npos, uint64_t nsamples, std::mt19937& rng)
  {
    const uint64_t n = std::min(npos, nsamples);
    // Sparse draw: rejection-sample. Dense draw: shuffle.
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

  // Packed keys for retained k-mers in [j0, j1); nvalid counts all valid k-mers.
  template<bool StrandAware>
  void collect_range_keys(const char* cseq,
                          uint64_t len,
                          uint64_t j0,
                          uint64_t j1,
                          const LSHF& lshf,
                          uint32_t nrows,
                          uint64_t mask_bp,
                          uint64_t mask_lr,
                          vec<uint64_t>& out_fw,
                          vec<uint64_t>& out_rc,
                          uint32_t& nvalid_fw,
                          uint32_t& nvalid_rc)
  {
    const uint32_t k = lshf.get_k();
    out_fw.clear();
    out_rc.clear();
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
      // l >= k implies i + 1 >= k; see the note in scan.hpp.
      const uint64_t j = i + 1 - k;
      if (l == k) {
        compute_encoding(cseq + j, cseq + i + 1, enc_lr, enc_bp);
      } else {
        update_encoding(cseq + i, enc_lr, enc_bp);
      }
      enc_bp &= mask_bp;
      enc_lr &= mask_lr;
      if (__builtin_expect(j >= j1, 0)) break;

      ++nvalid_fw;
      const uint64_t rc_bp = revcomp_bp64(enc_bp, static_cast<uint8_t>(k));
      if constexpr (StrandAware) {
        ++nvalid_rc;
        const uint32_t bix_fw = lshf.compute_hash_bp(enc_bp);
        if (bix_fw < nrows) out_fw.push_back(pack_key(bix_fw, lshf.drop_ppos_lr(enc_lr)));
        const uint32_t bix_rc = lshf.compute_hash_bp(rc_bp);
        if (bix_rc < nrows) out_rc.push_back(pack_key(bix_rc, lshf.drop_ppos_lr(bp64_to_lr64(rc_bp))));
      } else {
        uint32_t bix;
        enc_t enc;
        if (rc_bp < enc_bp) {
          bix = lshf.compute_hash_bp(enc_bp);
          enc = lshf.drop_ppos_lr(enc_lr);
        } else {
          bix = lshf.compute_hash_bp(rc_bp);
          enc = lshf.drop_ppos_lr(bp64_to_lr64(rc_bp));
        }
        if (bix < nrows) out_fw.push_back(pack_key(bix, enc));
      }
    }
  }

  // Sort by packed key, keeping win_ix aligned.
  void sort_pool(hash_pool_t& pool)
  {
    const size_t n = pool.hashes.size();
    if (n <= 1) return;
    vec<uint32_t> order(n);
    std::iota(order.begin(), order.end(), uint32_t(0));
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) { return pool.hashes[a] < pool.hashes[b]; });
    vec<uint64_t> hashes(n);
    vec<uint16_t> win_ix(n);
    for (size_t i = 0; i < n; ++i) {
      hashes[i] = pool.hashes[order[i]];
      win_ix[i] = pool.win_ix[order[i]];
    }
    pool.hashes = std::move(hashes);
    pool.win_ix = std::move(win_ix);
  }

  const char* view_pool(const char* p, const char* end, hash_pool_t& pool)
  {
    const uint64_t n = read_pod<uint64_t>(p, end, "pool size");
    pool.clear();
    pool.hashes_view = view_array<uint64_t>(p, end, n, "pool hashes");
    pool.win_ix_view = view_array<uint16_t>(p, end, n, "pool window indices");
    pool.n_view = n;
    return p;
  }

  void write_pool(std::ostream& os, const hash_pool_t& pool)
  {
    const uint64_t n = pool.hashes.size();
    write_pod(os, n);
    write_array(os, pool.hashes.data(), n);
    write_array(os, pool.win_ix.data(), n);
  }

} // namespace

bool sketch_config_t::compatible_with(const sketch_config_t& other) const
{
  if (k != other.k || w != other.w || h != other.h) return false;
  if (nrows != other.nrows || canonical != other.canonical) return false;
  if (tau != other.tau) return false;
  return ppos == other.ppos && npos == other.npos;
}

void seq_pack_t::unpack(uint64_t wi, str& cseq) const
{
  static const char bases[4] = {'A', 'C', 'G', 'T'};
  const uint64_t* off = base_off_ptr();
  const uint64_t* packed_p = packed_ptr();
  const uint64_t* nmask_p = nmask_ptr();
  const uint64_t b0 = off[wi];
  const uint64_t b1 = off[wi + 1];
  cseq.resize(static_cast<size_t>(b1 - b0));
  for (uint64_t b = b0; b < b1; ++b) {
    if (nmask_p && ((nmask_p[b >> 6] >> (b & 63)) & 1ull)) {
      cseq[static_cast<size_t>(b - b0)] = 'N';
      continue;
    }
    const uint64_t code = (packed_p[b >> 5] >> (2 * (b & 31))) & 3ull;
    cseq[static_cast<size_t>(b - b0)] = bases[code];
  }
}

Sketch::Sketch(const sketch_config_t& cfg, str rname, Buckets&& buckets, uint64_t nkmers, double rho)
  : rname(std::move(rname))
  , rho(rho)
  , nkmers(nkmers)
  , cfg(cfg)
  , lshf(std::make_shared<LSHF>(cfg.ppos, cfg.npos))
  , buckets(std::move(buckets))
{
}

void Sketch::canonicalize()
{
  if (cfg.canonical) return;
  const uint64_t mask_bp = std::numeric_limits<uint64_t>::max() >> ((32 - cfg.k) * 2);
  const uint32_t nrows = cfg.nrows;
  vec<uint64_t> keys;
  keys.reserve(static_cast<size_t>(buckets.get_nkmers()));
  for (uint32_t bix = 0; bix < nrows; ++bix) {
    const enc_t* beg = nullptr;
    const enc_t* end = nullptr;
    if (!buckets.range(bix, beg, end)) continue;
    const uint64_t bp_ppos = lshf->inv_ppos_bp(bix);
    for (; beg < end; ++beg) {
      const uint64_t bp_npos = lr64_to_bp64(lshf->inv_ppos_lr(*beg));
      const uint64_t fw_bp = (bp_ppos | bp_npos) & mask_bp;
      const uint64_t can_bp = std::max(fw_bp, revcomp_bp64(fw_bp, cfg.k));
      const uint32_t bix_new = lshf->compute_hash_bp(can_bp);
      if (bix_new >= nrows) continue;
      keys.push_back(pack_key(bix_new, lshf->drop_ppos_lr(bp64_to_lr64(can_bp))));
    }
  }
  Buckets rebuilt;
  rebuilt.build(nrows, std::move(keys));
  buckets = std::move(rebuilt);
  nkmers = buckets.get_nkmers();
  cfg.canonical = true;
}

SketchFile::SketchFile(std::filesystem::path path)
  : path(std::move(path))
{
  mapped = MappedFile::open(this->path);
  const char* p = mapped->begin();
  const char* end = mapped->end();

  const uint32_t magic = read_pod<uint32_t>(p, end, "magic");
  const uint32_t version = read_pod<uint32_t>(p, end, "version");
  if (magic != GDSK_MAGIC) error_exit("Not a gdiff sketch file: " + this->path.string());
  if (version != GDSK_VERSION) {
    error_exit(
      concat_msg("Sketch format version ", version, " != ", GDSK_VERSION, "; re-run `gdiff sketch`: ", this->path.string()));
  }
  const uint64_t nsketches = read_pod<uint64_t>(p, end, "nsketches");

  sketch_config_t& cfg = idx.cfg;
  cfg.seed = read_pod<uint64_t>(p, end, "seed");
  cfg.timestamp = read_pod<uint64_t>(p, end, "timestamp");
  cfg.k = read_pod<uint8_t>(p, end, "k");
  cfg.w = read_pod<uint8_t>(p, end, "w");
  cfg.h = read_pod<uint8_t>(p, end, "h");
  cfg.canonical = read_pod<uint8_t>(p, end, "canonical") != 0;
  cfg.win_repr = static_cast<WinRepr>(read_pod<uint8_t>(p, end, "win_repr"));
  p = need(p, end, 3, "config padding");
  p += 3;
  cfg.nrows = read_pod<uint32_t>(p, end, "nrows");
  p = need(p, end, 4, "config padding");
  p += 4;
  cfg.frac = read_pod<double>(p, end, "frac");
  cfg.tau = read_pod<uint64_t>(p, end, "tau");
  cfg.sample_size = read_pod<uint64_t>(p, end, "sample_size");

  if (cfg.h > cfg.k || cfg.k == 0) error_exit("Corrupt sketch LSH config: " + this->path.string());
  p = need(p, end, cfg.k, "LSH positions");
  cfg.ppos.assign(p, p + cfg.h);
  p += cfg.h;
  cfg.npos.assign(p, p + (cfg.k - cfg.h));
  p += (cfg.k - cfg.h);
  p = need(p, end, round_up_8(cfg.k) - cfg.k, "LSH position padding");
  p += round_up_8(cfg.k) - cfg.k;

  idx.records.resize(static_cast<size_t>(nsketches));
  for (record_entry_t& r : idx.records) {
    p = need(p, end, sizeof(record_entry_t), "index entry");
    std::memcpy(&r, p, sizeof(record_entry_t));
    p += sizeof(record_entry_t);
  }
}

void SketchFile::advise(uint64_t off, uint64_t len, int advice) const noexcept
{
  if (mapped) mapped->advise(off, len, advice);
}

Sketch SketchFile::open(uint32_t rec, SketchPart part) const
{
  if (rec >= idx.records.size()) error_exit("Sketch record index out of range: " + path.string());
  const record_entry_t& e = idx.records[rec];
  const char* base = mapped->begin();
  const char* end = mapped->end();
  if (e.offset + e.len > mapped->size()) error_exit("Sketch record past EOF: " + path.string());

  Sketch sk;
  sk.mapped = mapped;
  sk.cfg = idx.cfg;
  sk.lshf = std::make_shared<LSHF>(idx.cfg.ppos, idx.cfg.npos);

  const char* p = base + static_cast<ptrdiff_t>(e.offset);
  const uint64_t rname_len = read_pod<uint64_t>(p, end, "record name length");
  p = need(p, end, round_up_8(rname_len), "record name");
  sk.rname.assign(p, p + rname_len);
  p += round_up_8(rname_len);
  sk.timestamp = read_pod<uint64_t>(p, end, "record timestamp");
  sk.genome_bp = read_pod<uint64_t>(p, end, "genome length");
  sk.nvalid_bases = read_pod<uint64_t>(p, end, "valid base count");
  sk.card_est = static_cast<double>(read_pod<uint64_t>(p, end, "card_est"));
  sk.nkmers = read_pod<uint64_t>(p, end, "nkmers");
  sk.rho = read_pod<double>(p, end, "rho");

  const bool want_buckets = (static_cast<uint8_t>(part) & static_cast<uint8_t>(SketchPart::Buckets)) != 0;
  const bool want_windows = (static_cast<uint8_t>(part) & static_cast<uint8_t>(SketchPart::Windows)) != 0;

  if (want_buckets) {
    if (e.buckets_len == 0) error_exit("Sketch record has no buckets: " + sk.rname);
    const char* bp = base + static_cast<ptrdiff_t>(e.buckets_off);
    sk.buckets.view(bp, bp + e.buckets_len, idx.cfg.nrows);
  }
  if (want_windows && e.windows_len) {
    const char* wp = base + static_cast<ptrdiff_t>(e.windows_off);
    const char* wend = wp + e.windows_len;
    const uint64_t nwins = read_pod<uint64_t>(wp, wend, "nwins");
    sk.windows.wins.resize(static_cast<size_t>(nwins));
    for (window_t& win : sk.windows.wins) {
      const uint64_t qid_len = read_pod<uint64_t>(wp, wend, "window qid length");
      wp = need(wp, wend, round_up_8(qid_len), "window qid");
      win.qid.assign(wp, wp + qid_len);
      wp += round_up_8(qid_len);
      win.start = read_pod<uint64_t>(wp, wend, "window start");
      win.end = read_pod<uint64_t>(wp, wend, "window end");
      win.nvalid_fw = read_pod<uint32_t>(wp, wend, "window nvalid_fw");
      win.nvalid_rc = read_pod<uint32_t>(wp, wend, "window nvalid_rc");
    }
    if (idx.cfg.win_repr == WinRepr::Pool) {
      wp = view_pool(wp, wend, sk.windows.pool_fw);
      wp = view_pool(wp, wend, sk.windows.pool_rc);
    } else {
      seq_pack_t& packs = sk.windows.packs;
      const uint64_t packed_words = read_pod<uint64_t>(wp, wend, "packed words");
      packs.packed_view = view_array<uint64_t>(wp, wend, packed_words, "packed bases");
      packs.packed_words_view = packed_words;
      const uint64_t nmask_words = read_pod<uint64_t>(wp, wend, "nmask words");
      packs.nmask_view = view_array<uint64_t>(wp, wend, nmask_words, "N mask");
      packs.nmask_words_view = nmask_words;
      packs.base_off_view = view_array<uint64_t>(wp, wend, nwins + 1, "window base offsets");
    }
    sk.loaded_windows = true;
  }
  return sk;
}

bool is_sketch_file(const std::filesystem::path& path)
{
  std::ifstream in(path, std::ios::binary);
  if (!in) return false;
  uint32_t magic = 0;
  in.read(reinterpret_cast<char*>(&magic), sizeof(uint32_t));
  return in.good() && magic == GDSK_MAGIC;
}

void read_path_list(const std::filesystem::path& list_path, vec<str>& paths, vec<str>& names)
{
  std::ifstream in(list_path);
  check_fstream(in, "Cannot open input list", list_path.string());
  str line;
  while (std::getline(in, line)) {
    const size_t first = line.find_first_not_of(" \t\r");
    if (first == str::npos || line[first] == '#') continue;
    line.erase(0, first);
    while (!line.empty() && (line.back() == ' ' || line.back() == '\t' || line.back() == '\r'))
      line.pop_back();
    str name, path;
    const size_t tab = line.find('\t');
    if (tab != str::npos) {
      name = line.substr(0, tab);
      path = line.substr(tab + 1);
      while (!name.empty() && (name.back() == ' ' || name.back() == '\t'))
        name.pop_back();
      while (!path.empty() && (path.front() == ' ' || path.front() == '\t'))
        path.erase(0, 1);
    } else {
      path = line;
    }
    if (path.empty()) continue;
    paths.push_back(std::move(path));
    names.push_back(std::move(name));
  }
}

built_sketch_t build_buckets(const str& input_path, const lshf_sptr_t& lshf, uint8_t w, uint32_t nrows, bool canonical)
{
  built_sketch_t out;
  RSeq rs(input_path, lshf, w, nrows, canonical);
  vec<uint64_t> keys;
  while (rs.read_next_seq()) {
    if (!rs.set_curr_seq()) continue;
    // Reserve from expected minimizer yield (one per w-k+2 bases, times frac).
    if (keys.capacity() == 0) {
      const double keep = static_cast<double>(nrows) / static_cast<double>(uint64_t(1) << (2 * lshf->get_h()));
      const uint64_t denom = std::max<uint64_t>(1, static_cast<uint64_t>(w) - lshf->get_k() + 2);
      keys.reserve(static_cast<size_t>(2.0 * static_cast<double>(rs.get_len()) / static_cast<double>(denom) * keep) + 16);
    }
    rs.extract_mers(keys);
  }
  out.card_est = rs.get_card_est();
  out.buckets.build(nrows, std::move(keys));
  out.nkmers = out.buckets.get_nkmers();

  // rho = exact retained count / HLL card_est; clamped to (0, 1].
  const double raw = out.card_est > 0.0 ? static_cast<double>(out.nkmers) / out.card_est : 1.0;
  if (raw > 1.0) {
    warn_msg(concat_msg("rho > 1 for ",
                        input_path,
                        " (retained ",
                        out.nkmers,
                        " vs estimated ",
                        out.card_est,
                        " distinct k-mers); clamping to 1"));
  }
  out.rho = std::min(1.0, std::max(d_eps, raw));
  return out;
}

void BaseLSH::set_lshf() { lshf = std::make_shared<LSHF>(k, h); }

void BaseLSH::set_nrows()
{
  // Flat sampling: keep the prefix [0, T) of the H = 2^(2h) LSH space.
  const uint64_t H = uint64_t(1) << (2 * h);
  const uint64_t t = static_cast<uint64_t>(H * frac + 0.5);
  nrows = static_cast<uint32_t>(std::max<uint64_t>(1, std::min(t, H)));
}

bool SketchSC::validate_configuration()
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
    cerr_msg("-l (window length in k-mers) must be positive");
  }
  if (tau < k) {
    is_invalid = true;
    cerr_msg("-l must be at least k (-k); got ", tau);
  }
  if (sample_size == 0 || sample_size > 65535) {
    is_invalid = true;
    cerr_msg("--sample-size must be in [1, 65535]; got ", sample_size);
  }
  return !is_invalid;
}

sketch_config_t SketchSC::make_config(uint64_t timestamp) const
{
  sketch_config_t cfg;
  cfg.seed = seed;
  cfg.timestamp = timestamp;
  cfg.k = k;
  cfg.w = w;
  cfg.h = h;
  cfg.canonical = canonical;
  cfg.nrows = nrows;
  cfg.frac = frac;
  cfg.tau = tau;
  cfg.sample_size = sample_size;
  cfg.win_repr = win_repr;
  cfg.ppos.assign(lshf->ppos_data(), lshf->ppos_data() + h);
  cfg.npos.assign(lshf->npos_data(), lshf->npos_data() + (k - h));
  return cfg;
}

window_sample_t SketchSC::sample_windows(const str& input_path, uint64_t& genome_bp, uint64_t& nvalid_bases)
{
  const uint64_t xtau = tau + k - 1;
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  const uint64_t mask_bp = u64m >> ((32 - k) * 2);
  const uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));

  // Pass 1: sequence lengths only.
  vec<uint64_t> lens;
  genome_bp = 0;
  nvalid_bases = 0;
  {
    QSeq qs(input_path);
    bool more = true;
    while (more) {
      qs.clear();
      more = qs.read_next_batch();
      for (const qseq_t& q : qs.get_batch_v()) {
        lens.push_back(q.seq.size());
        genome_bp += q.seq.size();
        for (const char c : q.seq)
          if (SEQ_NT4_TABLE[static_cast<uint8_t>(c)] < 4) ++nvalid_bases;
      }
    }
  }

  uint64_t total_npos = 0;
  vec<uint64_t> lenc_v(lens.size());
  for (size_t i = 0; i < lens.size(); ++i) {
    if (lens[i] >= xtau) total_npos += lens[i] - k + 1 - tau + 1;
    lenc_v[i] = total_npos;
  }

  window_sample_t sample;
  if (total_npos == 0) {
    warn_msg(concat_msg("No eligible windows in ", input_path, " for -l=", tau));
    return sample;
  }

  vec<uint64_t> positions_v = sample_coords(total_npos, sample_size, gen);
  std::sort(positions_v.begin(), positions_v.end());
  sample.wins.reserve(positions_v.size());
  if (win_repr == WinRepr::Seq) sample.packs.base_off.push_back(0);

  // Pass 2: re-stream and extract only the sampled windows.
  QSeq qs(input_path);
  vec<uint64_t> tmp_fw, tmp_rc;
  str window_seq;
  size_t pidx = 0;
  size_t seq_ix = 0;
  bool more = true;
  while (more && pidx < positions_v.size()) {
    qs.clear();
    more = qs.read_next_batch();
    for (const qseq_t& q : qs.get_batch_v()) {
      const size_t si = seq_ix++;
      if (pidx >= positions_v.size()) continue;
      const uint64_t L = q.seq.size();
      if (L < xtau) continue;
      const uint64_t rend = lenc_v[si];
      const uint64_t roff = si == 0 ? 0 : lenc_v[si - 1];
      const char* cseq = q.seq.data();

      while (pidx < positions_v.size() && positions_v[pidx] < rend) {
        const uint64_t jx = positions_v[pidx] - roff;
        const uint64_t jy = jx + tau;
        const uint16_t wix = static_cast<uint16_t>(sample.wins.size());

        window_t win;
        win.qid = q.qid;
        win.start = jx;
        win.end = jy;
        if (canonical) {
          collect_range_keys<false>(
            cseq, L, jx, jy, *lshf, nrows, mask_bp, mask_lr, tmp_fw, tmp_rc, win.nvalid_fw, win.nvalid_rc);
          win.nvalid_rc = 0;
        } else {
          collect_range_keys<true>(
            cseq, L, jx, jy, *lshf, nrows, mask_bp, mask_lr, tmp_fw, tmp_rc, win.nvalid_fw, win.nvalid_rc);
        }

        if (win_repr == WinRepr::Pool) {
          for (const uint64_t hv : tmp_fw) {
            sample.pool_fw.hashes.push_back(hv);
            sample.pool_fw.win_ix.push_back(wix);
          }
          for (const uint64_t hv : tmp_rc) {
            sample.pool_rc.hashes.push_back(hv);
            sample.pool_rc.win_ix.push_back(wix);
          }
        } else {
          seq_pack_t& packs = sample.packs;
          const uint64_t nbases = jy - jx + k - 1;
          const uint64_t b0 = packs.base_off.back();
          packs.packed.resize(static_cast<size_t>((b0 + nbases + 31) / 32), 0);
          packs.nmask.resize(static_cast<size_t>((b0 + nbases + 63) / 64), 0);
          for (uint64_t t = 0; t < nbases; ++t) {
            const uint64_t b = b0 + t;
            const uint8_t code = SEQ_NT4_TABLE[static_cast<uint8_t>(cseq[jx + t])];
            if (code >= 4) {
              packs.nmask[static_cast<size_t>(b >> 6)] |= uint64_t(1) << (b & 63);
            } else {
              packs.packed[static_cast<size_t>(b >> 5)] |= uint64_t(code) << (2 * (b & 31));
            }
          }
          packs.base_off.push_back(b0 + nbases);
        }

        sample.wins.push_back(std::move(win));
        ++pidx;
      }
    }
  }

  if (win_repr == WinRepr::Pool) {
    sort_pool(sample.pool_fw);
    sort_pool(sample.pool_rc);
  } else if (std::all_of(sample.packs.nmask.begin(), sample.packs.nmask.end(), [](uint64_t x) { return x == 0; })) {
    // No ambiguous base anywhere: drop the mask instead of storing zeros.
    sample.packs.nmask.clear();
  }
  return sample;
}

void write_sketch_header(std::ostream& os, const sketch_config_t& cfg, uint64_t nsketches)
{
  write_pod(os, GDSK_MAGIC);
  write_pod(os, GDSK_VERSION);
  write_pod(os, nsketches);
  write_pod(os, cfg.seed);
  write_pod(os, cfg.timestamp);
  write_pod(os, cfg.k);
  write_pod(os, cfg.w);
  write_pod(os, cfg.h);
  write_pod(os, static_cast<uint8_t>(cfg.canonical ? 1 : 0));
  write_pod(os, static_cast<uint8_t>(cfg.win_repr));
  const char pad[4] = {};
  os.write(pad, 3);
  write_pod(os, cfg.nrows);
  os.write(pad, 4);
  write_pod(os, cfg.frac);
  write_pod(os, cfg.tau);
  write_pod(os, cfg.sample_size);
  os.write(reinterpret_cast<const char*>(cfg.ppos.data()), static_cast<std::streamsize>(cfg.ppos.size()));
  os.write(reinterpret_cast<const char*>(cfg.npos.data()), static_cast<std::streamsize>(cfg.npos.size()));
  pad_to_8(os);
}

void SketchSC::write_file_header(std::ostream& os, uint64_t nsketches)
{
  const uint64_t now =
    std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  write_sketch_header(os, make_config(now), nsketches);
}

void SketchSC::write_windows(std::ostream& os, const window_sample_t& sample)
{
  const uint64_t nwins = sample.wins.size();
  write_pod(os, nwins);
  for (const window_t& win : sample.wins) {
    const uint64_t qid_len = win.qid.size();
    write_pod(os, qid_len);
    os.write(win.qid.data(), static_cast<std::streamsize>(qid_len));
    pad_to_8(os);
    write_pod(os, win.start);
    write_pod(os, win.end);
    write_pod(os, win.nvalid_fw);
    write_pod(os, win.nvalid_rc);
  }
  if (win_repr == WinRepr::Pool) {
    write_pool(os, sample.pool_fw);
    write_pool(os, sample.pool_rc);
  } else {
    const seq_pack_t& packs = sample.packs;
    write_pod(os, static_cast<uint64_t>(packs.packed.size()));
    write_array(os, packs.packed.data(), packs.packed.size());
    write_pod(os, static_cast<uint64_t>(packs.nmask.size()));
    write_array(os, packs.nmask.data(), packs.nmask.size());
    write_array(os, packs.base_off.data(), packs.base_off.size());
  }
}

record_entry_t SketchSC::write_record(str& bytes,
                                      const str& rname,
                                      uint64_t timestamp,
                                      uint64_t genome_bp,
                                      uint64_t nvalid_bases,
                                      uint64_t nkmers,
                                      double card_est,
                                      double rho,
                                      const Buckets& buckets,
                                      const window_sample_t& sample)
{
  std::ostringstream os(std::ios::out | std::ios::binary);
  record_entry_t e;
  const uint64_t rname_len = rname.size();
  write_pod(os, rname_len);
  os.write(rname.data(), static_cast<std::streamsize>(rname_len));
  pad_to_8(os);
  write_pod(os, timestamp);
  write_pod(os, genome_bp);
  write_pod(os, nvalid_bases);
  write_pod(os, static_cast<uint64_t>(card_est + 0.5));
  write_pod(os, nkmers);
  write_pod(os, rho);

  e.buckets_off = static_cast<uint64_t>(os.tellp());
  buckets.save(os);
  e.buckets_len = static_cast<uint64_t>(os.tellp()) - e.buckets_off;
  if (!sample.wins.empty()) {
    e.windows_off = static_cast<uint64_t>(os.tellp());
    write_windows(os, sample);
    e.windows_len = static_cast<uint64_t>(os.tellp()) - e.windows_off;
  }
  pad_to_8(os);
  e.len = static_cast<uint64_t>(os.tellp());
  bytes = os.str();
  return e;
}

void SketchSC::process()
{
  if (!input_list_path.empty()) {
    vec<str> list_paths, list_names;
    read_path_list(input_list_path, list_paths, list_names);
    paths_v.insert(paths_v.begin(), list_paths.begin(), list_paths.end());
    rnames_v.insert(rnames_v.begin(), list_names.begin(), list_names.end());
  }
  if (paths_v.empty()) error_exit("No input files provided! Use -i <path> and/or --input-list <file>.");

  const uint64_t nsketches = paths_v.size();
  const uint32_t nthreads = std::max(1u, num_threads);
  cerr_msg(
    "Sketching ", nsketches, " file(s) w/ ", nthreads, " thread(s), windows=", win_repr == WinRepr::Pool ? "pool" : "seq");

  std::ofstream sketch_stream(sketch_path, std::ofstream::binary);
  check_fstream(sketch_stream, "Cannot open output sketch file", sketch_path.string());
  write_file_header(sketch_stream, nsketches);

  // Index placeholder, patched after every record is written.
  const std::streampos index_pos = sketch_stream.tellp();
  vec<record_entry_t> entries(nsketches);
  sketch_stream.write(reinterpret_cast<const char*>(entries.data()),
                      static_cast<std::streamsize>(sizeof(record_entry_t) * nsketches));
  pad_to_8(sketch_stream);

  struct pending_t
  {
    str bytes;
    record_entry_t rel; // section offsets relative to the record start
    uint64_t nwins = 0;
  };
  vec<pending_t> pending(nsketches);

  std::atomic<uint64_t> ndone{0};
  {
    ThreadPool pool(nthreads);
    pool.parallel_for(nsketches, 1, [&](uint64_t i) {
      // Seed per genome index so the bundle is independent of scheduling.
      init_thread_rng(static_cast<uint32_t>(1 + i));
      const str& input_path = paths_v[i];
      built_sketch_t built = build_buckets(input_path, lshf, w, nrows, canonical);

      uint64_t genome_bp = 0, nvalid_bases = 0;
      window_sample_t sample = sample_windows(input_path, genome_bp, nvalid_bases);

      str rname =
        (i < rnames_v.size() && !rnames_v[i].empty()) ? rnames_v[i] : std::filesystem::path(input_path).filename().string();
      const uint64_t timestamp =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      pending[i].rel = write_record(pending[i].bytes,
                                    rname,
                                    timestamp,
                                    genome_bp,
                                    nvalid_bases,
                                    built.nkmers,
                                    built.card_est,
                                    built.rho,
                                    built.buckets,
                                    sample);
      pending[i].nwins = sample.wins.size();
      const uint64_t n = ndone.fetch_add(1, std::memory_order_relaxed) + 1;
      std::cerr << "\rSketched [" << n << "/" << nsketches << "]..." << std::flush;
    });
  }
  std::cerr << std::endl;

  uint64_t pos = static_cast<uint64_t>(sketch_stream.tellp());
  for (uint64_t i = 0; i < nsketches; ++i) {
    pending_t& rec = pending[i];
    entries[i] = rec.rel;
    entries[i].offset = pos;
    entries[i].buckets_off += pos;
    if (entries[i].windows_len) entries[i].windows_off += pos;
    pos += rec.bytes.size();
    sketch_stream.write(rec.bytes.data(), static_cast<std::streamsize>(rec.bytes.size()));
    str().swap(rec.bytes);
  }

  sketch_stream.seekp(index_pos);
  sketch_stream.write(reinterpret_cast<const char*>(entries.data()),
                      static_cast<std::streamsize>(sizeof(record_entry_t) * nsketches));
  sketch_stream.seekp(0, std::ios::end);
  check_fstream(sketch_stream, "Failed to write the sketch", sketch_path.string());
  sketch_stream.close();
  cerr_msg("Sketch file saved to ", sketch_path.string(), " with ", nsketches, " record(s)");
}

SketchSC::SketchSC(CLI::App& sc)
{
  set_sketch_defaults();
  sc.add_option("-i,--input-path", paths_v, "Input FASTA/FASTQ file(s) <path> (or URL) (gzip compatible)")
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("--input-list",
                input_list_path,
                "Read input paths from a file, one per line "
                "(optional record name + TAB + path); combines with -i")
    ->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", sketch_path, "Path to store the resulting binary sketch file")->required();
  sc.add_option("-k,--mer-len", k, "Length of k-mers [27]")->check(CLI::Range(19, 31));
  sc.add_option("-w,--win-len", w, "Length of the minimizer window (w>=k) [k+6]")->check(CLI::PositiveNumber);
  sc.add_option("-h,--num-positions", h, "Number of positions for the LSH [k-16]")->check(CLI::PositiveNumber);
  sc.add_option("--frac", frac, "Keep a k-mer if LSH(x) < frac * 2^(2h); i.e., subsampling ratio [1.0]")
    ->check(CLI::Range(std::numeric_limits<double>::min(), 1.0));
  sc.add_flag(
    "--strand-agnostic,!--strand-aware", canonical, "A (canonical) strand-agnostic (default) or strand-aware sketch");
  sc.add_option("-l", tau, "Length of sampled windows in k-mers [500]")->check(CLI::PositiveNumber);
  sc.add_option("--sample-size", sample_size, "Windows sampled across each genome [1000]")->check(CLI::Range(1, 65535));
  sc.add_option("--window-repr", win_repr, "How to store sampled windows: pool (fast) or seq (compact) [pool]")
    ->transform(
      CLI::CheckedTransformer(std::map<str, WinRepr>{{"pool", WinRepr::Pool}, {"seq", WinRepr::Seq}}, CLI::ignore_case));
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
