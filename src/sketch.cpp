#include "sketch.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <numeric>
#include <sstream>

#include "common.hpp"
#include "container.hpp"
#include "distance.hpp"
#include "enc.hpp"
#include "random.hpp"
#include "rqseq.hpp"
#include "serialize.hpp"
#include "tpool.hpp"

namespace {

  class strbuf_t : public std::streambuf
  {
  public:
    explicit strbuf_t(str& out) noexcept
      : out_(out)
    {
    }

  protected:
    std::streamsize xsputn(const char* s, std::streamsize n) override
    {
      out_.append(s, static_cast<size_t>(n));
      return n;
    }
    int overflow(int c) override
    {
      if (c != traits_type::eof()) out_.push_back(static_cast<char>(c));
      return c;
    }
    // pad_to_word() asks the stream for its position; report the bytes written so far.
    pos_type seekoff(off_type off, std::ios_base::seekdir dir, std::ios_base::openmode which) override
    {
      if (!(which & std::ios_base::out)) return pos_type(off_type(-1));
      const off_type base = dir == std::ios_base::beg ? 0 : static_cast<off_type>(out_.size());
      return pos_type(base + off);
    }

  private:
    str& out_;
  };

  // Packed keys for retained k-mers in [j0, j1); nvalid counts all valid k-mers.
  template<bool Canonical>
  void collect_range_keys(const char* cseq,
                          uint64_t len,
                          uint64_t j0,
                          uint64_t j1,
                          const LSHF& lshf,
                          uint32_t nrows,
                          uint64_t mask_bp,
                          uint64_t mask_lr,
                          vec<uint64_t>& out_fw_v,
                          vec<uint64_t>& out_rc_v,
                          uint32_t& nvalid_fw,
                          uint32_t& nvalid_rc)
  {
    const uint32_t k = lshf.get_k();
    out_fw_v.clear();
    out_rc_v.clear();
    nvalid_fw = 0;
    nvalid_rc = 0;
    if (j1 <= j0 || len < k) return;
    out_fw_v.reserve(static_cast<size_t>(j1 - j0));
    if constexpr (!Canonical) out_rc_v.reserve(static_cast<size_t>(j1 - j0));

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
      if constexpr (!Canonical) {
        ++nvalid_rc;
        const uint32_t bix_fw = lshf.compute_hash_bp(enc_bp);
        if (bix_fw < nrows) out_fw_v.push_back(pack_key(bix_fw, lshf.drop_ppos_lr(enc_lr)));
        const uint32_t bix_rc = lshf.compute_hash_bp(rc_bp);
        if (bix_rc < nrows) out_rc_v.push_back(pack_key(bix_rc, lshf.drop_ppos_lr(bp64_to_lr64(rc_bp))));
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
        if (bix < nrows) out_fw_v.push_back(pack_key(bix, enc));
      }
    }
  }

  // Sort by packed key, keeping win_ix_v aligned.
  void sort_pool(hash_pool_t& pool)
  {
    const size_t n = pool.hashes_v.size();
    if (n <= 1) return;
    vec<uint32_t> order_v(n);
    std::iota(order_v.begin(), order_v.end(), uint32_t(0));
    std::sort(order_v.begin(), order_v.end(), [&](uint32_t a, uint32_t b) { return pool.hashes_v[a] < pool.hashes_v[b]; });
    vec<uint64_t> hashes_v(n);
    vec<uint16_t> win_ix_v(n);
    for (size_t i = 0; i < n; ++i) {
      hashes_v[i] = pool.hashes_v[order_v[i]];
      win_ix_v[i] = pool.win_ix_v[order_v[i]];
    }
    pool.hashes_v = std::move(hashes_v);
    pool.win_ix_v = std::move(win_ix_v);
  }

  const char* view_pool(const char* p, const char* end, hash_pool_t& pool)
  {
    const uint64_t n = read_trivial<uint64_t>(p, end, "pool size");
    pool.clear();
    pool.hashes_view = {view_array<uint64_t>(p, end, n, "pool hashes"), n};
    pool.win_ix_view = {view_array<uint16_t>(p, end, n, "pool window indices"), n};
    return p;
  }

  void write_pool(std::ostream& os, const hash_pool_t& pool)
  {
    const uint64_t n = pool.hashes_v.size();
    write_trivial(os, n);
    write_array(os, pool.hashes_v.data(), n);
    write_array(os, pool.win_ix_v.data(), n);
  }

} // namespace

Sketch::Sketch(const sketch_config_t& cfg, str rname, Buckets&& buckets, double card, double rho)
  : rname(std::move(rname))
  , rho(rho)
  , nkmers(buckets.get_nkmers())
  , card(card)
  , cfg(cfg)
  , lshf(std::make_shared<LSHF>(cfg.ppos_v, cfg.npos_v))
  , buckets(std::move(buckets))
{
  loaded_buckets = true;
}

void Sketch::canonicalize()
{
  if (cfg.canonical) return;
  const uint64_t mask_bp = get_lsh_masks(cfg.k).bp;
  const uint32_t nrows = cfg.nrows;
  vec<uint64_t> keys_v;
  keys_v.reserve(static_cast<size_t>(buckets.get_nkmers()));
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
      keys_v.push_back(pack_key(bix_new, lshf->drop_ppos_lr(bp64_to_lr64(can_bp))));
    }
  }
  Buckets rebuilt;
  rebuilt.build(nrows, std::move(keys_v));
  buckets = std::move(rebuilt);
  nkmers = buckets.get_nkmers();
  cfg.canonical = true;
}

Sketch Container::open(uint32_t rec, SketchLoad part) const
{
  if (rec >= entries_v.size()) error_exit("Sketch index out of range: " + path.string());
  const sketch_entry_t& e = entries_v[rec];
  const char* base = map->begin();
  const char* end = map->end();
  if (e.offset + e.len > map->size()) error_exit("Sketch past EOF: " + path.string());

  Sketch gs;
  gs.map = map;
  gs.cfg = cfg;
  gs.lshf = lshf;

  const char* p = base + static_cast<ptrdiff_t>(e.offset);
  const uint64_t rname_len = read_trivial<uint64_t>(p, end, "sketch name length");
  p = need(p, end, round_up_word(rname_len), "sketch name");
  gs.rname.assign(p, p + rname_len);
  p += round_up_word(rname_len);
  gs.timestamp = read_trivial<uint64_t>(p, end, "sketch timestamp");
  gs.ntotal_bp = read_trivial<uint64_t>(p, end, "total bp");
  gs.nvalid_bp = read_trivial<uint64_t>(p, end, "valid bp");
  gs.card = static_cast<double>(read_trivial<uint64_t>(p, end, "card"));
  gs.nkmers = read_trivial<uint64_t>(p, end, "nkmers");
  gs.rho = read_trivial<double>(p, end, "rho");

  const bool want_buckets = (static_cast<uint8_t>(part) & static_cast<uint8_t>(SketchLoad::Buckets)) != 0;
  const bool want_windows = (static_cast<uint8_t>(part) & static_cast<uint8_t>(SketchLoad::Windows)) != 0;

  if (want_buckets) {
    if (e.buckets_len == 0) error_exit("Sketch has no buckets: " + gs.rname);
    const char* bp = base + static_cast<ptrdiff_t>(e.buckets_off);
    gs.buckets.view(bp, bp + e.buckets_len, cfg.nrows);
    gs.loaded_buckets = true;
  }
  if (want_windows && e.windows_len) {
    const char* wp = base + static_cast<ptrdiff_t>(e.windows_off);
    const char* wend = wp + e.windows_len;
    const uint64_t nwins = read_trivial<uint64_t>(wp, wend, "nwins");
    gs.windows.wins_v.resize(static_cast<size_t>(nwins));
    for (window_t& win : gs.windows.wins_v) {
      const uint64_t qid_len = read_trivial<uint64_t>(wp, wend, "window qid length");
      wp = need(wp, wend, round_up_word(qid_len), "window qid");
      win.qid.assign(wp, wp + qid_len);
      wp += round_up_word(qid_len);
      win.start = read_trivial<uint64_t>(wp, wend, "window start");
      win.end = read_trivial<uint64_t>(wp, wend, "window end");
      win.nvalid_fw = read_trivial<uint32_t>(wp, wend, "window nvalid_fw");
      win.nvalid_rc = read_trivial<uint32_t>(wp, wend, "window nvalid_rc");
    }
    if (!cfg.keep_seq) {
      wp = view_pool(wp, wend, gs.windows.pool_fw);
      wp = view_pool(wp, wend, gs.windows.pool_rc);
    } else {
      seq_pack_t& packs = gs.windows.packs;
      const uint64_t packed_words = read_trivial<uint64_t>(wp, wend, "packed words");
      packs.packed_view = {view_array<uint64_t>(wp, wend, packed_words, "packed bases"), packed_words};
      const uint64_t nmask_words = read_trivial<uint64_t>(wp, wend, "nmask words");
      packs.nmask_view = {view_array<uint64_t>(wp, wend, nmask_words, "N mask"), nmask_words};
      packs.base_off_view = {view_array<uint64_t>(wp, wend, nwins + 1, "window base offsets"), nwins + 1};
    }
    gs.loaded_windows = true;
  }
  return gs;
}

void read_path_list(const std::filesystem::path& list_path, vec<str>& paths_v, vec<str>& names_v)
{
  std::ifstream in(list_path);
  check_fstream(in, "Cannot open input list", list_path.string());
  str line;
  while (std::getline(in, line)) {
    str name;
    std::filesystem::path entry;
    if (!parse_input_entry(line, name, entry)) continue;
    paths_v.push_back(entry.string());
    names_v.push_back(std::move(name));
  }
}

bool compatible_configs(const sketch_config_t& a, const sketch_config_t& b)
{
  if (a.k != b.k || a.w != b.w || a.h != b.h) return false;
  if (a.nrows != b.nrows || a.canonical != b.canonical) return false;
  if (a.tau != b.tau) return false;
  return a.ppos_v == b.ppos_v && a.npos_v == b.npos_v;
}

Sketch build_sketch(const str& input_path, const sketch_config_t& cfg, str rname)
{
  const lshf_sptr_t lshf = std::make_shared<LSHF>(cfg.ppos_v, cfg.npos_v);
  RSeq rs(input_path, *lshf, cfg.w, cfg.nrows, cfg.canonical);
  vec<uint64_t> keys_v;
  while (rs.read_next_seq()) {
    if (!rs.set_curr_seq()) continue;
    // Reserve from expected minimizer yield (one per w-k+2 bases, times frac).
    if (keys_v.capacity() == 0) {
      const double keep = static_cast<double>(cfg.nrows) / static_cast<double>(uint64_t(1) << (2 * lshf->get_h()));
      const uint64_t denom = std::max<uint64_t>(1, static_cast<uint64_t>(cfg.w) - lshf->get_k() + 2);
      keys_v.reserve(static_cast<size_t>(2.0 * static_cast<double>(rs.get_len()) / static_cast<double>(denom) * keep) + 16);
    }
    rs.extract_mers(keys_v);
  }
  const double card = rs.get_card();
  Buckets buckets;
  buckets.build(cfg.nrows, std::move(keys_v));
  const uint64_t nkmers = buckets.get_nkmers();

  // rho = exact retained count / HLL cardinality; clamped to (0, 1].
  const double raw = card > 0.0 ? static_cast<double>(nkmers) / card : 1.0;
  if (raw > 1.0) {
    warn_msg(concat_msg(
      "rho > 1 for ", input_path, " (retained ", nkmers, " vs estimated ", card, " distinct k-mers); clamping to 1"));
  }
  const double rho = std::min(1.0, std::max(d_eps, raw));
  return Sketch(cfg, std::move(rname), std::move(buckets), card, rho);
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
  if (params.tau == 0) {
    is_invalid = true;
    cerr_msg("-l (window length in k-mers) must be positive");
  }
  if (params.tau < k) {
    is_invalid = true;
    cerr_msg("-l must be at least k (-k); got ", params.tau);
  }
  if (params.sample_size == 0 || params.sample_size > 65535) {
    is_invalid = true;
    cerr_msg("--sample-size must be in [1, 65535]; got ", params.sample_size);
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
  cfg.tau = params.tau;
  cfg.sample_size = params.sample_size;
  cfg.keep_seq = params.keep_seq;
  cfg.ppos_v.assign(lshf->get_ppos_data(), lshf->get_ppos_data() + h);
  cfg.npos_v.assign(lshf->get_npos_data(), lshf->get_npos_data() + (k - h));
  return cfg;
}

window_sample_t SketchSC::sample_windows(const str& input_path, uint64_t& ntotal_bp, uint64_t& nvalid_bp)
{
  const lsh_masks_t masks = get_lsh_masks(k);

  // Pass 1: sequence lengths only.
  vec<uint64_t> len_v;
  ntotal_bp = 0;
  nvalid_bp = 0;
  {
    QSeq qs(input_path);
    bool more = true;
    while (more) {
      qs.clear();
      more = qs.read_next_batch();
      for (const qseq_t& q : qs.get_batch_v()) {
        len_v.push_back(q.seq.size());
        ntotal_bp += q.seq.size();
        for (const char c : q.seq)
          if (SEQ_NT4_TABLE[static_cast<uint8_t>(c)] < 4) ++nvalid_bp;
      }
    }
  }

  // Which windows to keep, spread over every sequence of this input.
  const window_plan_t wp = make_window_plan(len_v, k, params.tau, 0, params.sample_size, gen);
  const uint64_t nwins = wp.nwins();

  window_sample_t sample;
  if (nwins == 0) {
    warn_msg(concat_msg("No eligible windows in ", input_path, " for -l=", params.tau));
    return sample;
  }

  sample.wins_v.reserve(nwins);
  if (params.keep_seq) {
    // Exact payload size for full-length windows: reserve once instead of resizing per window.
    const uint64_t total_bases = nwins * (params.tau + k - 1);
    sample.packs.packed_v.reserve(static_cast<size_t>((total_bases + 31) / 32));
    sample.packs.nmask_v.reserve(static_cast<size_t>((total_bases + 63) / 64));
    sample.packs.base_off_v.reserve(nwins + 1);
    sample.packs.base_off_v.push_back(0);
  }

  // Pass 2: re-stream and extract only the sampled windows.
  QSeq qs(input_path);
  vec<uint64_t> tmp_fw_v, tmp_rc_v;
  size_t pidx = 0;
  size_t s_ix = 0;
  bool more = true;
  while (more && pidx < wp.sources_v.size()) {
    qs.clear();
    more = qs.read_next_batch();
    for (const qseq_t& q : qs.get_batch_v()) {
      const size_t scurr_ix = s_ix++;
      while (pidx < wp.sources_v.size() && wp.sources_v[pidx].bix < scurr_ix)
        ++pidx;
      if (pidx >= wp.sources_v.size()) break;
      if (wp.sources_v[pidx].bix != scurr_ix) continue;
      const uint64_t L = q.seq.size();
      const char* cseq = q.seq.data();

      for (const uint64_t jx : wp.sources_v[pidx].starts_v) {
        const uint64_t jy = jx + params.tau;
        const uint16_t wix = static_cast<uint16_t>(sample.wins_v.size());

        window_t win;
        win.qid = q.qid;
        win.start = jx;
        win.end = jy;
        if (canonical) {
          collect_range_keys<true>(
            cseq, L, jx, jy, *lshf, nrows, masks.bp, masks.lr, tmp_fw_v, tmp_rc_v, win.nvalid_fw, win.nvalid_rc);
          win.nvalid_rc = 0;
        } else {
          collect_range_keys<false>(
            cseq, L, jx, jy, *lshf, nrows, masks.bp, masks.lr, tmp_fw_v, tmp_rc_v, win.nvalid_fw, win.nvalid_rc);
        }

        if (!params.keep_seq) {
          for (const uint64_t hv : tmp_fw_v) {
            sample.pool_fw.hashes_v.push_back(hv);
            sample.pool_fw.win_ix_v.push_back(wix);
          }
          for (const uint64_t hv : tmp_rc_v) {
            sample.pool_rc.hashes_v.push_back(hv);
            sample.pool_rc.win_ix_v.push_back(wix);
          }
        } else {
          seq_pack_t& packs = sample.packs;
          const uint64_t nbases = jy - jx + k - 1;
          const uint64_t b0 = packs.base_off_v.back();
          packs.packed_v.resize(static_cast<size_t>((b0 + nbases + 31) / 32), 0);
          packs.nmask_v.resize(static_cast<size_t>((b0 + nbases + 63) / 64), 0);
          for (uint64_t t = 0; t < nbases; ++t) {
            const uint64_t b = b0 + t;
            const uint8_t code = SEQ_NT4_TABLE[static_cast<uint8_t>(cseq[jx + t])];
            if (code >= 4) {
              packs.nmask_v[static_cast<size_t>(b >> 6)] |= uint64_t(1) << (b & 63);
            } else {
              packs.packed_v[static_cast<size_t>(b >> 5)] |= uint64_t(code) << (2 * (b & 31));
            }
          }
          packs.base_off_v.push_back(b0 + nbases);
        }

        sample.wins_v.push_back(std::move(win));
        ++pidx;
      }
    }
  }

  if (!params.keep_seq) {
    sort_pool(sample.pool_fw);
    sort_pool(sample.pool_rc);
  } else if (std::all_of(sample.packs.nmask_v.begin(), sample.packs.nmask_v.end(), [](uint64_t x) { return x == 0; })) {
    // No ambiguous base anywhere: drop the mask instead of storing zeros.
    sample.packs.nmask_v.clear();
  }
  return sample;
}

void SketchSC::write_file_header(std::ostream& os, uint64_t nsketches)
{
  const uint64_t now =
    std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  write_container_header(os, make_config(now), nsketches);
}

void SketchSC::write_windows(std::ostream& os, const window_sample_t& sample)
{
  const uint64_t nwins = sample.wins_v.size();
  write_trivial(os, nwins);
  for (const window_t& win : sample.wins_v) {
    const uint64_t qid_len = win.qid.size();
    write_trivial(os, qid_len);
    os.write(win.qid.data(), static_cast<std::streamsize>(qid_len));
    pad_to_word(os);
    write_trivial(os, win.start);
    write_trivial(os, win.end);
    write_trivial(os, win.nvalid_fw);
    write_trivial(os, win.nvalid_rc);
  }
  if (!params.keep_seq) {
    write_pool(os, sample.pool_fw);
    write_pool(os, sample.pool_rc);
  } else {
    const seq_pack_t& packs = sample.packs;
    write_trivial(os, static_cast<uint64_t>(packs.packed_v.size()));
    write_array(os, packs.packed_v.data(), packs.packed_v.size());
    write_trivial(os, static_cast<uint64_t>(packs.nmask_v.size()));
    write_array(os, packs.nmask_v.data(), packs.nmask_v.size());
    write_array(os, packs.base_off_v.data(), packs.base_off_v.size());
  }
}

sketch_entry_t SketchSC::write_sketch(str& bytes,
                                      uint64_t timestamp,
                                      uint64_t ntotal_bp,
                                      uint64_t nvalid_bp,
                                      const Sketch& built,
                                      const window_sample_t& sample)
{
  bytes.clear();
  strbuf_t buf(bytes);
  std::ostream os(&buf);
  sketch_entry_t e;
  const str& rname = built.get_rname();
  const uint64_t rname_len = rname.size();
  write_trivial(os, rname_len);
  os.write(rname.data(), static_cast<std::streamsize>(rname_len));
  pad_to_word(os);
  write_trivial(os, timestamp);
  write_trivial(os, ntotal_bp);
  write_trivial(os, nvalid_bp);
  write_trivial(os, static_cast<uint64_t>(built.get_card() + 0.5));
  write_trivial(os, built.get_nkmers());
  write_trivial(os, built.get_rho());

  e.buckets_off = static_cast<uint64_t>(os.tellp());
  built.get_buckets().save(os);
  e.buckets_len = static_cast<uint64_t>(os.tellp()) - e.buckets_off;
  if (!sample.wins_v.empty()) {
    e.windows_off = static_cast<uint64_t>(os.tellp());
    write_windows(os, sample);
    e.windows_len = static_cast<uint64_t>(os.tellp()) - e.windows_off;
  }
  pad_to_word(os);
  e.len = static_cast<uint64_t>(os.tellp());
  return e;
}

void SketchSC::process()
{
  if (!input_list_path.empty()) {
    vec<str> list_paths_v, list_names_v;
    read_path_list(input_list_path, list_paths_v, list_names_v);
    paths_v.insert(paths_v.begin(), list_paths_v.begin(), list_paths_v.end());
    rnames_v.insert(rnames_v.begin(), list_names_v.begin(), list_names_v.end());
  }
  if (paths_v.empty()) error_exit("No input files provided! Use -i <path> and/or --input-list <file>.");

  const uint64_t nsketches = paths_v.size();
  const uint32_t nthreads = std::max(1u, num_threads);
  cerr_msg("Sketching ", nsketches, " file(s) w/ ", nthreads, " thread(s), windows=", params.keep_seq ? "seq" : "default");

  std::ofstream sketch_stream(sketch_path, std::ofstream::binary);
  check_fstream(sketch_stream, "Cannot open output container", sketch_path.string());
  write_file_header(sketch_stream, nsketches);

  // Index placeholder, patched after every sketch is written.
  const std::streampos index_pos = sketch_stream.tellp();
  vec<sketch_entry_t> entries_v(nsketches);
  sketch_stream.write(reinterpret_cast<const char*>(entries_v.data()),
                      static_cast<std::streamsize>(sizeof(sketch_entry_t) * nsketches));
  pad_to_word(sketch_stream);

  struct pending_t
  {
    str bytes;
    sketch_entry_t rel; // section offsets relative to the sketch start
  };
  vec<pending_t> pending_v(nsketches);

  std::atomic<uint64_t> ndone{0};
  {
    ThreadPool pool(nthreads);
    pool.parallel_for(nsketches, 1, [&](uint64_t i) {
      // Seed per genome index so the bundle is independent of scheduling.
      init_thread_rng(static_cast<uint32_t>(1 + i));
      const str& input_path = paths_v[i];
      str rname =
        (i < rnames_v.size() && !rnames_v[i].empty()) ? rnames_v[i] : std::filesystem::path(input_path).filename().string();
      const uint64_t timestamp =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      Sketch built = build_sketch(input_path, make_config(timestamp), rname);

      uint64_t ntotal_bp = 0, nvalid_bp = 0;
      window_sample_t sample = sample_windows(input_path, ntotal_bp, nvalid_bp);

      pending_v[i].rel = write_sketch(pending_v[i].bytes, timestamp, ntotal_bp, nvalid_bp, built, sample);
      const uint64_t n = ndone.fetch_add(1, std::memory_order_relaxed) + 1;
      progress("Sketched", n, nsketches);
    });
  }
  progress_done();

  uint64_t pos = static_cast<uint64_t>(sketch_stream.tellp());
  for (uint64_t i = 0; i < nsketches; ++i) {
    pending_t& rec = pending_v[i];
    entries_v[i] = rec.rel;
    entries_v[i].offset = pos;
    entries_v[i].buckets_off += pos;
    if (entries_v[i].windows_len) entries_v[i].windows_off += pos;
    pos += rec.bytes.size();
    sketch_stream.write(rec.bytes.data(), static_cast<std::streamsize>(rec.bytes.size()));
    str().swap(rec.bytes);
  }

  sketch_stream.seekp(index_pos);
  sketch_stream.write(reinterpret_cast<const char*>(entries_v.data()),
                      static_cast<std::streamsize>(sizeof(sketch_entry_t) * nsketches));
  sketch_stream.seekp(0, std::ios::end);
  check_fstream(sketch_stream, "Failed to write the sketch", sketch_path.string());
  sketch_stream.close();
  cerr_msg("Container saved to ", sketch_path.string(), " with ", nsketches, " sketch(es)");
}

SketchSC::SketchSC(CLI::App& sc)
{
  set_default_params();
  sc.add_option("-i,--input-path", paths_v, "Input FASTA/FASTQ file(s) <path> (or URL) (gzip compatible)")
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("--input-list",
                input_list_path,
                "Read input paths from a file, one per line "
                "(optional sketch name + TAB + path); combines with -i")
    ->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", sketch_path, "Path to store the resulting binary container")->required();
  sc.add_option("-k,--mer-len", k, "Length of k-mers [27]")->check(CLI::Range(19, 31));
  sc.add_option("-w,--win-len", w, "Length of the minimizer window (w>=k) [k+6]")->check(CLI::PositiveNumber);
  sc.add_option("-h,--num-positions", h, "Number of positions for the LSH [max(floor(k/2)-2, k-16)]")
    ->check(CLI::PositiveNumber);
  sc.add_option("--frac", frac, "Keep a k-mer if LSH(x) < frac * 2^(2h); i.e., subsampling ratio [1.0]")
    ->check(CLI::Range(std::numeric_limits<double>::min(), 1.0));
  sc.add_flag(
    "--strand-agnostic,!--strand-aware", canonical, "A (canonical) strand-agnostic (default) or strand-aware sketch");
  sc.add_option("-l", params.tau, "Length of sampled windows in k-mers [500]")->check(CLI::PositiveNumber);
  sc.add_option("--sample-size", params.sample_size, "Windows sampled across each genome [1000]")
    ->check(CLI::Range(1, 65535));
  sc.add_flag("--keep-seq,!--no-keep-seq",
              params.keep_seq,
              "Store sampled windows as 2-bit packed bases (compact) instead of the default pre-resolved keys");
  sc.callback([&]() {
    // -w and -h are independent: each falls back to its own k-derived default.
    if (!(sc.count("-w") + sc.count("--win-len"))) {
      w = k + 6;
    }
    if (!(sc.count("-h") + sc.count("--num-positions"))) {
      // floor(k/2) - 2, but never below k - 16: drop_ppos_lr packs 2(k - h) bits into enc_t.
      h = std::max<uint8_t>(static_cast<uint8_t>(k / 2 - 2), static_cast<uint8_t>(k > 16 ? k - 16 : 0));
    }
    if (!validate_configuration()) {
      error_exit("Invalid configuration!");
    }
  });
}
