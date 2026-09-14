#include "dist.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
// #include <limits>
#include <numeric>
#include <sstream>
#include <sys/mman.h>
#include <unistd.h>
#include <unordered_set>

#include "common.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "scan.hpp"

extern uint32_t num_threads;

namespace {

  LLH<double> make_llhf(const Sketch& sketch, uint32_t hdist_th)
  {
    return {sketch.get_k(), sketch.get_h(), sketch.get_rho(), hdist_th, 0.0, false};
  }

  vec<uint64_t> sample_random_coordinates(const uint64_t npos, const uint64_t nsamples, std::mt19937& rng)
  {
    assert(npos >= 1);
    const uint64_t n = std::min(npos, nsamples);
    // A full shuffle is O(npos) memory; reject-sample instead while n << npos.
    if (n < npos / 2) {
      vec<uint64_t> starts_v;
      starts_v.reserve(n);
      std::unordered_set<uint64_t> seen;
      seen.reserve(static_cast<size_t>(n) * 2);
      std::uniform_int_distribution<uint64_t> pick(0, npos - 1);
      while (starts_v.size() < n) {
        const uint64_t x = pick(rng);
        if (seen.insert(x).second) starts_v.push_back(x);
      }
      return starts_v;
    }
    vec<uint64_t> starts_v(npos);
    std::iota(starts_v.begin(), starts_v.end(), uint64_t(0));
    std::shuffle(starts_v.begin(), starts_v.end(), rng);
    starts_v.resize(n);
    return starts_v;
  }

  // Per-window match histograms for one strand of one direction.
  struct win_counts_t
  {
    vec<uint64_t> hist; // nwins * (hdist_bound + 1)
    vec<uint64_t> u;    // nwins

    void assign(size_t nwins)
    {
      hist.assign(nwins * (hdist_bound + 1), 0);
      u.assign(nwins, 0);
    }
    const uint64_t* row(size_t wix) const noexcept { return hist.data() + wix * (hdist_bound + 1); }
    uint64_t* row(size_t wix) noexcept { return hist.data() + wix * (hdist_bound + 1); }
  };

  // The pool is grouped by bucket, so each bucket is resolved once per sample.
  void accumulate_pool(const Buckets& buckets,
                       const hash_pool_t& pool,
                       uint32_t hdist_th,
                       size_t nwins,
                       win_counts_t& out) noexcept
  {
    const uint64_t* hashes = pool.hashes_ptr();
    const uint16_t* win_ix = pool.win_ix_ptr();
    const uint64_t n = pool.size();
    for (uint64_t i = 0; i < n;) {
      const uint32_t bix = key_bix(hashes[i]);
      uint64_t j = i + 1;
      while (j < n && key_bix(hashes[j]) == bix)
        ++j;
      const enc_t* beg = nullptr;
      const enc_t* end = nullptr;
      if (!buckets.range(bix, beg, end)) {
        i = j;
        continue;
      }
      __builtin_prefetch(beg, 0, 0);
      for (uint64_t e = i; e < j; ++e) {
        const size_t wix = win_ix[e];
        if (wix >= nwins) continue;
        const uint32_t hd = bucket_hdist_min(beg, end, key_enc(hashes[e]));
        if (hd <= hdist_th) ++out.row(wix)[hd];
      }
      i = j;
    }
  }

  // Per-window counts for one direction; rc stays empty in canonical mode.
  void query_windows(const Sketch& query,
                     const Sketch& reference,
                     uint32_t hdist_th,
                     uint64_t nmers_limit,
                     win_counts_t& fw,
                     win_counts_t& rc)
  {
    const vec<window_t>& wins = query.get_wins();
    const size_t nwins = wins.size();
    const bool canonical = query.is_canonical();
    fw.assign(nwins);
    if (!canonical) rc.assign(nwins);

    if (query.get_config().win_repr == WinRepr::Pool) {
      const Buckets& buckets = reference.get_buckets();
      accumulate_pool(buckets, query.get_windows().pool_fw, hdist_th, nwins, fw);
      if (!canonical) accumulate_pool(buckets, query.get_windows().pool_rc, hdist_th, nwins, rc);
      // A pool holds only retained k-mers; the rest are recovered as misses.
      for (size_t wix = 0; wix < nwins; ++wix) {
        uint64_t matched = 0;
        const uint64_t* row = fw.row(wix);
        for (uint32_t d = 0; d <= hdist_th; ++d)
          matched += row[d];
        fw.u[wix] = wins[wix].nvalid_fw > matched ? wins[wix].nvalid_fw - matched : 0;
        if (!canonical) {
          matched = 0;
          const uint64_t* rrow = rc.row(wix);
          for (uint32_t d = 0; d <= hdist_th; ++d)
            matched += rrow[d];
          rc.u[wix] = wins[wix].nvalid_rc > matched ? wins[wix].nvalid_rc - matched : 0;
        }
      }
      return;
    }

    // scan_mers_range already counts misses, so u comes out directly.
    const scan_ctx_t ctx = make_scan_ctx(reference, 0, hdist_th);
    const seq_pack_t& packs = query.get_windows().packs;
    str cseq;
    for (size_t wix = 0; wix < nwins; ++wix) {
      packs.unpack(wix, cseq);
      const uint64_t nmers = std::min(nmers_limit, wins[wix].end - wins[wix].start);
      if (canonical) {
        window_counts_t agg(hdist_th);
        scan_mers_range<false>(ctx, cseq.data(), 0, nmers, agg);
        std::copy(agg.hist(), agg.hist() + hdist_bound + 1, fw.row(wix));
        fw.u[wix] = agg.u;
      } else {
        swindow_counts_t agg(hdist_th);
        scan_mers_range<true>(ctx, cseq.data(), 0, nmers, agg);
        std::copy(agg.hist_fw(), agg.hist_fw() + hdist_bound + 1, fw.row(wix));
        fw.u[wix] = agg.u_fw;
        std::copy(agg.hist_rc(), agg.hist_rc() + hdist_bound + 1, rc.row(wix));
        rc.u[wix] = agg.u_rc;
      }
    }
  }

  // A pair is complete only once both directions ran, hence the halved counts.
  void print_pair_progress(uint64_t done_jobs, uint64_t total_jobs)
  {
    if (!isatty(STDERR_FILENO)) return; // keep redirected logs clean
    constexpr int bar_width = 40;
    const double frac = total_jobs ? static_cast<double>(done_jobs) / static_cast<double>(total_jobs) : 1.0;
    const int filled = static_cast<int>(frac * bar_width + 0.5);
    std::ostringstream os;
    os << "\rpairs [" << str(static_cast<size_t>(filled), '#') << str(static_cast<size_t>(bar_width - filled), '.') << "] "
       << std::setw(3) << static_cast<int>(frac * 100.0 + 0.5) << "% (" << std::fixed << std::setprecision(1)
       << 0.5 * static_cast<double>(done_jobs) << "/" << 0.5 * static_cast<double>(total_jobs) << ")" << std::flush;
    std::cerr << os.str();
  }

} // namespace

xy_t bracket_distance(const double d, const vec<double>& th_v)
{
  xy_t d_range{d_eps, d_ub};
  if (!is_valid_distance(d)) return d_range;
  const auto it = std::lower_bound(th_v.begin(), th_v.end(), d);
  if (it != th_v.begin()) d_range.first = *(it - 1);
  if (it != th_v.end()) d_range.second = *it;
  return d_range;
}

std::pair<double, char> select_strand_distance(const double d_fw, const double d_rc)
{
  const bool fw_valid = is_valid_distance(d_fw);
  const bool rc_valid = is_valid_distance(d_rc);
  if (!fw_valid && !rc_valid) return {nanx(), '.'};
  if (!rc_valid || (fw_valid && d_fw <= d_rc)) return {d_fw, '+'};
  return {d_rc, '-'};
}

dir_result_t
run_direction(const Sketch& query, const Sketch& reference, uint32_t hdist_th, uint64_t nmers_limit, bool output_samples)
{
  if (!query.get_config().compatible_with(reference.get_config())) {
    error_exit(concat_msg("Incompatible sketch pair: ",
                          query.get_rname(),
                          " vs ",
                          reference.get_rname(),
                          " (k/w/h/nrows/LSH/-l must match; use the same --seed)"));
  }
  if (!reference.has_buckets()) {
    error_exit(concat_msg("Reference sketch has no buckets: ", reference.get_rname()));
  }
  if (!query.has_windows()) {
    error_exit(concat_msg("Query sketch has no sampled windows: ", query.get_rname()));
  }

  // lr_ub depends on the reference's rho, so it differs between directions.
  const LLH<double> llhf = make_llhf(reference, hdist_th);
  const bool canonical = query.is_canonical();
  const vec<window_t>& wins = query.get_wins();
  const size_t nwins = wins.size();

  win_counts_t fw, rc;
  query_windows(query, reference, hdist_th, nmers_limit, fw, rc);

  dir_result_t out;
  out.rows.reserve(nwins);
  vec<double> d_v;
  d_v.reserve(nwins);

  // lr_bg needs this direction's median, so rows are emitted in a second pass.
  struct row_t
  {
    char strand;
    const uint64_t* hist;
    uint64_t u;
  };
  vec<row_t> rows_v;
  if (output_samples) rows_v.reserve(nwins);

  for (size_t wix = 0; wix < nwins; ++wix) {
    const uint64_t* hfw = fw.row(wix);
    uint64_t t_fw = 0;
    for (uint32_t d = 0; d <= hdist_th; ++d)
      t_fw += hfw[d];
    const double d_fw = t_fw == 0 ? nanx() : llhf.mle(hfw, fw.u[wix]);

    double d = d_fw;
    char strand = canonical ? '.' : '+';
    const uint64_t* hist_sel = hfw;
    uint64_t u_sel = fw.u[wix];
    if (!canonical) {
      const uint64_t* hrc = rc.row(wix);
      uint64_t t_rc = 0;
      for (uint32_t di = 0; di <= hdist_th; ++di)
        t_rc += hrc[di];
      const double d_rc = t_rc == 0 ? nanx() : llhf.mle(hrc, rc.u[wix]);
      const auto picked = select_strand_distance(d_fw, d_rc);
      d = picked.first;
      strand = picked.second;
      if (strand == '-') {
        hist_sel = hrc;
        u_sel = rc.u[wix];
      }
    }

    dpoint_t row;
    row.d = d;
    if (is_valid_distance(d)) {
      uint64_t n_total = u_sel;
      for (uint32_t di = 0; di <= hdist_th; ++di)
        n_total += hist_sel[di];
      row.lr_ub = compute_lr_ub(llhf, d, n_total);
      d_v.push_back(d);
    }
    out.rows.push_back(row);
    if (output_samples) rows_v.push_back({strand, hist_sel, u_sel});
  }

  std::sort(d_v.begin(), d_v.end());
  out.d_median = linear_quantile(d_v, 0.5);
  out.n_valid = d_v.size();
  out.n_unmapped = nwins - out.n_valid;

  if (output_samples) {
    strstream ss;
    set_precision(ss, 5);
    // set_precision(ss, std::numeric_limits<double>::max_digits10);
    for (size_t wix = 0; wix < nwins; ++wix) {
      const row_t& row = rows_v[wix];
      double lr_bg = nanx();
      if (is_valid_distance(out.rows[wix].d) && is_valid_distance(out.d_median)) {
        lr_bg =
          likelihood_ratio_statistic(llhf.nll(out.d_median, row.hist, row.u), llhf.nll(out.rows[wix].d, row.hist, row.u));
      }
      write_tsv(ss,
                "gdiff",
                query.get_rname(),
                reference.get_rname(),
                wins[wix].qid,
                wins[wix].start + 1,
                wins[wix].end + query.get_k() - 1,
                row.strand,
                reference.get_rname(),
                out.rows[wix].d,
                lr_bg,
                out.rows[wix].lr_ub)
        << '\n';
    }
    out.samples = ss.str();
  }
  return out;
}

DistanceSampler::DistanceSampler(const Sketch& sketch,
                                 const vec<qseq_t>& batch_v,
                                 uint64_t tau,
                                 uint64_t bin_shift,
                                 uint32_t hdist_th)
  : sketch(sketch)
  , batch_v(batch_v)
  , hdist_th(hdist_th)
  , tau(tau)
  , bin_shift(bin_shift)
  , llhf(make_llhf(sketch, hdist_th))
  , k(llhf.k)
{
  canonical = sketch.is_canonical();
  bin_size = uint64_t(1) << bin_shift;
  tau_bin = std::max<uint64_t>(1, (this->tau + bin_size - 1) >> bin_shift);
  nwinmers = tau_bin << bin_shift;
}

void DistanceSampler::run_for_all(uint64_t sample_size, bool keep_counts, ThreadPool& pool)
{
  build_for_all(sample_size, keep_counts);
  evaluate(pool);
}

void DistanceSampler::run_per_sequence(uint64_t sample_size, bool keep_counts, ThreadPool& pool)
{
  build_per_sequence(sample_size, keep_counts);
  evaluate(pool);
}

void DistanceSampler::build_for_all(uint64_t sample_size, bool keep_counts)
{
  schemes_v.clear();
  const uint64_t xtau = nwinmers + k - 1;
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
  if (total_npos == 0) return;

  vec<uint64_t> positions_v = sample_random_coordinates(total_npos, sample_size, gen);
  std::sort(positions_v.begin(), positions_v.end());

  schemes_v.reserve(batch_v.size());
  size_t pidx = 0;
  for (size_t bix = 0; bix < batch_v.size() && pidx < positions_v.size(); ++bix) {
    const uint64_t L = batch_v[bix].seq.size();
    if (L < xtau) continue;
    const uint64_t rend = lenc_v[bix];
    const uint64_t roff = bix == 0 ? 0 : lenc_v[bix - 1];
    vec<uint64_t> starts_v;
    while (pidx < positions_v.size() && positions_v[pidx] < rend) {
      starts_v.push_back(positions_v[pidx] - roff);
      ++pidx;
    }
    if (starts_v.empty()) continue;
    const uint64_t enmers = L - k + 1;
    const uint64_t nbins = (enmers + bin_size - 1) >> bin_shift;
    schemes_v.emplace_back(bix, starts_v.size(), enmers, nbins, std::move(starts_v), keep_counts);
  }
}

void DistanceSampler::build_per_sequence(const uint64_t sample_size, const bool keep_counts)
{
  schemes_v.clear();
  const uint64_t xtau = nwinmers + k - 1;
  schemes_v.reserve(batch_v.size());
  for (size_t bix = 0; bix < batch_v.size(); ++bix) {
    const uint64_t L = batch_v[bix].seq.size();
    if (L < xtau) continue;
    const uint64_t enmers = L - k + 1;
    const uint64_t nbins = (enmers + bin_size - 1) >> bin_shift;
    const uint64_t npos = (enmers - nwinmers) / bin_size + 1;
    vec<uint64_t> starts_v = sample_random_coordinates(npos, sample_size, gen);
    schemes_v.emplace_back(bix, starts_v.size(), enmers, nbins, std::move(starts_v), keep_counts);
  }
}

void DistanceSampler::evaluate(ThreadPool& pool)
{
  const uint32_t nworkers = pool.size();
  const scan_ctx_t ctx = make_scan_ctx(sketch, bin_shift, hdist_th);

  struct task_t
  {
    uint32_t six;
    uint64_t a, b;
  };
  vec<task_t> tasks_v;
  for (uint32_t six = 0; six < schemes_v.size(); ++six) {
    const scheme_t& scheme = schemes_v[six];
    const uint64_t nper = std::max<uint64_t>(4, (scheme.nsamples + nworkers * 4 - 1) / (nworkers * 4));
    for (uint64_t s0 = 0; s0 < scheme.nsamples; s0 += nper)
      tasks_v.push_back({six, s0, std::min(s0 + nper, scheme.nsamples)});
  }

  pool.parallel_for(tasks_v.size(), 1, [&](const uint64_t ti) {
    const task_t& t = tasks_v[ti];
    scheme_t& scheme = schemes_v[t.six];
    const char* cseq = batch_v[scheme.bix].seq.data();
    if (canonical) {
      window_counts_t agg(hdist_th);
      for (uint64_t s = t.a; s < t.b; ++s) {
        const uint64_t jx = scheme.starts_v[s] << bin_shift;
        const uint64_t jy = std::min(jx + nwinmers, scheme.enmers);
        agg.clear();
        scan_mers_range<false>(ctx, cseq, jx, jy, agg);
        const double d = llhf.mle(agg.hist(), agg.u);
        scheme.d_v[s] = d;
        scheme.strand_v[s] = '.';
        if (scheme.keep_counts && is_valid_distance(d)) {
          std::copy(agg.hist(), agg.hist() + hdist_bound + 1, scheme.hist_v.data() + s * (hdist_bound + 1));
          scheme.u_v[s] = agg.u;
        }
      }
    } else {
      swindow_counts_t agg(hdist_th);
      for (uint64_t s = t.a; s < t.b; ++s) {
        const uint64_t jx = scheme.starts_v[s] << bin_shift;
        const uint64_t jy = std::min(jx + nwinmers, scheme.enmers);
        agg.clear();
        scan_mers_range<true>(ctx, cseq, jx, jy, agg);
        const double d_fw = llhf.mle(agg.hist_fw(), agg.u_fw);
        const double d_rc = llhf.mle(agg.hist_rc(), agg.u_rc);
        const auto [d, strand] = select_strand_distance(d_fw, d_rc);
        scheme.d_v[s] = d;
        scheme.strand_v[s] = strand;
        if (scheme.keep_counts && is_valid_distance(d)) {
          uint64_t* dst = scheme.hist_v.data() + s * (hdist_bound + 1);
          if (strand == '-') {
            std::copy(agg.hist_rc(), agg.hist_rc() + hdist_bound + 1, dst);
            scheme.u_v[s] = agg.u_rc;
          } else {
            std::copy(agg.hist_fw(), agg.hist_fw() + hdist_bound + 1, dst);
            scheme.u_v[s] = agg.u_fw;
          }
        }
      }
    }
  });
}

void DistanceSampler::collect_distances(vec<double>& d_v) const
{
  for (const auto& scheme : schemes_v) {
    for (const double d : scheme.d_v) {
      if (is_valid_distance(d)) d_v.push_back(d);
    }
  }
}

void DistanceSampler::collect_distances(vec<vec<double>>& d_vvec) const
{
  for (const auto& scheme : schemes_v) {
    auto& dst = d_vvec[scheme.bix];
    for (const double d : scheme.d_v) {
      if (is_valid_distance(d)) dst.push_back(d);
    }
  }
}

dpoint_t DistanceSampler::sample_row(const scheme_t& sch, uint64_t s) const
{
  dpoint_t row;
  row.d = sch.d_v[s];
  if (!sch.keep_counts || !is_valid_distance(row.d)) return row;
  const uint64_t* hist = sch.hist_v.data() + s * (hdist_bound + 1);
  uint64_t n_total = sch.u_v[s];
  for (uint32_t d = 0; d <= llhf.hdist_th; ++d)
    n_total += hist[d];
  row.lr_ub = compute_lr_ub(llhf, row.d, n_total);
  return row;
}

void DistanceSampler::collect_samples(vec<dpoint_t>& rows) const
{
  for (const scheme_t& sch : schemes_v) {
    for (uint64_t s = 0; s < sch.nsamples; ++s)
      rows.push_back(sample_row(sch, s));
  }
}

void DistanceSampler::collect_samples(vec<vec<dpoint_t>>& rows_per_seq) const
{
  for (const scheme_t& sch : schemes_v) {
    auto& dst = rows_per_seq[sch.bix];
    for (uint64_t s = 0; s < sch.nsamples; ++s)
      dst.push_back(sample_row(sch, s));
  }
}

uint64_t DistanceSampler::get_nsamples() const
{
  uint64_t n = 0;
  for (const auto& scheme : schemes_v)
    n += scheme.nsamples;
  return n;
}

const SketchFile* DistSC::file_for(const std::filesystem::path& path)
{
  for (const auto& f : files) {
    if (f->get_path() == path) return f.get();
  }
  files.push_back(std::make_unique<SketchFile>(path));
  return files.back().get();
}

void DistSC::resolve_source(const std::filesystem::path& path, vec<uint32_t>& out)
{
  if (is_sketch_file(path)) {
    const SketchFile* file = file_for(path);
    for (uint32_t rec = 0; rec < file->size(); ++rec) {
      auto it =
        std::find_if(members.begin(), members.end(), [&](const member_t& m) { return m.file == file && m.rec == rec; });
      if (it != members.end()) {
        out.push_back(static_cast<uint32_t>(it - members.begin()));
        continue;
      }
      members.push_back({file, rec, str{}});
      out.push_back(static_cast<uint32_t>(members.size() - 1));
    }
    return;
  }

  // A list file of container paths, relative entries against its own directory.
  std::ifstream in(path);
  check_fstream(in, "Cannot open input list", path.string());
  const std::filesystem::path base_dir = path.parent_path();
  str line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    line.erase(line.find_last_not_of(" \t\r\n") + 1);
    if (line.empty()) continue;
    std::filesystem::path entry(line);
    if (entry.is_relative()) entry = base_dir / entry;
    resolve_source(entry, out);
  }
}

// A list's first entry resolves to a container (or another list). depth caps cycles.
static bool is_sketch_list(const std::filesystem::path& path, int depth = 0)
{
  if (depth > 8 || is_sketch_file(path)) return false;
  std::ifstream in(path);
  if (!in) return false;
  str line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line[0] == '>' || line[0] == '@') return false;
    line.erase(line.find_last_not_of(" \t\r\n") + 1);
    if (line.empty()) continue;
    std::filesystem::path entry(line);
    if (entry.is_relative()) entry = path.parent_path() / entry;
    if (!std::filesystem::exists(entry)) return false;
    return is_sketch_file(entry) || is_sketch_list(entry, depth + 1);
  }
  return false;
}

void DistSC::emit_header(std::ostream& os) const
{
  if (no_header) return;
  if (output_samples) {
    write_tsv(os, "config", "genome_a", "genome_b", "qid", "start", "end", "strand", "reference", "d", "lr_bg", "lr_ub")
      << '\n';
  } else if (symmetric) {
    write_tsv(os,
              "genome_a",
              "genome_b",
              "distance",
              "median",
              "num_filtered",
              "alternative_mean",
              "num_NA",
              "max_unfiltered_distance",
              "max_distance",
              "n_lr_zero",
              "n_ab",
              "d_ab",
              "n_ba",
              "d_ba")
      << '\n';
  } else {
    write_tsv(os, "genome_a", "genome_b", "n_ab", "d_ab") << '\n';
  }
}

void DistSC::write_pair_row(std::ostream& os, const pair_t& pr, const str& name_a, const str& name_b)
{
  if (!symmetric) {
    write_tsv(os, name_a, name_b, pr.ab.n_valid, pr.ab.d_median) << '\n';
    return;
  }
  const sym_merge_t rc = sym_merge(pr.ab.rows, pr.ba.rows);
  const sym_est_t est = sym_estimate(rc, lr_th, min_portion);
  write_tsv(os,
            name_a,
            name_b,
            est.distance,
            est.median,
            est.num_filtered,
            est.alternative_mean,
            est.num_na,
            est.max_unfiltered,
            est.max_distance,
            est.n_lr_zero,
            pr.ab.n_valid,
            pr.ab.d_median,
            pr.ba.n_valid,
            pr.ba.d_median)
    << '\n';
}

void DistSC::dist_sketches()
{
  if (set_a.empty()) error_exit("Query set is empty; nothing to compare");

  auto dedupe = [](vec<uint32_t>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  };
  dedupe(set_a);
  dedupe(set_b);

  const bool within = set_b.empty();
  cerr_msg("members=",
           members.size(),
           within ? " within-mode query=" : " cross-mode query=",
           set_a.size(),
           " ref=",
           within ? set_a.size() : set_b.size());

  // -l defaults to the stored window length; a Pool cannot subset a window.
  const sketch_config_t& cfg = members.front().file->get_config();
  if (cfg.tau == 0) error_exit("Sketch has no sampled windows; re-run `gdiff sketch` with -l");
  if (tau == 0) tau = cfg.tau;
  if (tau > cfg.tau || (tau != cfg.tau && cfg.win_repr == WinRepr::Pool)) {
    error_exit(concat_msg("-l ",
                          tau,
                          " is not usable with this sketch (-l ",
                          cfg.tau,
                          cfg.win_repr == WinRepr::Pool ? ", pool windows freeze -l; re-sketch with "
                                                          "--window-repr seq for runtime flexibility)"
                                                        : ", -l must not exceed the sketch's)"));
  }

  // Every member queries every other, so its windows stay resident throughout.
  ThreadPool pool(std::max(1u, num_threads));
  vec<Sketch> windows(members.size());
  pool.parallel_for(members.size(), 1, [&](uint64_t mi) {
    const member_t& m = members[static_cast<size_t>(mi)];
    windows[static_cast<size_t>(mi)] = m.file->open(m.rec, SketchPart::Windows);
  });
  for (size_t mi = 0; mi < members.size(); ++mi)
    members[mi].rname = windows[mi].get_rname();

  vec<pair_t> pairs;
  if (within) {
    pairs.reserve(set_a.size() * (set_a.size() - 1) / 2);
    for (size_t i = 0; i < set_a.size(); ++i)
      for (size_t j = i + 1; j < set_a.size(); ++j)
        pairs.push_back({set_a[i], set_a[j], {}, {}});
  } else {
    pairs.reserve(set_a.size() * set_b.size());
    for (const uint32_t q : set_a)
      for (const uint32_t r : set_b)
        pairs.push_back({q, r, {}, {}});
  }
  if (pairs.empty()) {
    cerr_msg("No pairs to compare");
    return;
  }

  // Reference-major, so each reference's pages are faulted in and dropped once.
  struct job_t
  {
    uint32_t pair_ix;
    uint32_t query_ix;
    uint32_t ref_ix;
    bool ba;
  };
  vec<job_t> jobs;
  jobs.reserve(pairs.size() * (symmetric ? 2 : 1));
  for (uint32_t p = 0; p < pairs.size(); ++p) {
    jobs.push_back({p, pairs[p].a, pairs[p].b, false});
    if (symmetric) jobs.push_back({p, pairs[p].b, pairs[p].a, true});
  }
  std::stable_sort(jobs.begin(), jobs.end(), [](const job_t& x, const job_t& y) { return x.ref_ix < y.ref_ix; });

  const uint64_t total_jobs = jobs.size();
  uint64_t done_jobs = 0;
  size_t batch_start = 0;
  while (batch_start < jobs.size()) {
    const uint32_t ref_ix = jobs[batch_start].ref_ix;
    size_t batch_end = batch_start;
    while (batch_end < jobs.size() && jobs[batch_end].ref_ix == ref_ix)
      ++batch_end;

    const member_t& ref = members[ref_ix];
    const Sketch ref_buckets = ref.file->open(ref.rec, SketchPart::Buckets);
    pool.parallel_for(batch_end - batch_start, 1, [&](uint64_t i) {
      const job_t& job = jobs[batch_start + static_cast<size_t>(i)];
      dir_result_t r = run_direction(windows[job.query_ix], ref_buckets, hdist_th, tau, output_samples);
      pair_t& pr = pairs[job.pair_ix];
      if (job.ba)
        pr.ba = std::move(r);
      else
        pr.ab = std::move(r);
    });

    const record_entry_t& e = ref.file->get_index().records[ref.rec];
    ref.file->advise(e.buckets_off, e.buckets_len, MADV_DONTNEED);

    done_jobs += batch_end - batch_start;
    print_pair_progress(done_jobs, total_jobs);
    batch_start = batch_end;
  }
  if (isatty(STDERR_FILENO)) std::cerr << '\n';

  std::ostream& os = *output_stream;
  emit_header(os);
  if (output_samples) {
    for (const pair_t& pr : pairs) {
      os << pr.ab.samples << pr.ba.samples;
    }
  } else {
    set_precision(os, 8);
    for (const pair_t& pr : pairs)
      write_pair_row(os, pr, members[pr.a].rname, members[pr.b].rname);
  }
  os.flush();
}

void DistSC::dist_fasta()
{
  if (!is_sketch_file(set_b_path)) {
    error_exit("With a FASTA query the second argument must be a sketch container: " + set_b_path.string());
  }
  const SketchFile* file = file_for(set_b_path);
  const sketch_config_t& cfg = file->get_config();
  if (tau == 0) tau = cfg.tau ? cfg.tau : 500;
  if (sample_size == 0) sample_size = cfg.sample_size ? cfg.sample_size : 1000;

  ThreadPool pool(std::max(1u, num_threads));
  init_thread_rng(1);

  // The reverse direction reuses one in-memory sketch of the query.
  Sketch query_sketch;
  if (symmetric) {
    if (cfg.tau == 0) {
      error_exit("Symmetric mode needs stored reference windows; re-sketch with -l or pass --no-symmetric");
    }
    const lshf_sptr_t lshf = std::make_shared<LSHF>(cfg.ppos, cfg.npos);
    built_sketch_t built = build_buckets(set_a_path.string(), lshf, cfg.w, cfg.nrows, cfg.canonical);
    query_sketch =
      Sketch(cfg, std::filesystem::path(set_a_path).filename().string(), std::move(built.buckets), built.nkmers, built.rho);
  }

  std::ostream& os = *output_stream;
  emit_header(os);
  set_precision(os, output_samples ? 5 : 8);

  const uint32_t nsketches = file->size();
  cerr_msg("Processing ", nsketches, " sketch(es) w/ ", pool.size(), " thread(s)...");

  for (uint32_t i = 0; i < nsketches; ++i) {
    const Sketch reference = file->open(i, symmetric ? SketchPart::All : SketchPart::Buckets);

    // Sample the query FASTA batch by batch so it never becomes fully resident.
    dir_result_t ab;
    vec<dpoint_t> fwd_rows;
    vec<double> fwd_d;
    strstream fwd_samples;
    if (output_samples) set_precision(fwd_samples, std::numeric_limits<double>::max_digits10);
    {
      QSeq qs(set_a_path.string());
      bool more = true;
      while (more) {
        qs.clear();
        more = qs.read_next_batch();
        if (qs.get_batch_v().empty()) break;
        DistanceSampler sampler(reference, qs.get_batch_v(), tau, bin_shift, hdist_th);
        sampler.run_for_all(sample_size, true, pool);
        sampler.collect_samples(fwd_rows);
        if (output_samples) {
          const vec<qseq_t>& batch_v = qs.get_batch_v();
          const uint64_t nwinmers = sampler.get_nwinmers();
          const uint32_t k = reference.get_k();
          sampler.for_each_counts(
            [&](uint64_t bix, uint64_t enmers, uint64_t start_bin, double d, char strand, const uint64_t* hist, uint64_t u) {
              double lr_ub = nanx();
              if (hist && is_valid_distance(d)) {
                uint64_t n_total = u;
                for (uint32_t di = 0; di <= hdist_th; ++di)
                  n_total += hist[di];
                lr_ub = compute_lr_ub(sampler.get_llhf(), d, n_total);
              }
              const uint64_t jx = start_bin << bin_shift;
              const uint64_t jy = std::min(jx + nwinmers, enmers);
              write_tsv(fwd_samples,
                        "gdiff",
                        set_a_path.filename().string(),
                        reference.get_rname(),
                        batch_v[bix].qid,
                        jx + 1,
                        jy + k - 1,
                        strand,
                        reference.get_rname(),
                        d,
                        nanx(),
                        lr_ub)
                << '\n';
            });
        }
      }
    }
    for (const dpoint_t& row : fwd_rows)
      if (is_valid_distance(row.d)) fwd_d.push_back(row.d);
    std::sort(fwd_d.begin(), fwd_d.end());
    ab.rows = std::move(fwd_rows);
    ab.d_median = linear_quantile(fwd_d, 0.5);
    ab.n_valid = fwd_d.size();
    ab.n_unmapped = ab.rows.size() - ab.n_valid;
    if (output_samples) ab.samples = fwd_samples.str();

    dir_result_t ba;
    if (symmetric) {
      ba = run_direction(reference, query_sketch, hdist_th, tau, output_samples);
    }

    pair_t pr;
    pr.ab = std::move(ab);
    pr.ba = std::move(ba);
    if (output_samples) {
      os << pr.ab.samples << pr.ba.samples;
    } else {
      write_pair_row(os, pr, set_a_path.filename().string(), reference.get_rname());
    }

    const record_entry_t& e = file->get_index().records[i];
    file->advise(e.offset, e.len, MADV_DONTNEED);
    std::cerr << "\rProcessed sketch " << i + 1 << "/" << nsketches << "..." << std::flush;
  }
  std::cerr << std::endl;
  os.flush();
}

void DistSC::dist()
{
  const bool have_pos_a = !set_a_path.empty();
  const bool have_pos_b = !set_b_path.empty();
  const bool have_list_q = !query_list_path.empty();
  const bool have_list_r = !reference_list_path.empty();

  if (have_pos_a && have_list_q) error_exit("Query set given both positionally and via --query-list; pick one");
  if (have_pos_b && have_list_r) {
    error_exit("Reference set given both positionally and via --reference-list; pick one");
  }

  // A FASTA has no records to enumerate, so only the sketch side is resolved.
  if (have_pos_a && have_pos_b && !is_sketch_file(set_a_path) && !is_sketch_list(set_a_path)) {
    dist_fasta();
    return;
  }

  if (have_list_q) {
    resolve_source(query_list_path, set_a);
  } else if (have_pos_a) {
    resolve_source(set_a_path, set_a);
  } else {
    error_exit("No set given: pass `dist <set>` (all pairs within), `dist <a> <b>` (cross), "
               "or `dist <query.fa> <ref.gdsk>`");
  }
  if (have_list_r) {
    resolve_source(reference_list_path, set_b);
  } else if (have_pos_b) {
    resolve_source(set_b_path, set_b);
  }
  dist_sketches();
}

DistSC::DistSC(CLI::App& sc)
{
  sc.add_option("query-set", set_a_path, "Sketch container, list file, or query FASTA/FASTQ")
    ->expected(0, 1)
    ->check(CLI::ExistingFile);
  sc.add_option(
      "reference-set", set_b_path, "Reference sketch container or list file (omit for all pairs within the query set)")
    ->expected(0, 1)
    ->check(CLI::ExistingFile);
  sc.add_option("--query-list", query_list_path, "Query set as a list file (alternative to the first positional)")
    ->excludes("query-set")
    ->check(CLI::ExistingFile);
  sc.add_option(
      "--reference-list", reference_list_path, "Reference set as a list file (alternative to the second positional)")
    ->excludes("reference-set")
    ->check(CLI::ExistingFile);
  sc.add_option("-l", tau, "Length of sampled windows in k-mers [the sketch's]")->check(CLI::PositiveNumber);
  sc.add_option("--sample-size", sample_size, "Windows sampled from a FASTA query [the sketch's]")
    ->check(CLI::PositiveNumber);
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("--lr-th", lr_th, "Likelihood-ratio cut for the symmetric reconciliation filter [3.841]")
    ->check(CLI::NonNegativeNumber);
  sc.add_option(
      "--min-portion", min_portion, "Report the filtered mean only if this fraction of windows survives --lr-th [0.66]")
    ->check(CLI::Range(0.0, 1.0));
  sc.add_flag("--symmetric,!--no-symmetric",
              symmetric,
              "Reconcile both directions of every pair (default) or report one direction only");
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_flag("--output-samples", output_samples, "Write per-window sample rows instead of per-pair summaries");
  sc.add_flag("--no-header", no_header, "Omit the header line");
  sc.callback([&]() {
    if (!validate_configuration()) error_exit("Invalid configuration!");
    if (!output_path.empty()) {
      output_file.open(output_path);
      check_fstream(output_file, "Cannot open output file", output_path.string());
      output_stream = &output_file;
    }
  });
}

bool validate_binning(const uint64_t bin_shift, const uint64_t tau)
{
  bool is_invalid = false;
  if (bin_shift > 16) {
    is_invalid = true;
    cerr_msg("--bin-shift must be less than or equal to 16; got ", bin_shift);
  }
  const uint64_t bin_size = (bin_shift <= 16) ? (uint64_t(1) << bin_shift) : 0;
  if (tau && bin_size > tau) {
    is_invalid = true;
    cerr_msg("--bin-shift gives bin_size=", bin_size, ", which exceeds -l=", tau);
  }
  return !is_invalid;
}

bool DistSC::validate_configuration()
{
  bool is_invalid = false;
  if (!validate_binning(bin_shift, tau)) is_invalid = true;
  if (hdist_th > hdist_bound) {
    is_invalid = true;
    cerr_msg("--hdist-th must be in [0, ", hdist_bound, "]; got ", hdist_th);
  }
  return !is_invalid;
}
