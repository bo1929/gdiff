#include "estimator.hpp"

#include <algorithm>
#include <sstream>

#include "msg.hpp"
#include "random.hpp"
#include "scan.hpp"
#include "tsv.hpp"
#include "windows.hpp"

namespace {

  // One strand's match histograms for every window of a direction, flattened.
  struct strand_counts_t
  {
    vec<uint64_t> hist_v; // nwins * (hdist_bound + 1)
    vec<uint64_t> u_v;    // nwins

    void assign(size_t nwins)
    {
      hist_v.assign(nwins * (hdist_bound + 1), 0);
      u_v.assign(nwins, 0);
    }
    const uint64_t* row(size_t wix) const noexcept { return hist_v.data() + wix * (hdist_bound + 1); }
    uint64_t* row(size_t wix) noexcept { return hist_v.data() + wix * (hdist_bound + 1); }
  };

  // The pool is grouped by bucket, so each bucket is resolved once per sample.
  void accumulate_pool(const Buckets& buckets,
                       const hash_pool_t& pool,
                       uint32_t hdist_th,
                       size_t nwins,
                       strand_counts_t& out) noexcept
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
                     strand_counts_t& fw,
                     strand_counts_t& rc)
  {
    const vec<window_t>& wins_v = query.get_windows().wins_v;
    const size_t nwins = wins_v.size();
    const bool canonical = query.is_canonical();
    fw.assign(nwins);
    if (!canonical) rc.assign(nwins);

    if (!query.get_config().keep_seq) {
      const Buckets& buckets = reference.get_buckets();
      accumulate_pool(buckets, query.get_windows().pool_fw, hdist_th, nwins, fw);
      if (!canonical) accumulate_pool(buckets, query.get_windows().pool_rc, hdist_th, nwins, rc);
      // A pool holds only retained k-mers; the rest are recovered as misses.
      for (size_t wix = 0; wix < nwins; ++wix) {
        const uint64_t matched = hist_total(fw.row(wix), hdist_th);
        fw.u_v[wix] = wins_v[wix].nvalid_fw > matched ? wins_v[wix].nvalid_fw - matched : 0;
        if (!canonical) {
          const uint64_t rmatched = hist_total(rc.row(wix), hdist_th);
          rc.u_v[wix] = wins_v[wix].nvalid_rc > rmatched ? wins_v[wix].nvalid_rc - rmatched : 0;
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
      const uint64_t nmers = std::min(nmers_limit, wins_v[wix].end - wins_v[wix].start);
      if (canonical) {
        window_counts_t agg(hdist_th);
        scan_mers_range<true>(ctx, cseq.data(), 0, nmers, agg);
        std::copy(agg.hist(), agg.hist() + hdist_bound + 1, fw.row(wix));
        fw.u_v[wix] = agg.u;
      } else {
        swindow_counts_t agg(hdist_th);
        scan_mers_range<false>(ctx, cseq.data(), 0, nmers, agg);
        std::copy(agg.hist_fw(), agg.hist_fw() + hdist_bound + 1, fw.row(wix));
        fw.u_v[wix] = agg.u_fw;
        std::copy(agg.hist_rc(), agg.hist_rc() + hdist_bound + 1, rc.row(wix));
        rc.u_v[wix] = agg.u_rc;
      }
    }
  }

} // namespace

direction_t run_direction(const Sketch& query, const Sketch& reference, uint32_t hdist_th, bool output_samples, bool is_ba)
{
  if (!compatible_configs(query.get_config(), reference.get_config())) {
    error_exit(concat_msg("Incompatible sketch pair: ",
                          query.get_rname(),
                          " vs ",
                          reference.get_rname(),
                          " (k/w/h/nrows/LSH/window length must match; use the same --seed)"));
  }
  if (!reference.has_buckets()) {
    error_exit(concat_msg("Reference sketch has no buckets: ", reference.get_rname()));
  }
  if (!query.has_windows()) {
    error_exit(concat_msg("Query sketch has no sampled windows: ", query.get_rname()));
  }

  // s depends on the reference's rho, so it differs between directions.
  const LLH<double> llhf = make_llhf(reference, hdist_th);
  const bool canonical = query.is_canonical();
  const vec<window_t>& wins_v = query.get_windows().wins_v;
  const size_t nwins = wins_v.size();
  // Compatible sketches share a window length, so the query's is authoritative.
  const uint64_t nmers_limit = query.get_config().tau;

  strand_counts_t fw, rc;
  query_windows(query, reference, hdist_th, nmers_limit, fw, rc);

  direction_t out;
  out.rows_v.reserve(nwins);
  vec<double> d_v;
  d_v.reserve(nwins);

  // One scored window, plus what its sample line needs. lr_bg depends on this
  // direction's median, so the sample lines are rendered in a second pass.
  struct scored_t
  {
    ds_t row;
    char strand;
    const uint64_t* hist;
    uint64_t u;
  };
  vec<scored_t> scored_v;
  scored_v.reserve(nwins);

  for (size_t wix = 0; wix < nwins; ++wix) {
    const uint64_t* hfw = fw.row(wix);
    const uint64_t t_fw = hist_total(hfw, hdist_th);
    const double d_fw = t_fw == 0 ? nanx() : llhf.mle(hfw, fw.u_v[wix]);

    scored_t s;
    s.row.d = d_fw;
    s.strand = canonical ? '.' : '+';
    s.hist = hfw;
    s.u = fw.u_v[wix];
    if (!canonical) {
      const uint64_t* hrc = rc.row(wix);
      const uint64_t t_rc = hist_total(hrc, hdist_th);
      const double d_rc = t_rc == 0 ? nanx() : llhf.mle(hrc, rc.u_v[wix]);
      const auto picked = select_strand_distance(d_fw, d_rc);
      s.row.d = picked.first;
      s.strand = picked.second;
      if (s.strand == '-') {
        s.hist = hrc;
        s.u = rc.u_v[wix];
      }
    }

    if (is_valid_distance(s.row.d)) {
      s.row.s = compute_lr_ub(llhf, s.row.d, s.u + hist_total(s.hist, hdist_th));
      d_v.push_back(s.row.d);
    }
    scored_v.push_back(s);
  }
  for (const scored_t& s : scored_v)
    out.rows_v.push_back(s.row);

  std::sort(d_v.begin(), d_v.end());
  out.d_median = linear_quantile(d_v, 0.5);
  out.n_valid = d_v.size();

  if (output_samples) {
    strstream ss;
    set_precision(ss, 5);
    for (size_t wix = 0; wix < nwins; ++wix) {
      const scored_t& s = scored_v[wix];
      double lr_bg = nanx();
      if (is_valid_distance(s.row.d) && is_valid_distance(out.d_median)) {
        lr_bg = likelihood_ratio_statistic(llhf.nll(out.d_median, s.hist, s.u), llhf.nll(s.row.d, s.hist, s.u));
      }
      write_tsv(ss,
                "gdiff",
                is_ba ? reference.get_rname() : query.get_rname(),
                is_ba ? query.get_rname() : reference.get_rname(),
                wins_v[wix].qid,
                wins_v[wix].start + 1,
                wins_v[wix].end + query.get_k() - 1,
                s.strand,
                is_ba ? "ba" : "ab",
                s.row.d,
                lr_bg,
                s.row.s)
        << '\n';
    }
    out.samples = ss.str();
  }
  return out;
}

WindowEstimator::WindowEstimator(const Sketch& sketch,
                                 const vec<qseq_t>& batch_v,
                                 uint64_t tau,
                                 uint64_t bin_shift,
                                 uint32_t hdist_th)
  : sketch(sketch)
  , batch_v(batch_v)
  , hdist_th(hdist_th)
  , tau(tau)
  , bin_shift(bin_shift)
  , nwinmers(window_nmers(tau, bin_shift))
  , llhf(make_llhf(sketch, hdist_th))
{
  canonical = sketch.is_canonical();
}

void WindowEstimator::estimate_all(uint64_t sample_size, ThreadPool& pool)
{
  plan_for_all(sample_size);
  evaluate(pool);
}

void WindowEstimator::estimate_per_sequence(uint64_t sample_size, ThreadPool& pool)
{
  plan_per_sequence(sample_size);
  evaluate(pool);
}

void WindowEstimator::plan_for_all(uint64_t sample_size)
{
  schemes_v.clear();
  vec<uint64_t> source_lens_v(batch_v.size());
  for (size_t i = 0; i < batch_v.size(); ++i)
    source_lens_v[i] = batch_v[i].seq.size();

  // One global sample, spread over the whole batch.
  window_plan_t wp = make_window_plan(source_lens_v, sketch.get_k(), tau, bin_shift, sample_size, gen);
  schemes_v.reserve(wp.sources_v.size());
  for (window_plan_t::source_t& src : wp.sources_v)
    schemes_v.emplace_back(src.bix, src.starts_v.size(), src.enmers, std::move(src.starts_v));
}

void WindowEstimator::plan_per_sequence(uint64_t sample_size)
{
  schemes_v.clear();
  schemes_v.reserve(batch_v.size());
  for (size_t bix = 0; bix < batch_v.size(); ++bix) {
    // A sample per source: one draw each, in source order.
    const vec<uint64_t> source_lens_v{batch_v[bix].seq.size()};
    window_plan_t wp = make_window_plan(source_lens_v, sketch.get_k(), tau, bin_shift, sample_size, gen);
    for (window_plan_t::source_t& src : wp.sources_v)
      schemes_v.emplace_back(bix, src.starts_v.size(), src.enmers, std::move(src.starts_v));
  }
}

void WindowEstimator::evaluate(ThreadPool& pool)
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

  pool.parallel_for(tasks_v.size(), 1, [&](uint64_t ti) {
    const task_t& t = tasks_v[ti];
    scheme_t& scheme = schemes_v[t.six];
    const char* cseq = batch_v[scheme.bix].seq.data();
    if (canonical) {
      window_counts_t agg(hdist_th);
      for (uint64_t s = t.a; s < t.b; ++s) {
        const uint64_t jx = scheme.starts_v[s] << bin_shift;
        const uint64_t jy = std::min(jx + nwinmers, scheme.enmers);
        agg.clear();
        scan_mers_range<true>(ctx, cseq, jx, jy, agg);
        const double d = llhf.mle(agg.hist(), agg.u);
        scheme.d_v[s] = d;
        scheme.strand_v[s] = '.';
        if (is_valid_distance(d)) {
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
        scan_mers_range<false>(ctx, cseq, jx, jy, agg);
        const double d_fw = llhf.mle(agg.hist_fw(), agg.u_fw);
        const double d_rc = llhf.mle(agg.hist_rc(), agg.u_rc);
        const auto [d, strand] = select_strand_distance(d_fw, d_rc);
        scheme.d_v[s] = d;
        scheme.strand_v[s] = strand;
        if (is_valid_distance(d)) {
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

ds_t WindowEstimator::make_row(double d, const uint64_t* hist, uint64_t u) const
{
  ds_t row;
  row.d = d;
  if (!hist || !is_valid_distance(d)) return row;
  row.s = compute_lr_ub(llhf, d, u + hist_total(hist, llhf.hdist_th));
  return row;
}

uint64_t WindowEstimator::nwins() const
{
  uint64_t n = 0;
  for (const auto& scheme : schemes_v)
    n += scheme.nsamples;
  return n;
}
