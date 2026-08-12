#include "detect.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>

#include "common.hpp"
#include "dist.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "scan.hpp"

extern uint32_t num_threads;

namespace {

  // Counts of the sampled background window closest to a target distance;
  // represents a typical window's information content for the threshold screen.
  struct median_window_t
  {
    std::array<uint64_t, hdist_bound + 1> hist{};
    uint64_t u = 0;
    double dev = std::numeric_limits<double>::infinity();
    bool valid = false;

    void consider(double d, const uint64_t* hist_in, uint64_t u_in, double target)
    {
      if (hist_in == nullptr || u_in == 0 || !is_valid_distance(d) || !is_valid_distance(target)) return;
      const double cand = std::abs(d - target);
      if (cand >= dev) return;
      std::copy(hist_in, hist_in + hdist_bound + 1, hist.begin());
      u = u_in;
      dev = cand;
      valid = true;
    }
  };

  // Single pass over the samples: per query with a finite dmed_v[bix], copy the
  // counts of the window whose distance is closest to that query's median.
  vec<median_window_t> find_median_windows(const DistanceSampler& sampler, const vec<double>& dmed_v)
  {
    vec<median_window_t> out_v(dmed_v.size());
    sampler.for_each_sample_counts([&](uint64_t bix, uint64_t, uint64_t, double d, char, const uint64_t* hist, uint64_t u) {
      if (bix >= dmed_v.size()) return;
      out_v[bix].consider(d, hist, u, dmed_v[bix]);
    });
    return out_v;
  }

  // Pooled-fit variant: one target median shared across all sampled windows.
  median_window_t find_median_window(const DistanceSampler& sampler, double d_median)
  {
    median_window_t out;
    if (!is_valid_distance(d_median)) return out;
    sampler.for_each_sample_counts([&](uint64_t, uint64_t, uint64_t, double d, char, const uint64_t* hist, uint64_t u) {
      out.consider(d, hist, u, d_median);
    });
    return out;
  }

  constexpr uint32_t lane_bit(size_t ix) { return uint32_t{1} << ix; }

  // chi-square(1) 99% critical value; same cut dist uses for median filtering.
  constexpr double lr_th_99 = 6.63;

  // True when d_median cannot be told apart from max_estimable_distance on a
  // window with n_total observed k-mers (lr_ub below lr_th). That ceiling is the
  // weakest homology the model can estimate — the practical lower limit of
  // usable signal — so a background parked there has nothing left to detect.
  bool median_at_estimable_floor(const LLH<double>& llhf, double d_median, uint64_t n_total, double* lr_out = nullptr)
  {
    const double lr = compute_lr_ub(llhf, d_median, n_total);
    if (lr_out) *lr_out = lr;
    return std::isfinite(lr) && lr < lr_th_99;
  }

  uint64_t window_n_total(const median_window_t& win, uint32_t hdist_th)
  {
    uint64_t n = win.u;
    for (uint32_t d = 0; d <= hdist_th; ++d)
      n += win.hist[d];
    return n;
  }

  // Sample mean and sample standard deviation; NaN components when empty.
  xy_t sample_mean_sd(const vec<double>& d_v)
  {
    if (d_v.empty()) return {nanx(), nanx()};
    const double mean = std::accumulate(d_v.begin(), d_v.end(), 0.0) / static_cast<double>(d_v.size());
    if (d_v.size() == 1) return {mean, 0.0};
    double sum_sq = 0.0;
    for (const double d : d_v) {
      const double delta = d - mean;
      sum_sq += delta * delta;
    }
    return {mean, std::sqrt(sum_sq / static_cast<double>(d_v.size() - 1))};
  }

  // SIMD params and LLH for one threshold set. Pooled mode builds this once and
  // shares it across queries (LLH is const after construction; safe to share).
  struct detect_ctx_t
  {
    params_t<cm512_t> params;
    llh_sptr_t<cm512_t> llhf;

    detect_ctx_t(const arr<double, RWIDTH>& extrema,
                 uint32_t k,
                 uint32_t h,
                 double rho,
                 uint32_t hdist_th,
                 uint64_t tau,
                 double chisq,
                 uint64_t bin_shift,
                 bool canonical)
      : params(extrema, hdist_th, tau, chisq, bin_shift, 0, canonical, false)
      , llhf(std::make_shared<LLH<cm512_t>>(k, h, rho, hdist_th, extrema))
    {
    }
  };

  // Interval extracted by one or more threshold slots (bit i = slot i). The
  // shared likelihood estimate is filled in after deduplication.
  struct candidate_t
  {
    uint64_t a_bin;
    uint64_t b_bin;
    uint32_t mask;
    likelihood_estimate_t est;
  };

  inline bool operator<(const candidate_t& x, const candidate_t& y)
  {
    return x.a_bin < y.a_bin || (x.a_bin == y.a_bin && x.b_bin < y.b_bin);
  }

  // Merge same-interval candidates from nested lanes; OR their lane masks.
  void merge_candidates(vec<candidate_t>& cv)
  {
    if (cv.empty()) return;
    std::sort(cv.begin(), cv.end());
    size_t w = 0;
    for (size_t r = 1; r < cv.size(); ++r) {
      if (cv[r].a_bin == cv[w].a_bin && cv[r].b_bin == cv[w].b_bin) {
        cv[w].mask |= cv[r].mask;
      } else {
        cv[++w] = cv[r];
      }
    }
    cv.resize(w + 1);
  }

  // LR screen for the high threshold only: drops high-side detection when a
  // typical background window cannot tell t_high apart from the fitted median.
  // Low-side detection is left alone — near-floor backgrounds still need it.
  bool screen_high_side(const clvl_t& level,
                        const LLH<double>& llhf,
                        const uint64_t* win_hist,
                        const uint64_t win_u,
                        const double d_median)
  {
    if (win_hist == nullptr || win_u == 0 || !is_valid_distance(d_median)) return true;
    const double crit = GammaModel::quantile(1.0 - level.alpha, 0.5, 2.0); // chi-square(1) == Gamma(1/2, 2)
    const double nll_med = llhf.nll(d_median, win_hist, win_u);
    const double lr_high = likelihood_ratio_statistic(llhf.nll(level.t_high, win_hist, win_u), nll_med);
    if (!std::isfinite(lr_high) || lr_high < crit) {
      warn_msg(concat_msg("level ",
                          level.alpha,
                          " high threshold indistinguishable from the background median (LR: high=",
                          lr_high,
                          ", chi2 crit=",
                          crit,
                          "); dropping high-side detection at this level"));
      return false;
    }
    return true;
  }

  // Writes one TSV row for an extracted interval (outlier or unmapped segment).
  void write_outlier(strstream& os,
                     const str& qid,
                     uint64_t L,
                     const interval_t& seq_iv,
                     const str& rname,
                     bool canonical,
                     bool is_rc,
                     double d_diff,
                     double d_q,
                     const likelihood_estimate_t& est,
                     bool high_side,
                     const vec<double>& side_thresholds_v,
                     uint32_t mask,
                     double alpha)
  {
    const char* side = est.has_hits ? (high_side ? "high" : "low") : "unmapped";
    const xy_t d_range = bracket_distance(est.d, side_thresholds_v);
    std::ostringstream d_bin;
    d_bin.flags(os.flags());
    d_bin.precision(os.precision());
    d_bin << '(' << d_range.first << ", " << d_range.second << ')';

    if (canonical) {
      write_tsv(
        os, qid, L, seq_iv.a, seq_iv.b, rname, est.d, side, mask, alpha, d_bin.str(), d_q, est.I, est.lr_bg, est.lr_ub)
        << '\n';
    } else {
      write_tsv(os,
                qid,
                L,
                seq_iv.a,
                seq_iv.b,
                report_strand(is_rc, d_diff),
                static_cast<uint32_t>(is_rc),
                rname,
                est.d,
                side,
                mask,
                alpha,
                d_bin.str(),
                d_q,
                d_diff,
                est.I,
                est.lr_bg,
                est.lr_ub)
        << '\n';
    }
  }

} // namespace

Detector::Detector(const sketch_sptr_t& sketch,
                   const vec<qseq_t>& batch_v,
                   const uint64_t tau,
                   const uint64_t bin_shift,
                   const uint32_t hdist_th,
                   const double chisq,
                   const uint64_t sample_size,
                   const vec<double>& levels,
                   const vec<double>& fit_quantiles,
                   const bool per_sequence,
                   const uint32_t verbosity)
  : sketch(sketch)
  , batch_v(batch_v)
  , tau(tau)
  , bin_shift(bin_shift)
  , hdist_th(hdist_th)
  , chisq(chisq)
  , sample_size(sample_size)
  , levels(levels)
  , fit_quantiles(fit_quantiles)
  , per_sequence(per_sequence)
  , verbosity(verbosity)
{
}

bggamma_t Detector::fit(const vec<double>& d_v) const
{
  bggamma_t fit;
  const auto prepared = GammaModel::prepare_samples(d_v, d_eps);
  fit.ndropped = prepared.ndropped;
  fit.nfloored = prepared.nfloored;
  fit.nsamples = prepared.x.size();
  if (prepared.x.size() < GammaModel::min_nsamples) return fit;
  if (fit.nfloored * 20 > fit.nsamples) {
    warn_msg(concat_msg(
      "more than 5% floored windows (", fit.nfloored, "/", fit.nsamples, "); low-side thresholds may be meaningless"));
  }

  GammaModel::Config cfg{};
  cfg.quantile_probs = {fit_quantiles[0], fit_quantiles[1], fit_quantiles[2]};
  fit.params = GammaModel::fit_from_samples(prepared.x, cfg, &fit.objective, &fit.niter);
  if (!GammaModel::validate_params(fit.params)) {
    fit.params = GammaModel::moments_estimate(prepared.x);
    if (!GammaModel::validate_params(fit.params)) return fit;
    // Diagnostics from the failed Nelder-Mead run do not describe this estimate.
    fit.objective = nanx();
    fit.niter = 0;
    warn_msg("Nelder-Mead gamma fit failed; falling back to the moment estimate");
  }
  if (fit.params.shape > 1e6) {
    fit.params.scale *= fit.params.shape / 1e6;
    fit.params.shape = 1e6;
    warn_msg("near-degenerate background (all samples nearly identical); clamping gamma shape while preserving its mean");
  }
  const double mass_above = 1.0 - GammaModel::cdf(d_ub, fit.params.shape, fit.params.scale);
  if (mass_above > 1e-3) {
    warn_msg(concat_msg(
      "fitted gamma puts ", mass_above * 100.0, "% of its mass above the maximum distance; tail thresholds are biased low"));
  }
  fit.ok = true;
  return fit;
}

thcfg_t Detector::thresholds_for(const bggamma_t& fit,
                                 const LLH<double>& llhf,
                                 const uint64_t* win_hist,
                                 const uint64_t win_u,
                                 const double d_median,
                                 const bool disable_high) const
{
  vec<double> sorted_levels = levels;
  std::sort(sorted_levels.begin(), sorted_levels.end(), std::greater{});

  thcfg_t thresholds;
  vec<double> used;
  used.reserve(2 * sorted_levels.size());
  for (const double alpha : sorted_levels) {
    clvl_t level{alpha,
                 GammaModel::quantile(alpha / 2.0, fit.params.shape, fit.params.scale),
                 GammaModel::quantile(1.0 - alpha / 2.0, fit.params.shape, fit.params.scale),
                 true};
    level.t_low = std::clamp(level.t_low, d_eps, d_ub - d_eps);
    level.t_high = std::clamp(level.t_high, d_eps, d_ub - d_eps);
    if (!std::isfinite(level.t_low) || !std::isfinite(level.t_high) || !(level.t_low < level.t_high)) {
      warn_msg(concat_msg("degenerate thresholds at level ", alpha, "; dropping the level"));
      continue;
    }
    level.high = !disable_high && screen_high_side(level, llhf, win_hist, win_u, d_median);
    const bool duplicate = std::any_of(used.begin(), used.end(), [&](double t) {
      return std::abs(t - level.t_low) < 1e-12 || (level.high && std::abs(t - level.t_high) < 1e-12);
    });
    if (duplicate) {
      warn_msg(concat_msg("duplicate thresholds at level ", level.alpha, " after clipping; dropping the level"));
      continue;
    }
    used.push_back(level.t_low);
    if (level.high) used.push_back(level.t_high);
    thresholds.levels.push_back(level);
  }
  if (!thresholds.empty()) thresholds.pack();
  return thresholds;
}

vec<thcfg_t> Detector::plan(const DistanceSampler& sampler, const vvec<double>& d_per_seq, const LLH<double>& llhf) const
{
  const str rname = sketch->get_rname();
  const size_t nseq = batch_v.size();

  uint64_t nmapped = 0;
  for (const auto& d_v : d_per_seq)
    nmapped += d_v.size();
  const uint64_t nunmapped = sampler.get_nsamples() - nmapped;

  if (!per_sequence) {
    vec<double> d_all;
    d_all.reserve(nmapped);
    for (const auto& d_v : d_per_seq)
      d_all.insert(d_all.end(), d_v.begin(), d_v.end());

    const bggamma_t bg = fit(d_all);
    if (!bg.ok) {
      error_exit(concat_msg("gamma fit failed for sketch ",
                            rname,
                            " (",
                            bg.nsamples,
                            " usable samples; need at least ",
                            GammaModel::min_nsamples,
                            ")"));
    }

    const double d_median = GammaModel::quantile(0.5, bg.params.shape, bg.params.scale);
    const median_window_t win = find_median_window(sampler, d_median);
    bool disable_high = false;
    if (win.valid) {
      double lr_ub = nanx();
      if (median_at_estimable_floor(llhf, d_median, window_n_total(win, hdist_th), &lr_ub)) {
        disable_high = true;
        warn_msg(concat_msg("background median indistinguishable from max estimable distance for sketch ",
                            rname,
                            " (d_median=",
                            d_median,
                            ", lr_ub=",
                            lr_ub,
                            ", chi2 crit=",
                            lr_th_99,
                            "); disabling high-side detection"));
      }
    }

    thcfg_t thresholds =
      thresholds_for(bg, llhf, win.valid ? win.hist.data() : nullptr, win.u, d_median, disable_high);
    if (thresholds.empty()) {
      error_exit(concat_msg("no usable confidence levels for sketch ", rname));
    }

    if (verbosity >= 1) {
      const auto [mean, sd] = sample_mean_sd(d_all);
      report_fit(bg, thresholds, mean, sd, nunmapped, nmapped);
    }

    return {std::move(thresholds)};
  }

  // Per-sequence: fit all sequences first so every median is known, then gather
  // each sequence's representative (median-closest) window in one pass.
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  const uint64_t nwinmers = ((tau + bin_size - 1) >> bin_shift) << bin_shift;
  const uint64_t min_len = nwinmers + sketch->get_lshf_sptr()->get_k() - 1;

  vec<bggamma_t> fits(nseq);
  vec<double> dmed_v(nseq, nanx());
  for (size_t bix = 0; bix < nseq; ++bix) {
    fits[bix] = fit(d_per_seq[bix]);
    if (fits[bix].ok) dmed_v[bix] = GammaModel::quantile(0.5, fits[bix].params.shape, fits[bix].params.scale);
  }
  const vec<median_window_t> win_v = find_median_windows(sampler, dmed_v);

  vec<thcfg_t> sets(nseq);
  uint64_t nfit = 0;
  uint64_t nskipped = 0;
  for (size_t bix = 0; bix < nseq; ++bix) {
    if (fits[bix].ok) {
      const median_window_t& win = win_v[bix];
      bool disable_high = false;
      if (win.valid && median_at_estimable_floor(llhf, dmed_v[bix], window_n_total(win, hdist_th))) {
        disable_high = true;
        if (batch_v[bix].seq.size() >= min_len) {
          warn_pmsg(batch_v[bix].qid,
                    "background median indistinguishable from max estimable distance; disabling high-side detection");
        }
      }
      sets[bix] =
        thresholds_for(fits[bix], llhf, win.valid ? win.hist.data() : nullptr, win.u, dmed_v[bix], disable_high);
    }
    if (sets[bix].empty() || sets[bix].nlevels() != levels.size()) {
      sets[bix] = {};
      if (batch_v[bix].seq.size() >= min_len) {
        warn_pmsg(batch_v[bix].qid, "gamma fit failed or degenerate; detection skipped for this sequence");
        ++nskipped;
      }
      continue;
    }
    ++nfit;
  }

  if (verbosity >= 1) {
    cerr_msg("[",
             rname,
             "] per-sequence fits: ",
             nfit,
             " ok, ",
             nskipped,
             " skipped",
             " (unmapped windows: ",
             nunmapped,
             ", hit: ",
             nmapped,
             ")");
  }
  return sets;
}

void Detector::extract_batch(const vec<thcfg_t>& sets,
                             const bool per_sequence,
                             const LLH<double>& llhf,
                             vec<lvlstat_t>& stats,
                             uint64_t& unmapped_iv,
                             uint64_t& unmapped_bp,
                             strstream& sout,
                             ThreadPool& pool) const
{
  const lshf_sptr_t lshf = sketch->get_lshf_sptr();
  const uint32_t k = lshf->get_k();
  const bool canonical = sketch->is_canonical();
  const str rname = sketch->get_rname();
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  const uint64_t tau_bin = std::max<uint64_t>(1, (tau + bin_size - 1) >> bin_shift);

  const scan_ctx_t scan_ctx = make_scan_ctx(*sketch, bin_shift, hdist_th);

  const size_t nseq = batch_v.size();
  vec<strstream> out_v(nseq);
  vvec<lvlstat_t> stats_v(nseq);
  vec<uint64_t> unmapped_iv_v(nseq, 0);
  vec<uint64_t> unmapped_bp_v(nseq, 0);

  // Pooled mode shares one threshold set across all queries.
  std::optional<detect_ctx_t> pooled_ctx;
  if (!per_sequence) {
    assert(!sets.empty() && !sets[0].empty());
    pooled_ctx.emplace(sets[0].extrema, k, lshf->get_h(), sketch->get_rho(), hdist_th, tau, chisq, bin_shift, canonical);
  }

  auto process = [&](const uint64_t bix) {
    const thcfg_t& thresholds = sets[per_sequence ? bix : 0];
    if (thresholds.empty()) return;

    const char* cseq = batch_v[bix].seq.data();
    const uint64_t len = batch_v[bix].seq.size();
    if (len < static_cast<uint64_t>(k)) {
      warn_pmsg(batch_v[bix].qid, "skipped: sequence shorter than k-mer length ", "(len=", len, ", k=", k, ")");
      return;
    }
    const uint64_t enmers = len - k + 1;
    const uint64_t nbins = (enmers + bin_size - 1) >> bin_shift;
    if (nbins < 2) {
      warn_pmsg(batch_v[bix].qid, "skipped: fewer than two bins after binning ", "(len=", len, ", bin_size=", bin_size, ")");
      return;
    }
    const uint64_t tau_eff = std::min(tau_bin, nbins) - 1;

    // Per-sequence mode: thresholds differ per query, so params/LLH are built here.
    std::optional<detect_ctx_t> sequence_ctx;
    if (!pooled_ctx) {
      sequence_ctx.emplace(
        thresholds.extrema, k, lshf->get_h(), sketch->get_rho(), hdist_th, tau, chisq, bin_shift, canonical);
    }
    const detect_ctx_t& ctx = pooled_ctx ? *pooled_ctx : *sequence_ctx;

    strstream& os = out_v[bix];
    set_precision(os, 5);
    stats_v[bix].assign(thresholds.nlanes(), {});

    DIM<cm512_t> dim_fw(ctx.params, ctx.llhf, nbins, enmers);
    std::optional<DIM<cm512_t>> dim_rc;
    if (!canonical) dim_rc.emplace(ctx.params, ctx.llhf, nbins, enmers);
    dim_agg_t<cm512_t> agg{dim_fw, canonical ? nullptr : &*dim_rc};
    if (canonical)
      scan_mers_range<false>(scan_ctx, cseq, 0, enmers, agg);
    else
      scan_mers_range<true>(scan_ctx, cseq, 0, enmers, agg);

    // Uniform per-strand pipeline: prefix sums, whole-query distance, then
    // per-threshold extraction. Canonical mode runs a single forward "strand".
    struct strand_t
    {
      DIM<cm512_t>* dim;
      bool is_rc;
      double d_q = nanx();
    };
    arr<strand_t, 2> strands{{strand_t{&dim_fw, false}, strand_t{canonical ? nullptr : &*dim_rc, true}}};
    const size_t nstrands = canonical ? 1 : 2;

    for (size_t si = 0; si < nstrands; ++si) {
      strand_t& strand = strands[si];
      strand.dim->inclusive_scan();
      strand.dim->compute_prefhistsum();
      strand.dim->extrema_scan();
      vec<uint64_t> v_q;
      uint64_t u_q = 0, t_q = 0;
      strand.dim->total_histogram(v_q, u_q, t_q);
      strand.d_q = llhf.mle(v_q.data(), u_q);
    }
    const double d_diff = canonical ? nanx() : strand_diff(strands[0].d_q, strands[1].d_q);

    // Collect candidates from every lane, merge nested duplicates, score once.
    arr<vec<candidate_t>, 2> candidates_v;
    for (size_t ix = 0; ix < thresholds.nlanes(); ++ix) {
      if (thresholds.is_high_side(ix) && !thresholds.high_enabled(ix)) continue;
      for (size_t si = 0; si < nstrands; ++si) {
        DIM<cm512_t>& dim = *strands[si].dim;
        dim.extract_intervals_mx(tau_eff, 1, nbins, ix);
        dim.expand_intervals(ctx.params.chisq, ix);
        for (const auto& iv : dim.get_intervals_v(ix))
          candidates_v[si].push_back({iv.a, iv.b + 1, lane_bit(ix), {}});
      }
    }

    for (size_t si = 0; si < nstrands; ++si) {
      vec<candidate_t>& cv = candidates_v[si];
      merge_candidates(cv);

      vec<uint64_t> scratch_v;
      uint64_t u = 0, t = 0;
      for (auto& c : cv) {
        strands[si].dim->extract_histogram(c.a_bin - 1, c.b_bin - 1, scratch_v, u, t);
        c.est = compute_likelihood_estimate(llhf, scratch_v.data(), u, t, strands[si].d_q);
        if (!c.est.has_hits) {
          // Count unmapped intervals once per unique (strand, a_bin, b_bin).
          const interval_t seq_iv = get_coordinates({c.a_bin, c.b_bin}, bin_shift, enmers, k);
          unmapped_iv_v[bix] += 1;
          unmapped_bp_v[bix] += seq_iv.b - seq_iv.a + 1;
        }
      }
    }

    // Emit from merged candidates by lane mask (lane-major, then strand).
    const uint64_t L = enmers + k - 1;
    for (size_t ix = 0; ix < thresholds.nlanes(); ++ix) {
      if (thresholds.is_high_side(ix) && !thresholds.high_enabled(ix)) continue;
      const uint32_t bit = lane_bit(ix);
      const bool high_side = thresholds.is_high_side(ix);
      const vec<double>& side_thresholds_v = high_side ? thresholds.high_v : thresholds.low_v;
      for (size_t si = 0; si < nstrands; ++si) {
        for (const candidate_t& c : candidates_v[si]) {
          if ((c.mask & bit) == 0) continue;
          const interval_t seq_iv = get_coordinates({c.a_bin, c.b_bin}, bin_shift, enmers, k);
          if (c.est.has_hits) {
            stats_v[bix][ix].nintervals += 1;
            stats_v[bix][ix].bp_covered += seq_iv.b - seq_iv.a + 1;
          }
          write_outlier(os,
                        batch_v[bix].qid,
                        L,
                        seq_iv,
                        rname,
                        canonical,
                        strands[si].is_rc,
                        d_diff,
                        strands[si].d_q,
                        c.est,
                        high_side,
                        side_thresholds_v,
                        bit,
                        thresholds.alpha(ix));
        }
      }
    }
  };

  pool.parallel_for(nseq, 1, process);

  unmapped_iv = 0;
  unmapped_bp = 0;
  for (size_t bix = 0; bix < nseq; ++bix) {
    if (out_v[bix].tellp() > 0) sout << out_v[bix].rdbuf();
    assert(stats_v[bix].size() <= stats.size());
    for (size_t ix = 0; ix < stats_v[bix].size(); ++ix) {
      stats[ix].nintervals += stats_v[bix][ix].nintervals;
      stats[ix].bp_covered += stats_v[bix][ix].bp_covered;
    }
    unmapped_iv += unmapped_iv_v[bix];
    unmapped_bp += unmapped_bp_v[bix];
  }
}

void Detector::report_fit(const bggamma_t& fit,
                          const thcfg_t& thresholds,
                          const double mean,
                          const double sd,
                          const uint64_t nunmapped,
                          const uint64_t nmapped) const
{
  const str rname = sketch->get_rname();
  cerr_msg("[",
           rname,
           "] background windows: n=",
           fit.nsamples,
           " (floored=",
           fit.nfloored,
           ", dropped=",
           fit.ndropped,
           ", unmapped=",
           nunmapped,
           ", hit=",
           nmapped,
           ")",
           " mean=",
           mean,
           " sd=",
           sd);
  if (verbosity >= 2) {
    cerr_msg("[",
             rname,
             "] fit: shape=",
             fit.params.shape,
             " scale=",
             fit.params.scale,
             " objective=",
             fit.objective,
             " niter=",
             fit.niter,
             " median=",
             GammaModel::quantile(0.5, fit.params.shape, fit.params.scale));
  } else {
    cerr_msg("[", rname, "] fit: shape=", fit.params.shape, " scale=", fit.params.scale);
  }
  for (const auto& level : thresholds.levels) {
    if (level.high)
      cerr_msg("[", rname, "] level ", level.alpha, ": t_low=", level.t_low, " t_high=", level.t_high);
    else
      cerr_msg("[", rname, "] level ", level.alpha, ": t_low=", level.t_low, " t_high=skipped");
  }
}

void Detector::report_stats(const thcfg_t& thresholds,
                            const vec<lvlstat_t>& stats,
                            const uint64_t unmapped_iv,
                            const uint64_t unmapped_bp) const
{
  const str rname = sketch->get_rname();
  for (size_t j = 0; j < thresholds.nlevels(); ++j) {
    if (thresholds.high_enabled(j)) {
      cerr_msg("[",
               rname,
               "] level ",
               thresholds.levels[j].alpha,
               ": high: ",
               stats[j].nintervals,
               " interval(s), ",
               stats[j].bp_covered,
               " bp | low: ",
               stats[thresholds.nlevels() + j].nintervals,
               " interval(s), ",
               stats[thresholds.nlevels() + j].bp_covered,
               " bp");
    } else {
      cerr_msg("[",
               rname,
               "] level ",
               thresholds.levels[j].alpha,
               ": high: skipped | low: ",
               stats[thresholds.nlevels() + j].nintervals,
               " interval(s), ",
               stats[thresholds.nlevels() + j].bp_covered,
               " bp");
    }
  }
  cerr_msg("[", rname, "] unmapped (no k-mer hits): ", unmapped_iv, " interval(s), ", unmapped_bp, " bp");
}

void Detector::run(std::ostream& out, ThreadPool& pool)
{
  const lshf_sptr_t lshf = sketch->get_lshf_sptr();
  // Scalar LLH for the threshold screen and the per-interval estimates.
  const LLH<double> llhf(lshf->get_k(), lshf->get_h(), sketch->get_rho(), hdist_th, 0.0, false);

  // Pass 1: sample background windows. Counts are kept so that a representative
  // window per fit can feed the threshold screen.
  DistanceSampler sampler(sketch, batch_v, tau, bin_shift, hdist_th);
  if (per_sequence)
    sampler.run_per_sequence(sample_size, true, pool);
  else
    sampler.run_for_all(sample_size, true, pool);
  vvec<double> d_per_seq(batch_v.size());
  sampler.collect_distances(d_per_seq);

  const vec<thcfg_t> sets = plan(sampler, d_per_seq, llhf);

  // Pass 2: extract outlier intervals per threshold.
  size_t max_lanes = 0;
  for (const auto& set : sets)
    max_lanes = std::max(max_lanes, set.nlanes());
  vec<lvlstat_t> stats(max_lanes);
  uint64_t unmapped_iv = 0, unmapped_bp = 0;
  strstream sout;
  set_precision(sout, 5);
  extract_batch(sets, per_sequence, llhf, stats, unmapped_iv, unmapped_bp, sout, pool);
  if (sout.tellp() > 0) out << sout.rdbuf();

  if (verbosity >= 1 && !per_sequence && !sets.empty()) {
    report_stats(sets[0], stats, unmapped_iv, unmapped_bp);
  }
}

DetectSC::DetectSC(CLI::App& sc)
{
  sc.add_option("target-path", target_path, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible)")
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("sketch-path", sketch_path, "Reference sketch file <path>")->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_option("-l", tau, "Minimum interval length in k-mers; also the background window length")
    ->required()
    ->check(CLI::PositiveNumber);
  sc.add_option("-b,--bin-shift", bin_shift, "Group consecutive k-mers into bins of size 2^b [0]")->check(CLI::Range(0, 16));
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("--chisq", chisq, "Chi-square threshold for merging overlapping intervals [33.00051]")
    ->check(CLI::NonNegativeNumber);
  sc.add_option("--sample-size", sample_size, "Background windows sampled per fit [1000]")->check(CLI::PositiveNumber);
  sc.add_option("--levels", levels, "Two-sided confidence level(s) - provide exactly 1 or 4 values [0.05 0.01 0.001 0.0001]")
    ->expected(1, 4);
  sc.add_flag("--per-sequence", per_sequence, "Fit and sample the background per query sequence [per-sketch]");
  sc.add_option("--fit-quantiles", fit_quantiles, "Three central quantile probabilities for the gamma fit [0.2 0.4 0.6]")
    ->expected(3);
  sc.add_option("--verbosity", verbosity, "Statistics and report detail level (0-3) [1]")->check(CLI::Range(0, 3));
  sc.callback([&]() {
    if (levels.empty()) levels = {0.05, 0.01, 0.001, 0.0001};
    if (fit_quantiles.empty()) fit_quantiles = {0.2, 0.4, 0.6};
    if (!validate_configuration()) {
      error_exit("Invalid configuration!");
    }
    if (!output_path.empty()) {
      output_file.open(output_path);
      check_fstream(output_file, "Cannot open output file", output_path.string());
      output_stream = &output_file;
    }
  });
}

bool DetectSC::validate_configuration()
{
  bool is_invalid = false;
  if (levels.size() != 1 && levels.size() != 4) {
    is_invalid = true;
    cerr_msg("--levels requires exactly 1 or 4 levels; got ", levels.size());
  }
  for (size_t i = 0; i < levels.size(); ++i) {
    if (!std::isfinite(levels[i]) || levels[i] <= 1e-8 || levels[i] >= 0.5) {
      is_invalid = true;
      cerr_msg("--levels[", i, "] must be in (1e-8, 0.5): ", levels[i]);
    }
  }
  {
    auto slevels = levels;
    std::sort(slevels.begin(), slevels.end());
    if (const auto it = std::adjacent_find(slevels.begin(), slevels.end()); it != slevels.end()) {
      is_invalid = true;
      cerr_msg("--levels values must be unique; duplicate: ", *it);
    }
  }
  if (!std::is_sorted(fit_quantiles.begin(), fit_quantiles.end()) ||
      std::adjacent_find(fit_quantiles.begin(), fit_quantiles.end()) != fit_quantiles.end()) {
    is_invalid = true;
    cerr_msg("--fit-quantiles must be strictly increasing");
  }
  for (size_t i = 0; i < fit_quantiles.size(); ++i) {
    if (!std::isfinite(fit_quantiles[i]) || fit_quantiles[i] <= 0.0 || fit_quantiles[i] >= 1.0) {
      is_invalid = true;
      cerr_msg("--fit-quantiles[", i, "] must be in (0, 1): ", fit_quantiles[i]);
    }
  }
  if (!validate_binning(bin_shift, tau)) {
    is_invalid = true;
  }
  const uint64_t bin_size = (bin_shift <= 16) ? (uint64_t(1) << bin_shift) : 0;
  const uint64_t tau_bin = (bin_size > 0) ? ((tau + bin_size - 1) >> bin_shift) : 0;
  if (tau_bin < 2) {
    is_invalid = true;
    cerr_msg("-l must span at least two bins after binning ", "(tau=", tau, ", bin_size=", bin_size, ")");
  }
  return !is_invalid;
}

void DetectSC::detect()
{
  set_precision(*output_stream, 5);

  qseq_sptr_t qs = std::make_shared<QSeq>(target_path);
  while (qs->read_next_batch()) {
  }
  const auto& batch_v = qs->get_batch_v();

  const vec<uint64_t> sketch_offsets = read_sketch_offsets(sketch_path);
  const uint32_t nsketches = static_cast<uint32_t>(sketch_offsets.size());

  ThreadPool pool(num_threads);
  cerr_msg("Processing ", nsketches, " sketches w/ ", pool.size(), " thread(s)...");

  init_thread_rng(1);

  std::ifstream sin(sketch_path, std::ifstream::binary);
  check_fstream(sin, "Cannot open sketch file for reading", sketch_path.string());

  for (uint32_t i = 0; i < nsketches; ++i) {
    sketch_sptr_t sketch = std::make_shared<Sketch>(sketch_path);
    sketch->load_from_offset(sin, sketch_offsets[i]);
    Detector detector(
      sketch, batch_v, tau, bin_shift, hdist_th, chisq, sample_size, levels, fit_quantiles, per_sequence, verbosity);
    detector.run(*output_stream, pool);

    if (verbosity >= 1) {
      std::cerr << "\rProcessed sketch " << i + 1 << "/" << nsketches << "..." << std::flush;
      if (i + 1 == nsketches) std::cerr << std::endl;
    }
  }
}
