#include "map.hpp"

#include "common.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "records.hpp"
#include "scan.hpp"
#include "windows.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <boost/math/distributions/gamma.hpp>
#include <simde/x86/avx512.h>
#include <sys/mman.h>

namespace {
  // The accumulator is one 512-bit register wide, so the histogram must have exactly 8 lanes.
  static_assert(hdist_bound + 1 == 8, "add_to_acc assumes an 8-lane histogram");

  inline void add_to_acc(vec<uint64_t>& acc_v, uint64_t& u_acc, const vec<uint64_t>& source_v, uint64_t u)
  {
    simde__m512i s = simde_mm512_loadu_si512(acc_v.data());
    s = simde_mm512_add_epi64(s, simde_mm512_loadu_si512(source_v.data()));
    simde_mm512_storeu_si512(acc_v.data(), s);
    u_acc += u;
  }

  vec<double> sorted_unique(vec<double> v)
  {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
  }

  // Threshold values as one scalar or a padded SIMD lane vector.
  template<typename T>
  T make_dist_th(const vec<double>& th_v)
  {
    if constexpr (std::is_same_v<T, double>) {
      return th_v.front();
    } else {
      cmlane_t lanes{};
      lanes.fill(d_eps);
      for (size_t i = 0; i < th_v.size(); ++i)
        lanes[i] = th_v[i];
      return lanes;
    }
  }

  inline void write_map_header(std::ostream& sout, bool canonical)
  {
    if (canonical) {
      write_tsv(sout,
                "seq",
                "seq_len",
                "start",
                "end",
                "reference",
                "d",
                "mask",
                "d_bin",
                "d_q",
                "d_acc",
                "pvalue",
                "fold",
                "qvalue",
                "fisher_info",
                "lr_ub")
        << '\n';
    } else {
      write_tsv(sout,
                "seq",
                "seq_len",
                "start",
                "end",
                "strand",
                "is_rc",
                "reference",
                "d",
                "mask",
                "d_bin",
                "d_q",
                "d_diff",
                "d_acc",
                "pvalue",
                "fold",
                "qvalue",
                "fisher_info",
                "lr_ub")
        << '\n';
    }
  }

  // One output line; canonical mode drops the strand columns. The column order is a published
  // contract (`plot.py` reads it positionally), so keep the two branches in step.
  inline void write_map_line(std::ostream& sout,
                             const str& qid,
                             bool canonical,
                             const str& rname,
                             const record_t& r,
                             mask_t mask,
                             const bracket_t& d_bin,
                             double d_acc)
  {
    if (canonical) {
      write_tsv(sout,
                qid,
                r.L,
                r.seq_iv.a,
                r.seq_iv.b,
                rname,
                r.d,
                mask,
                d_bin,
                r.d_q,
                d_acc,
                r.pvalue,
                r.fold,
                r.qvalue,
                r.fisher_info,
                r.lr_ub)
        << '\n';
    } else {
      write_tsv(sout,
                qid,
                r.L,
                r.seq_iv.a,
                r.seq_iv.b,
                report_strand(r.is_rc, r.d_diff),
                static_cast<uint32_t>(r.is_rc),
                rname,
                r.d,
                mask,
                d_bin,
                r.d_q,
                r.d_diff,
                d_acc,
                r.pvalue,
                r.fold,
                r.qvalue,
                r.fisher_info,
                r.lr_ub)
        << '\n';
    }
  }
} // namespace

// Empirical two-sided thresholds for exactly four levels over a sorted background pool.
// Each alpha/2 quantile is read from the full pool (floored samples included), then its
// value is lifted to the next distinct background distance above d_eps; each later level
// takes the next distinct distance after the previous threshold. The four upper quantiles
// are handled symmetrically. Succeeds only when this yields eight distinct in-range values.
vec<double> lr_filtered_sample(const vec<sample_point_t>& points_v, double lr_th, double min_portion)
{
  // The same rule `dist` applies to its reconciled windows: a window's likelihood-ratio bound is
  // its evidence, and the filtered set is only used when enough windows carry it.
  vec<double> all_v;
  vec<double> kept_v;
  for (const sample_point_t& p : points_v) {
    if (!is_valid_distance(p.d)) continue;
    all_v.push_back(p.d);
    if (!std::isnan(p.s) && p.s > lr_th) kept_v.push_back(p.d);
  }
  const double portion = all_v.empty() ? 0.0 : static_cast<double>(kept_v.size()) / static_cast<double>(all_v.size());
  if (portion > min_portion) return kept_v;
  return all_v;
}

bool fit_gamma_mom(const vec<double>& sample_v, gamma_fit_t& fit)
{
  fit = {};
  if (sample_v.size() < min_gamma_samples) return false;

  const double n = static_cast<double>(sample_v.size());
  double mean = 0.0;
  for (const double x : sample_v)
    mean += x;
  mean /= n;
  if (!(mean > 0.0)) return false;

  double var = 0.0;
  for (const double x : sample_v) {
    const double dx = x - mean;
    var += dx * dx;
  }
  var /= n - 1.0;
  if (!(var > 0.0)) return false;

  fit.shape = (mean * mean) / var;
  fit.scale = var / mean;
  return std::isfinite(fit.shape) && std::isfinite(fit.scale) && fit.shape > 0.0 && fit.scale > 0.0;
}

namespace {
  // Quantiles of the fitted gamma, lifted to the nearest observed distance when a tail reaches the
  // estimable floor. `lifted` reports whether that happened, so the caller warns only when the
  // gamma path is the one actually accepted.
  bool gamma_thresholds(const gamma_fit_t& fit,
                        const vec<double>& sample_v,
                        const vec<double>& levels_v,
                        vec<double>& out_v,
                        bool& lifted)
  {
    const boost::math::gamma_distribution<double> gamma(fit.shape, fit.scale);
    vec<double> levels = levels_v;
    std::sort(levels.begin(), levels.end()); // tightest alpha first

    lifted = false;
    double prev_low = d_eps;
    double prev_high = d_ub - d_eps;
    for (const double alpha : levels) {
      double t_low = 0.0;
      double t_high = 0.0;
      try {
        t_low = boost::math::quantile(gamma, alpha / 2.0);
        t_high = boost::math::quantile(gamma, 1.0 - alpha / 2.0);
      } catch (const std::exception&) {
        return false;
      }
      if (!std::isfinite(t_low) || !std::isfinite(t_high)) return false;

      if (!(t_low > prev_low)) {
        const auto it = std::upper_bound(sample_v.begin(), sample_v.end(), prev_low);
        if (it == sample_v.end()) return false;
        t_low = *it;
        lifted = true;
      }
      if (!(t_high < prev_high)) {
        const auto it = std::lower_bound(sample_v.begin(), sample_v.end(), prev_high);
        if (it == sample_v.begin()) return false;
        t_high = *std::prev(it);
        lifted = true;
      }
      if (!(t_low < t_high) || !(t_low > d_eps) || !(t_high < d_ub - d_eps)) return false;
      out_v.push_back(t_low);
      out_v.push_back(t_high);
      prev_low = t_low;
      prev_high = t_high;
    }
    std::sort(out_v.begin(), out_v.end());
    return true;
  }
} // namespace

bool empirical_thresholds_from_levels(const vec<double>& sample_v,
                                      const vec<double>& levels_v,
                                      vec<double>& out_v,
                                      size_t min_samples,
                                      bool warn_lift)
{
  out_v.clear();
  if (sample_v.size() < min_samples || levels_v.size() != 4) return false;

  vec<double> levels = levels_v;
  std::sort(levels.begin(), levels.end()); // tightest alpha first

  bool lifted = false;
  double prev_low = d_eps;
  double prev_high = d_ub - d_eps;
  for (const double alpha : levels) {
    double t_low = linear_quantile(sample_v, alpha / 2.0);
    if (!(t_low > prev_low)) {
      const auto it = std::upper_bound(sample_v.begin(), sample_v.end(), prev_low);
      if (it == sample_v.end()) return false;
      t_low = *it;
      lifted = true;
    }
    double t_high = linear_quantile(sample_v, 1.0 - alpha / 2.0);
    if (!(t_high < prev_high)) {
      const auto it = std::lower_bound(sample_v.begin(), sample_v.end(), prev_high);
      if (it == sample_v.begin()) return false;
      t_high = *std::prev(it);
      lifted = true;
    }
    if (!(t_low < t_high) || !(t_low > d_eps) || !(t_high < d_ub - d_eps)) return false;
    out_v.push_back(t_low);
    out_v.push_back(t_high);
    prev_low = t_low;
    prev_high = t_high;
  }
  std::sort(out_v.begin(), out_v.end());
  if (lifted && warn_lift)
    warn_msg("--levels: some quantiles reach the background floor; using the nearest higher distinct distances");
  return true;
}

bool thresholds_from_levels(const vec<double>& sample_v,
                            const vec<double>& levels_v,
                            vec<double>& out_v,
                            size_t min_samples,
                            gamma_fit_t* fit_out)
{
  out_v.clear();
  if (levels_v.size() != 4) return false;

  // NOTE (estimator): the fit is method of moments -- closed form, no iteration; its job is only to
  // smooth the background, the robustness comes from the lr/min-portion filter upstream. An MLE or
  // trimmed alternative would live inside fit_gamma_mom.
  // NOTE (floor): floored windows (d below d_eps) are a point mass a continuous gamma cannot
  // represent, and they still take part in this fit, which pulls the shape down and the lower
  // quantiles towards the floor; the lift below catches those.
  gamma_fit_t fit;
  bool lifted = false;
  if (fit_gamma_mom(sample_v, fit) && gamma_thresholds(fit, sample_v, levels_v, out_v, lifted)) {
    if (lifted) warn_msg("--levels: some gamma quantiles reach the background floor; using the nearest observed distance");
    if (fit_out != nullptr) *fit_out = fit;
    return true;
  }

  // Rejected (too few samples, no variance, or quantiles that cannot be placed): fall back to the
  // empirical quantiles of the same sample rather than skipping the sketch.
  out_v.clear();
  if (!empirical_thresholds_from_levels(sample_v, levels_v, out_v, min_samples)) return false;
  if (fit_out != nullptr) *fit_out = gamma_fit_t{}; // NaN marks the empirical path for the caller
  return true;
}

template<typename T>
IntMap<T>::IntMap(const map_opts& opts, const Sketch& sketch, const vec<qseq_t>& batch_v)
  : opts(opts)
  , sketch(sketch)
  , batch_v(batch_v)
  , k(sketch.get_k())
  , hpos(sketch.get_h())
  , sketch_canonical(sketch.is_canonical())
{
  scratch_v.assign(hdist_bound + 1, 0);
  acc_v.assign(hdist_bound + 1, 0);
}

template<typename T>
void IntMap<T>::build_null(ThreadPool& pool)
{
  // One (possibly empty) pool per query sequence, so `null_seq_v[r.bix]` is in range even when the
  // significance test is switched off with --sample-size 0.
  null_seq_v.assign(batch_v.size(), {});
  if (opts.sample_size == 0) return;

  BackgroundSampler sampler(sketch, batch_v, opts.tau, opts.bin_shift, opts.hdist_th);
  const vec<sample_point_t> points = sampler.sample(opts.sample_size, opts.per_sequence, pool);

  null_v.clear();
  for (const sample_point_t& p : points) {
    if (!is_valid_distance(p.d)) continue;
    if (p.d < d_eps) ++null_floored;
    null_seq_v[p.bix].push_back(p.d);
    null_v.push_back(p.d);
  }
  std::sort(null_v.begin(), null_v.end());
  for (auto& v : null_seq_v)
    std::sort(v.begin(), v.end());

  // The --levels gamma is fitted to the lr-filtered sample; the empirical percentile above keeps
  // using the full pool.
  fit_v = lr_filtered_sample(points, opts.lr_th, opts.min_portion);
  std::sort(fit_v.begin(), fit_v.end());
}

template<typename T>
void IntMap<T>::map_sequences(std::ostream& sout, const str& rname, ThreadPool& pool)
{
  build_null(pool);
  // A record is tested only against a pool of at least a tenth of the requested --sample-size.
  const size_t min_samples = std::max<size_t>(1, opts.sample_size / 10);
  if (opts.levels_v.empty()) {
    th_v = sorted_unique(opts.thresholds_v);
    // Also guards direct (non-CLI) construction; the CLI checks this in validate_configuration.
    if (th_v.size() != 1 && th_v.size() != WIDTH)
      error_exit(concat_msg("-d needs exactly 1 or ", WIDTH, " thresholds; got ", th_v.size()));
  } else {
    // Two ways to place the low lanes. The empirical quantiles of the full background are the
    // faithful local estimate and the default; the gamma smooths and extrapolates, but its lower
    // tail can sit tighter than the data supports (see --bg-model).
    if (opts.bg_model == "gamma") {
      gamma_fit_t fit;
      if (!thresholds_from_levels(fit_v, opts.levels_v, th_v, min_samples, &fit)) {
        // Neither the gamma nor the empirical fallback could place WIDTH distinct quantiles for this
        // sketch; skip it rather than aborting the whole run. Other sketches are unaffected.
        warn_pmsg(rname, "skipped: --levels do not resolve into 8 distinct thresholds; the background is too coarse");
        return;
      }
      gamma_fit = fit;
      if (!std::isfinite(gamma_fit.shape))
        warn_pmsg(rname, "--levels: gamma fit rejected; fell back to empirical quantiles");
    } else {
      // The gate below reports and drops any lane that lands on the floor, so the generic lift
      // warning is only useful when the gate is off.
      if (!empirical_thresholds_from_levels(null_v, opts.levels_v, th_v, min_samples, opts.no_lane_gate)) {
        warn_pmsg(rname, "skipped: --levels do not resolve into 8 distinct thresholds; the background is too coarse");
        return;
      }
      gamma_fit = gamma_fit_t{}; // NaN marks the empirical path
    }

    // A threshold outside the likelihood-ratio bounds cannot separate anything: below the floor it
    // is indistinguishable from perfect identity, at or above the ceiling from the weakest
    // detectable match. Drop those lanes rather than lifting them to an observed distance.
    if (!opts.no_lane_gate) {
      const LLH<double> probe = make_llhf(sketch, opts.hdist_th);
      const uint64_t n_win = get_nwinmers(opts.tau, opts.bin_shift);
      ar = admissible_distance_range<double>(probe, n_win, opts.lr_th);
      size_t n_active = 0;
      for (size_t i = 0; i < th_v.size(); ++i) {
        const bool is_ok = ar.is_ok && th_v[i] >= ar.t_min && th_v[i] <= ar.t_max;
        active_v[i] = is_ok;
        if (is_ok) {
          ++n_active;
        } else if (opts.verbosity >= 1) {
          warn_pmsg(rname,
                    "--levels: dropping threshold ",
                    th_v[i],
                    ar.is_ok ? concat_msg(" outside [", ar.t_min, ", ", ar.t_max, "]")
                             : std::string(" (no admissible range for this window length)"));
        }
      }
      if (n_active == 0) {
        warn_pmsg(rname, "skipped: none of the ", th_v.size(), " thresholds survives the likelihood-ratio bounds");
        return;
      }
    }
  }

  const T dist_th = make_dist_th<T>(th_v);
  map_params<T> dparams(dist_th, opts.hdist_th, opts.tau, opts.bin_shift, opts.chisq);
  dparams.active = active_v;
  const LLH<T> llhf(k, hpos, sketch.get_rho(), opts.hdist_th, dist_th);

  for (size_t bix = 0; bix < batch_v.size(); ++bix)
    scan_sequence(dparams, llhf, bix);

  d_acc = llhf.mle(acc_v.data(), u_acc);
  for (record_t& r : records_v) {
    const vec<double>& pool_v = opts.per_sequence ? null_seq_v[r.bix] : null_v;
    apply_empirical_significance(r, pool_v, min_samples);
  }
  benjamini_hochberg_correction(records_v);

  if (opts.verbosity >= 1) report_null(rname);
  report_contiguous(sout, rname, llhf);
}

template<typename T>
void IntMap<T>::scan_sequence(const map_params<T>& params, const LLH<T>& llhf, size_t bix)
{
  const char* cseq = batch_v[bix].seq.data();
  const uint64_t len = batch_v[bix].seq.size();
  if (len < static_cast<uint64_t>(k)) {
    warn_pmsg(batch_v[bix].qid, "skipped: sequence shorter than k-mer length ", "(len=", len, ", k=", k, ")");
    return;
  }

  const uint64_t enmers = len - k + 1;
  const uint64_t nbins = ceil_bins(enmers, params.bin_shift);
  if (nbins < 2) {
    warn_pmsg(
      batch_v[bix].qid, "skipped: fewer than two bins after binning ", "(len=", len, ", bin_size=", params.bin_size, ")");
    return;
  }
  if (params.tau_bin > nbins)
    warn_pmsg(batch_v[bix].qid, "requested minimum exceeds the query; using the full query as the effective minimum");

  // Minimum interval span in bins: extraction keeps a run only when `right - left >= tau_eff`, and
  // the run then covers exactly the bins its score difference covers.
  const uint64_t tau_eff = std::min(params.tau_bin, nbins);
  const size_t srprev = records_v.size();
  const scan_ctx_t ctx = make_scan_ctx(sketch, params.bin_shift, params.hdist_th);

  if (sketch_canonical) {
    IntExt<T> ext(params, llhf, nbins, enmers);
    scan_mers_range<true>(ctx, cseq, 0, enmers, intext_agg_t<T>{ext, nullptr});
    ext.inclusive_scan();
    ext.compute_prefhistsum();

    vec<uint64_t> v_q_v;
    uint64_t u_q = 0, t_q = 0;
    ext.total_histogram(v_q_v, u_q, t_q);
    const double d_q = llhf.mle(v_q_v.data(), u_q);
    ext.set_query_distance(d_q);
    ext.extrema_scan();

    if (opts.iterative_removal)
      extract_iterative_intervals(ext, llhf, false, tau_eff, bix);
    else
      extract_lane_intervals(ext, llhf, false, tau_eff, bix);

    for (size_t ri = srprev; ri < records_v.size(); ++ri) {
      records_v[ri].d_q = d_q;
      records_v[ri].d_diff = nanx();
    }
    add_to_acc(acc_v, u_acc, v_q_v, u_q);
    return;
  }

  IntExt<T> ext_fw(params, llhf, nbins, enmers);
  IntExt<T> ext_rc(params, llhf, nbins, enmers);
  scan_mers_range<false>(ctx, cseq, 0, enmers, intext_agg_t<T>{ext_fw, &ext_rc});
  for (auto* ext : {&ext_fw, &ext_rc}) {
    ext->inclusive_scan();
    ext->compute_prefhistsum();
  }

  vec<uint64_t> v_q_fw_v, v_q_rc;
  uint64_t u_q_fw = 0, u_q_rc = 0, t_q_fw = 0, t_q_rc = 0;
  ext_fw.total_histogram(v_q_fw_v, u_q_fw, t_q_fw);
  ext_rc.total_histogram(v_q_rc, u_q_rc, t_q_rc);
  const double d_q_fw = llhf.mle(v_q_fw_v.data(), u_q_fw);
  const double d_q_rc = llhf.mle(v_q_rc.data(), u_q_rc);
  const double d_diff = strand_diff(d_q_fw, d_q_rc);

  ext_fw.set_query_distance(d_q_fw);
  ext_rc.set_query_distance(d_q_rc);
  for (auto* ext : {&ext_fw, &ext_rc})
    ext->extrema_scan();

  if (opts.iterative_removal) {
    extract_iterative_intervals(ext_fw, llhf, false, tau_eff, bix);
    extract_iterative_intervals(ext_rc, llhf, true, tau_eff, bix);
  } else {
    extract_lane_intervals(ext_fw, llhf, false, tau_eff, bix);
    extract_lane_intervals(ext_rc, llhf, true, tau_eff, bix);
  }

  for (size_t ri = srprev; ri < records_v.size(); ++ri) {
    record_t& r = records_v[ri];
    r.d_q = r.is_rc ? d_q_rc : d_q_fw;
    r.d_diff = d_diff;
  }

  // The lower-distance strand is the reference for the genome-wide accumulator.
  const bool is_rc = d_diff > 0.0;
  add_to_acc(acc_v, u_acc, is_rc ? v_q_rc : v_q_fw_v, is_rc ? u_q_rc : u_q_fw);
}

// Default mode: every lane enumerates its own maximal runs over the whole sequence.
template<typename T>
void IntMap<T>::extract_lane_intervals(IntExt<T>& ext, const LLH<T>& llhf, bool is_rc, uint64_t tau_eff, size_t bix)
{
  const uint64_t nbins = ext.get_nbins();
  // thrank_v holds the lanes that survived the likelihood-ratio bounds, so an inadmissible
  // threshold is never extracted in either mode.
  for (size_t ix : ext.get_thrank_v()) {
    // The extractor's indices are the 0-based bin boundaries of the prefix sums: a run [a, b)
    // covers bins a .. b-1, and `emit_record` takes the same range. Runs start at boundary 1
    // because boundary 0 is the score origin (the record-high search looks at `fdpmax[a - 1]`), so
    // an interval never includes the sequence's first k-mer; the fallback below does cover it.
    // `nbins` is the last boundary.
    ext.extract_intervals_mx(tau_eff, 1, nbins, ix);
    ext.expand_intervals(opts.chisq, ix);
    for (const auto& iv : ext.get_intervals_v(ix))
      emit_record(ext, llhf, bix, iv.a, iv.b, ix, is_rc);
  }
}

// --iterative-removal: lanes are walked from the most extreme threshold inward, and each lane
// extracts only in the gaps the earlier lanes left, so records do not repeat across thresholds.
template<typename T>
void IntMap<T>::extract_iterative_intervals(IntExt<T>& ext, const LLH<T>& llhf, bool is_rc, uint64_t tau_eff, size_t bix)
{
  const uint64_t nbins = ext.get_nbins();
  bp_v.clear();

  const auto& thrank = ext.get_thrank_v();
  if (thrank.empty()) return;

  auto run_range = [&](uint64_t a, uint64_t b, size_t ix) { ext.extract_intervals_mx(tau_eff, a, b, ix); };

  size_t sbprev = 0;
  auto merge_from = [&](size_t start) {
    if (bp_v.size() <= start) return;
    auto nit = std::next(bp_v.begin(), start);
    auto cmp = [](const bp_t& a, const bp_t& b) { return a.a_bin < b.a_bin; };
    std::sort(nit, bp_v.end(), cmp);
    std::inplace_merge(bp_v.begin(), nit, bp_v.end(), cmp);
  };

  for (size_t ix : thrank) {
    merge_from(sbprev);

    // `bp_v` holds half-open 0-based bin ranges [a, b), the same range a record covers. Each gap
    // is the unclaimed span between two claims, and only gaps that can still hold a full interval
    // are searched. `prev` starts at boundary 1 for the same reason the extractor does.
    uint64_t prev = 1;
    for (const auto& s : bp_v) {
      if (s.a_bin >= prev + tau_eff) run_range(prev, s.a_bin, ix);
      prev = std::max(prev, s.b_bin);
    }
    if (nbins >= prev + tau_eff) run_range(prev, nbins, ix);
    ext.expand_intervals(opts.chisq, ix);

    const auto& iv_v = ext.get_intervals_v(ix);
    const size_t nprev = bp_v.size();
    for (const auto& iv : iv_v)
      bp_v.push_back({iv.a, iv.b, ix});
    sbprev = nprev;
  }
  merge_from(sbprev);

  if (bp_v.empty()) {
    if (!ext.is_wskip()) {
      // No interval crossed a threshold; report the full query.
      emit_record(ext, llhf, bix, 0, nbins, no_threshold, is_rc);
    } else {
      // One record per skip-free segment of at least tau_eff bins. `x` counts the 1-based bins the
      // histogram uses, so a segment is the 0-based range [a - 1, x - 1).
      uint64_t a = 1;
      for (uint64_t x = 1; x <= nbins; ++x) {
        if (ext.is_skip(x)) {
          if (x - a >= tau_eff) emit_record(ext, llhf, bix, a - 1, x - 1, no_threshold, is_rc);
          a = x + 1;
        }
      }
      if (nbins + 1 - a >= tau_eff) emit_record(ext, llhf, bix, a - 1, nbins, no_threshold, is_rc);
    }
  } else {
    for (const auto& s : bp_v)
      emit_record(ext, llhf, bix, s.a_bin, s.b_bin, s.ix, is_rc);
  }
}

template<typename T>
void IntMap<T>::emit_record(IntExt<T>& ext,
                            const LLH<T>& llhf,
                            size_t bix,
                            uint64_t a_bin,
                            uint64_t b_bin,
                            size_t th_ix,
                            bool is_rc)
{
  const uint64_t L = batch_v[bix].seq.size();
  const uint64_t enmers = L - k + 1;

  uint64_t u = 0, t = 0;
  // The scored range [a_bin, b_bin) is also the reported range, so the distance, the Fisher
  // information and the likelihood ratios describe exactly the interval the extractor selected.
  ext.extract_histogram(a_bin, b_bin, scratch_v, u, t);
  if (t == 0) ++nunmapped;

  // d / fisher_info / lr_ub for this interval; NaN without hits or a valid MLE.
  const double d = t == 0 ? nanx() : llhf.mle(scratch_v.data(), u);
  double fisher_info = nanx();
  double lr_ub = nanx();
  if (is_valid_distance(d)) {
    const double raw_info = llhf.compute_fisher_info(scratch_v.data(), u, d);
    fisher_info = (std::isfinite(raw_info) && raw_info > 0.0) ? raw_info : nanx();
    lr_ub = compute_lr_ub(llhf, d, t + u);
  }

  // get_coordinates takes the 1-based inclusive bin range covering the same k-mers.
  const interval_t bin_iv{a_bin + 1, b_bin + 1};
  const interval_t seq_iv = get_coordinates(bin_iv, opts.bin_shift, enmers, k);
  records_v.emplace_back(bix, L, seq_iv, is_rc, d, fisher_info, th_ix, lr_ub);
}

template<typename T>
xy_t IntMap<T>::get_distance_bin(const record_t& r, const LLH<T>& llhf) const
{
  xy_t d_range{d_eps, d_ub};
  if (r.th_ix != no_threshold) {
    // The bracket names the band the record was *called* in, so it is derived from the lane's
    // threshold and never from the reported distance, which a merged interval may have drifted
    // away from. A low lane reports [inner, t) and a high lane (t, outer]: with
    // --iterative-removal the inner/outer end is the next threshold of the ladder, because those
    // lanes share the sequence between them, while the default mode has every lane on its own and
    // reports from the domain bound instead.
    const double t_i = lane_at(llhf.get_extrema(), r.th_ix);
    const bool is_low = std::isnan(r.d_q) || t_i <= r.d_q;
    const size_t pos = static_cast<size_t>(std::lower_bound(th_v.begin(), th_v.end(), t_i) - th_v.begin());
    if (is_low) {
      d_range.first = (opts.iterative_removal && pos > 0) ? th_v[pos - 1] : d_eps;
      d_range.second = t_i;
    } else {
      d_range.first = t_i;
      d_range.second = (opts.iterative_removal && pos + 1 < th_v.size()) ? th_v[pos + 1] : d_ub;
    }
  } else {
    const double d = is_valid_distance(r.d) ? r.d : r.d_q;
    if (is_valid_distance(d)) d_range = bracket_distance(d, th_v);
  }
  return d_range;
}

template<typename T>
void IntMap<T>::report_null(const str& rname) const
{
  cerr_msg("[",
           rname,
           "] background windows: n=",
           null_v.size(),
           " (floored=",
           null_floored,
           ") median=",
           linear_quantile(null_v, 0.5));
  if (opts.levels_v.empty()) {
    // -d supplied the thresholds; there is no background model to report.
  } else if (opts.bg_model == "gamma" && std::isfinite(gamma_fit.shape)) {
    cerr_msg("[",
             rname,
             "] --levels gamma fit: shape=",
             gamma_fit.shape,
             " scale=",
             gamma_fit.scale,
             " (n=",
             fit_v.size(),
             ", floored=",
             null_floored,
             ")");
  } else {
    const bool is_gamma_fallback = (opts.bg_model == "gamma");
    cerr_msg("[",
             rname,
             is_gamma_fallback ? "] --levels: gamma fit rejected, empirical quantiles over n="
                               : "] --levels: empirical quantiles of the full pool, n=",
             is_gamma_fallback ? fit_v.size() : null_v.size(),
             " (floored=",
             null_floored,
             ")");
  }
  if (!opts.no_lane_gate && !opts.levels_v.empty()) {
    if (ar.is_ok)
      cerr_msg("[", rname, "] threshold bounds (lr_th=", opts.lr_th, "): [", ar.t_min, ", ", ar.t_max, "]");
    else
      cerr_msg("[", rname, "] threshold bounds (lr_th=", opts.lr_th, "): none resolvable");
  }
  for (size_t i = 0; i < th_v.size(); ++i)
    cerr_msg("[", rname, "] threshold=", th_v[i], active_v[i] ? "" : " (dropped)");
}

template<typename T>
void IntMap<T>::report_contiguous(std::ostream& sout, const str& rname, const LLH<T>& llhf) const
{
  // Writes every record for this sketch, one TSV line each (the name predates the mode split).
  for (const auto& r : records_v) {
    const mask_t mask = (r.th_ix != no_threshold) ? static_cast<mask_t>(1u << r.th_ix) : 0;
    const auto d_range = get_distance_bin(r, llhf);
    const bracket_t d_bin{d_range.first, d_range.second};
    write_map_line(sout, batch_v[r.bix].qid, sketch_canonical, rname, r, mask, d_bin, d_acc);
  }
}

template class IntMap<double>;
template class IntMap<cmlane_t>;

// LLH is used by dist.cpp, intext.cpp and roll.cpp as well, but there is no llh.cpp, so its only
// explicit instantiation lives here and every TU links against these two definitions.
template class LLH<double>;
template class LLH<cmlane_t>;

BackgroundSampler::BackgroundSampler(const Sketch& sketch,
                                     const vec<qseq_t>& batch_v,
                                     uint64_t tau,
                                     uint64_t bin_shift,
                                     uint32_t hdist_th)
  : sketch(sketch)
  , batch_v(batch_v)
  , hdist_th(hdist_th)
  , tau(tau)
  , bin_shift(bin_shift)
  , nwinmers(get_nwinmers(tau, bin_shift))
  , llhf(make_llhf(sketch, hdist_th))
{
  sketch_canonical = sketch.is_canonical();
}

vec<sample_point_t> BackgroundSampler::sample(uint64_t sample_size, bool per_sequence, ThreadPool& pool)
{
  plan(sample_size, per_sequence);
  evaluate(pool);
  return std::move(points_v);
}

void BackgroundSampler::plan(uint64_t sample_size, bool per_sequence)
{
  points_v.clear();
  batches_v.clear();

  auto append = [&](uint64_t bix, const window_plan_t::source_t& src) {
    batch_t b;
    b.bix = bix;
    b.enmers = src.enmers;
    b.point0 = points_v.size();
    b.nsamples = src.starts_v.size();
    for (const uint64_t start : src.starts_v) {
      sample_point_t p;
      p.bix = bix;
      p.start = start;
      points_v.push_back(p);
    }
    batches_v.push_back(b);
  };

  if (!per_sequence) {
    vec<uint64_t> source_lens_v(batch_v.size());
    for (size_t i = 0; i < batch_v.size(); ++i)
      source_lens_v[i] = batch_v[i].seq.size();

    // One global sample, spread over the whole batch.
    const window_plan_t wp = make_window_plan(source_lens_v, sketch.get_k(), tau, bin_shift, sample_size, gen);
    for (const window_plan_t::source_t& src : wp.sources_v)
      append(src.bix, src);
    return;
  }

  vec<uint64_t> source_lens_v(1);
  for (size_t bix = 0; bix < batch_v.size(); ++bix) {
    // A sample per source: one draw each, in source order.
    source_lens_v[0] = batch_v[bix].seq.size();
    const window_plan_t wp = make_window_plan(source_lens_v, sketch.get_k(), tau, bin_shift, sample_size, gen);
    for (const window_plan_t::source_t& src : wp.sources_v)
      append(bix, src);
  }
}

void BackgroundSampler::evaluate(ThreadPool& pool)
{
  const uint32_t nworkers = pool.size();
  const scan_ctx_t ctx = make_scan_ctx(sketch, bin_shift, hdist_th);

  struct task_t
  {
    uint32_t batch;
    uint64_t a, b;
  };
  vec<task_t> tasks_v;
  for (uint32_t bi = 0; bi < batches_v.size(); ++bi) {
    const batch_t& b = batches_v[bi];
    const uint64_t nper = std::max<uint64_t>(4, (b.nsamples + nworkers * 4 - 1) / (nworkers * 4));
    for (uint64_t s0 = 0; s0 < b.nsamples; s0 += nper)
      tasks_v.push_back({bi, s0, std::min(s0 + nper, b.nsamples)});
  }

  pool.parallel_for(tasks_v.size(), 1, [&](uint64_t ti) {
    const task_t& t = tasks_v[ti];
    const batch_t& b = batches_v[t.batch];
    const char* cseq = batch_v[b.bix].seq.data();
    if (sketch_canonical) {
      window_counts_t agg(hdist_th);
      for (uint64_t s = t.a; s < t.b; ++s) {
        sample_point_t& p = points_v[b.point0 + s];
        const uint64_t jx = p.start << bin_shift;
        const uint64_t jy = std::min(jx + nwinmers, b.enmers);
        agg.clear();
        scan_mers_range<true>(ctx, cseq, jx, jy, agg);
        p.d = llhf.mle(agg.hist(), agg.u);
        p.s = compute_lr_ub(llhf, p.d, agg.u + hist_total(agg.hist(), hdist_th));
      }
    } else {
      swindow_counts_t agg(hdist_th);
      for (uint64_t s = t.a; s < t.b; ++s) {
        sample_point_t& p = points_v[b.point0 + s];
        const uint64_t jx = p.start << bin_shift;
        const uint64_t jy = std::min(jx + nwinmers, b.enmers);
        agg.clear();
        scan_mers_range<false>(ctx, cseq, jx, jy, agg);
        const double d_fw = llhf.mle(agg.hist_fw(), agg.u_fw);
        const double d_rc = llhf.mle(agg.hist_rc(), agg.u_rc);
        const auto picked = select_strand_distance(d_fw, d_rc);
        p.d = picked.first;
        const bool is_rc = picked.second == '-';
        p.s = compute_lr_ub(
          llhf, p.d, (is_rc ? agg.u_rc : agg.u_fw) + hist_total(is_rc ? agg.hist_rc() : agg.hist_fw(), hdist_th));
      }
    }
  });
}

bool MapSC::validate_configuration()
{
  bool is_invalid = false;
  const bool have_d = !params.thresholds_v.empty();
  const bool have_levels = !params.levels_v.empty();
  if (have_d == have_levels) {
    is_invalid = true;
    cerr_msg("provide exactly one of -d/--dist-th or --levels");
  }
  if (have_d) {
    if (params.thresholds_v.size() != 1 && params.thresholds_v.size() != rwidth) {
      is_invalid = true;
      cerr_msg("-d accepts exactly 1 or ", rwidth, " thresholds; got ", params.thresholds_v.size());
    }
    for (size_t i = 0; i < params.thresholds_v.size(); ++i) {
      if (params.thresholds_v[i] <= 0.0) {
        is_invalid = true;
        cerr_msg("--dist-th[", i, "] must be positive: ", params.thresholds_v[i]);
      }
    }
    const auto sorted_v = sorted_unique(params.thresholds_v);
    if (sorted_v.size() != params.thresholds_v.size()) {
      is_invalid = true;
      cerr_msg("--dist-th values must be unique");
    }
  }
  if (have_levels) {
    if (params.levels_v.size() != 4) {
      is_invalid = true;
      cerr_msg("--levels requires exactly 4 levels; got ", params.levels_v.size());
    }
    for (size_t i = 0; i < params.levels_v.size(); ++i) {
      if (!std::isfinite(params.levels_v[i]) || params.levels_v[i] <= 1e-8 || params.levels_v[i] >= 0.5) {
        is_invalid = true;
        cerr_msg("--levels[", i, "] must be in (1e-8, 0.5): ", params.levels_v[i]);
      }
    }
    if (sorted_unique(params.levels_v).size() != params.levels_v.size()) {
      is_invalid = true;
      cerr_msg("--levels values must be unique");
    }
    for (const double alpha : params.levels_v) {
      const double expected = static_cast<double>(params.sample_size) * alpha / 2.0;
      if (expected < 1.0) {
        warn_msg(concat_msg("--levels ",
                            alpha,
                            ": only ~",
                            expected,
                            " background windows expected in this tail with --sample-size ",
                            params.sample_size,
                            "; raise --sample-size"));
      }
    }
  }
  if (params.sample_size != 0 && params.sample_size < min_sample_size) {
    is_invalid = true;
    cerr_msg(
      "--sample-size must be 0 (skip the significance test) or at least ", min_sample_size, "; got ", params.sample_size);
  }
  if (params.hdist_th > hdist_bound) {
    is_invalid = true;
    cerr_msg("--hdist-th must be in [0, ", hdist_bound, "] with the current SIMD histogram layout; got ", params.hdist_th);
  }
  if (!validate_binning(params.bin_shift, params.tau)) {
    is_invalid = true;
  }
  const uint64_t bin_size = (params.bin_shift <= 16) ? (uint64_t(1) << params.bin_shift) : 0;
  const uint64_t tau_bin = (bin_size > 0) ? ceil_bins(params.tau, params.bin_shift) : 0;
  if (tau_bin < 2) {
    is_invalid = true;
    cerr_msg("-l must span at least two bins after binning ", "(tau=", params.tau, ", bin_size=", bin_size, ")");
  }
  return !is_invalid;
}

void MapSC::map()
{
  set_precision(*output_stream, 5);
  write_provenance(*output_stream);

  QSeq qs(query_path, std::numeric_limits<uint64_t>::max());

  // read_next_batch appends the last batch before returning false.
  while (qs.read_next_batch()) {
  }
  total_qseq = qs.get_batch_v().size();
  const vec<qseq_t>& batch_v = qs.get_batch_v();

  const Container file(sketch_path);
  const uint32_t nsketches = file.size();
  write_map_header(*output_stream, file.get_config().canonical);

  ThreadPool pool(num_threads);
  cerr_msg("Processing ", nsketches, " sketch(es) w/ ", pool.size(), " thread(s)...");
  // One RNG stream for the whole container: sketch k's null continues where sketch k-1 left off.
  init_thread_rng(1);

  const bool scalar = params.levels_v.empty() && params.thresholds_v.size() == 1;
  uint64_t nunmapped = 0;
  for (uint32_t i = 0; i < nsketches; ++i) {
    const Sketch sketch = file.open(i, SketchLoad::Buckets);

    if (scalar) {
      IntMap<double> intmap(params, sketch, batch_v);
      intmap.map_sequences(*output_stream, sketch.get_rname(), pool);
      nunmapped += intmap.get_nunmapped();
    } else {
      IntMap<cmlane_t> intmap(params, sketch, batch_v);
      intmap.map_sequences(*output_stream, sketch.get_rname(), pool);
      nunmapped += intmap.get_nunmapped();
    }

    const scentry& e = file.get_entry(i);
    file.advise(e.offset, e.len, MADV_DONTNEED);
    if (params.verbosity >= 1) progress("Processed sketch", i + 1, nsketches);
  }
  if (params.verbosity >= 1) progress_done();

  if (nunmapped > 0) cerr_msg("Unmapped intervals (no k-mer hits): ", nunmapped, " (distance reported as NA)");
}

MapSC::MapSC(CLI::App& sc)
{
  // Variable-width vector options keep consuming following tokens while the required positionals
  // still look outnumbered, so `--dist-th 0.05 <query> <ref> -o out` swallowed both paths, and so
  // did `--levels 0.1 0.05 0.01 0.005 <query> <ref>`. Validating the optional extras stops the scan
  // at the first token that is not a number; the per-option checks make that validation reject
  // paths.
  sc.validate_optional_arguments();
  // The loosest level sets how much of the genome is called: 0.1 put the loosest low lane at the
  // 0.05 quantile and covered ~18 % of mouse chr7 per reference. 0.02 covers ~5-7 % and scores a
  // higher F1 against the known chr7 introgression, because the loosest lanes were mostly noise.
  // The tightest level needs ~1/(alpha/2) background windows to be estimable, hence --sample-size.
  params.levels_v = {0.02, 0.01, 0.005, 0.004};
  sc.add_option("query-path", query_path, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible)")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("sketch-path", sketch_path, "Reference container <path>")->required()->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_option("-l", params.tau, "Minimum interval length in k-mers; also the background window length")
    ->required()
    ->check(CLI::PositiveNumber);
  sc.add_option("-b,--bin-shift", params.bin_shift, "Group consecutive k-mers into bins of size 2^b [0]")
    ->check(CLI::Range(0, 16));
  sc.add_option("--hdist-th", params.hdist_th, "Maximum Hamming distance for a k-mer to match [3]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("--chisq", params.chisq, "Chi-square threshold for merging overlapping intervals [33.00051]")
    ->check(CLI::NonNegativeNumber);
  sc.add_option("-d,--dist-th", params.thresholds_v, "Distance threshold(s): exactly 1 or 8 - overrides --levels")
    ->expected(1, rwidth)
    ->check(CLI::PositiveNumber);
  sc.add_option("--levels",
                params.levels_v,
                "Two-sided tail probabilities: exactly 4; thresholds are quantiles of the background "
                "[0.02 0.01 0.005 0.004]")
    ->expected(4)
    ->check(CLI::Range(1e-8, 0.5));
  sc.add_option("--lr-th",
                params.lr_th,
                "A background window is kept for the --levels gamma fit when its lr_ub exceeds this, and a "
                "--levels threshold must be resolvable at this likelihood ratio [6.635]")
    ->check(CLI::NonNegativeNumber);
  sc.add_option("--min-portion",
                params.min_portion,
                "Fit --levels on the lr-filtered background iff at least this portion survives [0.66]")
    ->check(CLI::Range(0.0, 1.0));
  sc.add_option("--bg-model",
                params.bg_model,
                "Where --levels thresholds come from: 'empirical' quantiles of the full background pool, "
                "or a 'gamma' fitted to the lr-filtered pool [empirical]")
    ->check(CLI::IsMember({"empirical", "gamma"}));
  sc.add_option("--sample-size",
                params.sample_size,
                "Background windows sampled per reference; 0 skips the test, otherwise at least 100 [500]")
    ->check(CLI::NonNegativeNumber);
  sc.add_flag(
    "--no-lane-gate", params.no_lane_gate, "Keep --levels thresholds that fall outside the likelihood-ratio bounds [off]");
  sc.add_flag("--iterative-removal",
              params.iterative_removal,
              "Let each threshold claim its intervals and search only the remaining gaps; off, every "
              "threshold enumerates its own intervals [off]");
  sc.add_flag("--per-sequence", params.per_sequence, "Measure significance against a per-query background [per-reference]");
  sc.add_option("--verbosity", params.verbosity, "Background and threshold report detail (0-1) [1]")->check(CLI::Range(0, 1));
  sc.callback([&]() {
    const bool has_d = sc.count("-d") != 0;
    const bool has_levels = sc.count("--levels") != 0;
    if (has_d && has_levels) {
      cerr_msg("-d/--dist-th and --levels are mutually exclusive");
      error_exit("Invalid configuration!");
    }
    if (has_d) params.levels_v.clear();
    if (!validate_configuration()) {
      error_exit("Invalid configuration!");
    }
    open_output(output_file, output_path, output_stream);
  });
}
