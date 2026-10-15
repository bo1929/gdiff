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
                "d_seq",
                "d_acc",
                "pvalue",
                "fold",
                "qvalue",
                "fisher_info",
                "p_lb",
                "p_ub")
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
                "d_seq",
                "d_diff",
                "d_acc",
                "pvalue",
                "fold",
                "qvalue",
                "fisher_info",
                "p_lb",
                "p_ub")
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
                r.d_seq,
                d_acc,
                r.pvalue,
                r.fold,
                r.qvalue,
                r.fisher_info,
                lr_pvalue(r.lr_lb),
                lr_pvalue(r.lr_ub))
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
                r.d_seq,
                r.d_diff,
                d_acc,
                r.pvalue,
                r.fold,
                r.qvalue,
                r.fisher_info,
                lr_pvalue(r.lr_lb),
                lr_pvalue(r.lr_ub))
        << '\n';
    }
  }
} // namespace

bool fit_gamma_mom(const vec<double>& sample_v, gamma_fit_t& fit, size_t min_samples)
{
  fit = {};
  if (sample_v.size() < min_samples) return false;

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
  struct gamma_suff_t
  {
    double n = 0.0;
    double sum = 0.0;
    double sum_log = 0.0;
  };

  gamma_suff_t gamma_suff(const vec<double>& v)
  {
    gamma_suff_t s;
    s.n = static_cast<double>(v.size());
    for (const double x : v) {
      s.sum += x;
      s.sum_log += std::log(x);
    }
    return s;
  }

  double gamma_loglik(const gamma_suff_t& s, double shape, double scale)
  {
    if (!(shape > 0.0) || !(scale > 0.0)) return ninf();
    return ((shape - 1.0) * s.sum_log) - (s.sum / scale) - (s.n * (std::lgamma(shape) + (shape * std::log(scale))));
  }

  // Shape that maximises the gamma likelihood at the profile scale mean/shape. The MLE equation
  // log(shape) - digamma(shape) = log(mean) - mean(log x) has a monotone left side, so bisection
  // settles it; a non-positive right side has no interior solution and is reported as NaN.
  double gamma_ml_shape(const gamma_suff_t& s)
  {
    const double k = std::log(s.sum / s.n) - (s.sum_log / s.n);
    if (!(k > 0.0) || !std::isfinite(k)) return nanx();
    double lo = 1e-3;
    double hi = 1e6;
    for (int it = 0; it < 60; ++it) {
      const double mid = std::sqrt(lo * hi);
      if ((std::log(mid) - boost::math::digamma(mid)) > k)
        lo = mid;
      else
        hi = mid;
    }
    return std::sqrt(lo * hi);
  }

  // Hurdle quantile: the spike below `p0`, the fitted gamma above it. A fit without a spike - every
  // hand-built or degenerate case - is the plain gamma.
  double hurdle_quantile(const gamma_fit_t& fit, double p)
  {
    if (fit.p0 > 0.0 && !(p > fit.p0)) return fit.d_lo;
    const boost::math::gamma_distribution<double> gamma(fit.shape, fit.scale);
    return boost::math::quantile(gamma, (p - fit.p0) / (1.0 - fit.p0));
  }

  // Hurdle CDF. Unlike a fit truncated to the admissible range it is a distribution over the whole
  // distance domain, so it is the p-value scale as well as the threshold scale.
  double hurdle_cdf(const gamma_fit_t& fit, double d)
  {
    if (fit.p0 > 0.0 && !(d > fit.d_lo)) return fit.p0;
    const boost::math::gamma_distribution<double> gamma(fit.shape, fit.scale);
    return fit.p0 + ((1.0 - fit.p0) * boost::math::cdf(gamma, d));
  }

  // The fitted gamma's two-sided thresholds for every level, tightest level first.
  bool gamma_thresholds(const gamma_fit_t& fit, const vec<double>& levels_v, vec<double>& out_v)
  {
    vec<double> levels = levels_v;
    std::sort(levels.begin(), levels.end()); // tightest alpha first

    for (const double alpha : levels) {
      double t_low = 0.0;
      double t_high = 0.0;
      try {
        t_low = hurdle_quantile(fit, alpha / 2.0);
        t_high = hurdle_quantile(fit, 1.0 - (alpha / 2.0));
      } catch (const std::exception&) {
        return false;
      }
      if (!std::isfinite(t_low) || !std::isfinite(t_high)) return false;
      out_v.push_back(t_low);
      out_v.push_back(t_high);
    }
    std::sort(out_v.begin(), out_v.end());
    return true;
  }
} // namespace

// The background pool's own two-sided thresholds for every level, tightest level first. A pool with
// a point mass at the floor returns equal (or near-equal) values there; placing them collapses
// them.
bool empirical_thresholds_from_levels(const vec<double>& sample_v,
                                      const vec<double>& levels_v,
                                      vec<double>& out_v,
                                      size_t min_samples)
{
  out_v.clear();
  if (levels_v.empty() || sample_v.size() < min_samples) return false;

  vec<double> levels = levels_v;
  std::sort(levels.begin(), levels.end()); // tightest alpha first
  for (const double alpha : levels) {
    out_v.push_back(linear_quantile(sample_v, alpha / 2.0));
    out_v.push_back(linear_quantile(sample_v, 1.0 - (alpha / 2.0)));
  }
  std::sort(out_v.begin(), out_v.end());
  return std::all_of(out_v.begin(), out_v.end(), [](double t) { return std::isfinite(t); });
}

// The --levels lane set: the two likelihood-ratio bounds, plus the model's own thresholds, sorted and
// reduced to the distinct values a lane can carry. Nothing is clamped: a threshold the background
// cannot resolve at `-l` length is still a lane, and the per-record p_lb/p_ub say whether that
// record's own length resolved it. Two lanes with the same threshold are one band, so the later
// duplicate is dropped - that is the only threshold a lane can do without.
size_t compose_lanes(vec<double>& th_v, const a_range_t& ar)
{
  th_v.push_back(ar.t_min);
  th_v.push_back(ar.t_max);
  th_v = sorted_unique(std::move(th_v));
  return th_v.size();
}

void apply_gamma_significance(record_t& r, const gamma_fit_t& fit)
{
  if (!is_valid_distance(r.d) || !(fit.shape > 0.0) || !(fit.scale > 0.0)) return;
  if (!std::isfinite(fit.shape) || !std::isfinite(fit.scale)) return;

  double cdf = nanx();
  double median = nanx();
  try {
    cdf = hurdle_cdf(fit, r.d);
    median = hurdle_quantile(fit, 0.5);
  } catch (const std::exception&) {
    return; // outside the fitted distribution domain: leave the row untested
  }
  if (!std::isfinite(cdf)) return;

  r.pvalue = std::clamp(2.0 * std::min(cdf, 1.0 - cdf), 0.0, 1.0);
  if (std::isfinite(median) && median > eps) r.fold = r.d / median;
}

bool fit_background_gamma(const vec<double>& sample_v, gamma_fit_t& fit, size_t min_samples)
{
  fit = {};
  if (sample_v.empty()) return false;

  // sample_v is sorted, and windows whose histogram repeats get an identical MLE, so the floor is
  // an exact point mass and `p0` counts it rather than modelling it with a continuous density - the
  // reason a bare gamma fit degenerates on a floored pool. Only the windows above the spike carry
  // the gamma, so `min_samples` counts those.
  fit.d_lo = sample_v.front();
  const size_t n_at = static_cast<size_t>(std::upper_bound(sample_v.begin(), sample_v.end(), fit.d_lo) - sample_v.begin());
  fit.p0 = static_cast<double>(n_at) / static_cast<double>(sample_v.size());
  const size_t n_above = sample_v.size() - n_at;
  if (n_above < min_samples) return false;

  const vec<double> cont(sample_v.begin() + static_cast<std::ptrdiff_t>(n_at), sample_v.end());
  gamma_fit_t mom;
  if (!fit_gamma_mom(cont, mom, min_samples)) return false;
  fit.shape = mom.shape;
  fit.scale = mom.scale;

  // ML, keeping the moment fit unless the likelihood improves on it.
  const gamma_suff_t suff = gamma_suff(cont);
  const double shape_ml = gamma_ml_shape(suff);
  if (std::isfinite(shape_ml)) {
    const double scale_ml = (suff.sum / suff.n) / shape_ml;
    if (gamma_loglik(suff, shape_ml, scale_ml) > gamma_loglik(suff, mom.shape, mom.scale)) {
      fit.shape = shape_ml;
      fit.scale = scale_ml;
      fit.via_ml = true;
    }
  }
  return std::isfinite(fit.shape) && std::isfinite(fit.scale) && fit.shape > 0.0 && fit.scale > 0.0;
}

bool thresholds_from_levels(const vec<double>& sample_v,
                            const vec<double>& levels_v,
                            vec<double>& out_v,
                            size_t min_samples,
                            gamma_fit_t* fit_out)
{
  out_v.clear();
  if (levels_v.empty()) return false;

  // The gamma is fitted to the pool as a hurdle (see fit_background_gamma): ML with a moment
  // fallback, and the floor spike carried separately. Thresholds the fit places stay lanes of their
  // own; `compose_lanes` neither clamps them nor invents a replacement.
  gamma_fit_t fit;
  if (fit_background_gamma(sample_v, fit, min_samples) && gamma_thresholds(fit, levels_v, out_v)) {
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
  if (opts.sample_size == 0) return;

  BackgroundSampler sampler(sketch, batch_v, opts.tau, opts.bin_shift, opts.hdist_th);
  const vec<sample_point_t> points = sampler.sample(opts.sample_size, pool);

  null_v.clear();
  for (const sample_point_t& p : points) {
    if (!is_valid_distance(p.d)) continue;
    if (p.d < d_eps) ++null_floored;
    null_v.push_back(p.d);
  }
  std::sort(null_v.begin(), null_v.end());
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
    n_active = th_v.size();
  } else {
    // Two model thresholds per level, plus the two likelihood-ratio bounds. The empirical quantiles
    // of the full background are the faithful local estimate and the default; the gamma smooths and
    // extrapolates, but a point mass at the floor pulls its lower quantiles below anything
    // observable (see --bg-model). Either way the bounds are always lanes: they are the claims a
    // `-l` window can always support, one near perfect identity and one near the no-homology
    // ceiling, so the run still reports both ends when the model has nothing to add.
    const LLH<double> probe = make_llhf(sketch, opts.hdist_th);
    const uint64_t n_win = get_nwinmers(opts.tau, opts.bin_shift);
    ar = admissible_distance_range<double>(probe, n_win, opts.lr_th);
    if (!ar.is_ok) {
      warn_pmsg(
        rname, "skipped: no distance is resolvable for a -l window at --lr-th ", opts.lr_th, " (raise -l or lower --lr-th)");
      return;
    }

    bool have_interior = false;
    if (opts.bg_model == "gamma") {
      gamma_fit_t fit;
      have_interior = thresholds_from_levels(null_v, opts.levels_v, th_v, min_samples, &fit);
      gamma_fit = have_interior ? fit : gamma_fit_t{};
      if (!have_interior)
        warn_pmsg(rname, "--levels: no interior threshold resolves from the background; running the two bounds only");
      else if (!std::isfinite(gamma_fit.shape))
        warn_pmsg(rname, "--levels: gamma fit rejected; fell back to empirical quantiles");
    } else {
      have_interior = empirical_thresholds_from_levels(null_v, opts.levels_v, th_v, min_samples);
      gamma_fit = gamma_fit_t{}; // NaN marks the empirical path
      if (!have_interior)
        warn_pmsg(rname, "--levels: no interior threshold resolves from the background; running the two bounds only");
    }
    if (!have_interior) th_v.clear();

    n_active = compose_lanes(th_v, ar);
    // 2 * levels + 2 bounds fits the lanes; the CLI caps the levels, so this only guards callers
    // that build map_opts directly.
    if (n_active > rwidth)
      error_exit(concat_msg("--levels yields ", n_active, " thresholds; at most ", rwidth, " fit the lanes"));
  }

  const T dist_th = make_dist_th<T>(th_v);
  map_params<T> dparams(dist_th, opts.hdist_th, opts.tau, opts.bin_shift, opts.chisq);
  dparams.n_active = static_cast<uint8_t>(n_active);
  dparams.lo_bound = ar.t_min;
  dparams.hi_bound = ar.t_max;
  const LLH<T> llhf(k, hpos, sketch.get_rho(), opts.hdist_th, dist_th);

  for (size_t bix = 0; bix < batch_v.size(); ++bix)
    scan_sequence(dparams, llhf, bix);

  d_acc = llhf.mle(acc_v.data(), u_acc);
  // One model per run: with an accepted gamma fit the p-value and fold come from that fit - the
  // same model that placed the thresholds - otherwise from the empirical pool.
  const bool gamma_significance =
    opts.bg_model == "gamma" && std::isfinite(gamma_fit.shape) && std::isfinite(gamma_fit.scale);
  for (record_t& r : records_v) {
    if (gamma_significance)
      apply_gamma_significance(r, gamma_fit);
    else
      apply_empirical_significance(r, null_v, min_samples);
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

    vec<uint64_t> v_seq_v;
    uint64_t u_seq = 0, t_seq = 0;
    ext.total_histogram(v_seq_v, u_seq, t_seq);
    const double d_seq = llhf.mle(v_seq_v.data(), u_seq);
    ext.set_seq_distance(d_seq);
    ext.extrema_scan();

    if (opts.iterative_removal)
      extract_iterative_intervals(ext, llhf, false, tau_eff, bix);
    else
      extract_lane_intervals(ext, llhf, false, tau_eff, bix);

    for (size_t ri = srprev; ri < records_v.size(); ++ri) {
      records_v[ri].d_seq = d_seq;
      records_v[ri].d_diff = nanx();
    }
    add_to_acc(acc_v, u_acc, v_seq_v, u_seq);
    return;
  }

  IntExt<T> ext_fw(params, llhf, nbins, enmers);
  IntExt<T> ext_rc(params, llhf, nbins, enmers);
  scan_mers_range<false>(ctx, cseq, 0, enmers, intext_agg_t<T>{ext_fw, &ext_rc});
  for (auto* ext : {&ext_fw, &ext_rc}) {
    ext->inclusive_scan();
    ext->compute_prefhistsum();
  }

  vec<uint64_t> v_seq_fw_v, v_seq_rc;
  uint64_t u_seq_fw = 0, u_seq_rc = 0, t_seq_fw = 0, t_seq_rc = 0;
  ext_fw.total_histogram(v_seq_fw_v, u_seq_fw, t_seq_fw);
  ext_rc.total_histogram(v_seq_rc, u_seq_rc, t_seq_rc);
  const double d_seq_fw = llhf.mle(v_seq_fw_v.data(), u_seq_fw);
  const double d_seq_rc = llhf.mle(v_seq_rc.data(), u_seq_rc);
  const double d_diff = strand_diff(d_seq_fw, d_seq_rc);

  ext_fw.set_seq_distance(d_seq_fw);
  ext_rc.set_seq_distance(d_seq_rc);
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
    r.d_seq = r.is_rc ? d_seq_rc : d_seq_fw;
    r.d_diff = d_diff;
  }

  // The lower-distance strand is the reference for the genome-wide accumulator.
  const bool is_rc = d_diff > 0.0;
  add_to_acc(acc_v, u_acc, is_rc ? v_seq_rc : v_seq_fw_v, is_rc ? u_seq_rc : u_seq_fw);
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
    // covers bins a .. b-1, and `emit_record` takes the same range. Boundary 0 is allowed, so a run
    // can cover the sequence's first k-mer; `nbins` is the last boundary.
    ext.extract_intervals_mx(tau_eff, 0, nbins, ix);
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
    // are searched. `prev` starts at boundary 0, so the first bin is searched too.
    uint64_t prev = 0;
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

  // d / fisher_info and both identifiability contrasts for this interval; NaN without hits or a
  // valid MLE. `lr_lb` asks whether the interval is distinguishable from perfect identity, `lr_ub`
  // whether it is distinguishable from the weakest detectable match.
  const double d = t == 0 ? nanx() : llhf.mle(scratch_v.data(), u);
  double fisher_info = nanx();
  double lr_ub = nanx();
  double lr_lb = nanx();
  if (is_valid_distance(d)) {
    const double raw_info = llhf.compute_fisher_info(scratch_v.data(), u, d);
    fisher_info = (std::isfinite(raw_info) && raw_info > 0.0) ? raw_info : nanx();
    lr_lb = compute_lr_lb(llhf, d, t + u);
    lr_ub = compute_lr_ub(llhf, d, t + u);
  }

  // get_coordinates takes the 1-based inclusive bin range covering the same k-mers.
  const interval_t bin_iv{a_bin + 1, b_bin + 1};
  const interval_t seq_iv = get_coordinates(bin_iv, opts.bin_shift, enmers, k);
  records_v.emplace_back(bix, L, seq_iv, is_rc, d, fisher_info, th_ix, lr_ub, lr_lb);
}

template<typename T>
xy_t IntMap<T>::get_distance_bin(const record_t& r, const LLH<T>& llhf) const
{
  xy_t d_range{d_eps, d_ub};
  if (r.th_ix != no_threshold) {
    // The bracket names the band the record was *called* in, so it is derived from the lane's
    // threshold and never from the reported distance, which a merged interval may have drifted
    // away from. A low lane reports [inner, t) and a high lane (t, outer]: with
    // --iterative-removal the inner/outer end is the neighbouring threshold, because those
    // lanes share the sequence between them, while the default mode has every lane on its own and
    // reports from the domain bound instead.
    const double t_i = lane_at(llhf.get_extrema(), r.th_ix);
    const bool is_low = std::isnan(r.d_seq) || t_i <= r.d_seq;
    const size_t pos = static_cast<size_t>(std::lower_bound(th_v.begin(), th_v.end(), t_i) - th_v.begin());
    if (is_low) {
      d_range.first = (opts.iterative_removal && pos > 0) ? th_v[pos - 1] : d_eps;
      d_range.second = t_i;
    } else {
      d_range.first = t_i;
      d_range.second = (opts.iterative_removal && pos + 1 < th_v.size()) ? th_v[pos + 1] : d_ub;
    }
  } else {
    const double d = is_valid_distance(r.d) ? r.d : r.d_seq;
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
             "] --levels gamma fit (thresholds and p-values): shape=",
             gamma_fit.shape,
             " scale=",
             gamma_fit.scale,
             " p0=",
             gamma_fit.p0,
             " (",
             gamma_fit.via_ml ? "ml" : "mom",
             ", n=",
             null_v.size(),
             ", floored=",
             null_floored,
             ")");
  } else {
    const bool is_gamma_fallback = (opts.bg_model == "gamma");
    cerr_msg("[",
             rname,
             is_gamma_fallback ? "] --levels: gamma fit rejected, empirical quantiles over n="
                               : "] --levels: empirical quantiles (thresholds and p-values) of the pool, n=",
             null_v.size(),
             " (floored=",
             null_floored,
             ")");
  }
  if (!opts.levels_v.empty()) {
    cerr_msg("[",
             rname,
             "] threshold bounds (lr_th=",
             opts.lr_th,
             "): [",
             ar.t_min,
             ", ",
             ar.t_max,
             "] + ",
             n_active > 2 ? n_active - 2 : 0,
             " interior lane(s)");
    if (opts.iterative_removal)
      cerr_msg("[", rname, "] iterative removal: the two bounds claim first, then the interior lanes inward");
  }
  for (size_t i = 0; i < n_active; ++i) {
    const bool is_bound = (th_v[i] == ar.t_min) || (th_v[i] == ar.t_max);
    cerr_msg("[", rname, "] threshold=", th_v[i], is_bound ? " (bound)" : "");
  }
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

vec<sample_point_t> BackgroundSampler::sample(uint64_t sample_size, ThreadPool& pool)
{
  plan(sample_size);
  evaluate(pool);
  return std::move(points_v);
}

void BackgroundSampler::plan(uint64_t sample_size)
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

  vec<uint64_t> source_lens_v(batch_v.size());
  for (size_t i = 0; i < batch_v.size(); ++i)
    source_lens_v[i] = batch_v[i].seq.size();

  // One global sample, spread over the whole batch.
  const window_plan_t wp = make_window_plan(source_lens_v, sketch.get_k(), tau, bin_shift, sample_size, gen);
  for (const window_plan_t::source_t& src : wp.sources_v)
    append(src.bix, src);
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
      }
    }
  });
}

bool MapSC::validate_configuration()
{
  bool is_invalid = false;
  if (!check_input_type(query_path, input_type::sequence, "map query-path (first positional)")) is_invalid = true;
  if (!check_input_type(sketch_path, input_type::container, "map sketch-path (second positional)")) is_invalid = true;
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
    if (params.levels_v.size() > max_levels) {
      is_invalid = true;
      cerr_msg("--levels accepts at most ",
               max_levels,
               " levels (two thresholds each, plus the two bounds); got ",
               params.levels_v.size());
    }
    // --levels thresholds are quantiles of the background sample, so a run without a pool has
    // nothing to place: refuse it here instead of skipping every reference at run time.
    if (params.sample_size == 0) {
      is_invalid = true;
      cerr_msg("--levels needs a background sample; --sample-size must be at least ",
               min_sample_size,
               " with --levels (use -d for a pool-free run)");
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
      if (params.sample_size != 0 && expected < 1.0) {
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
  // Intervals are bin-aligned, so a minimum that is not a whole number of bins is rounded up.
  if (bin_size > 1 && tau_bin * bin_size != params.tau)
    warn_msg(concat_msg(
      "-l ", params.tau, " is rounded up to ", tau_bin * bin_size, " k-mers: intervals span whole bins of ", bin_size));
  if (tau_bin < 2) {
    is_invalid = true;
    cerr_msg("-l must span at least two bins after binning ", "(tau=", params.tau, ", bin_size=", bin_size, ")");
  }
  return !is_invalid;
}

void MapSC::map()
{
  // The reference is opened first: a --levels run places its thresholds inside a range a `-l` window
  // can resolve, and a configuration that cannot resolve any distance (too few k-mers per window for
  // the --lr-th cut) is a mistake to report once, before anything is read or written.
  const Container file(sketch_path);
  const uint32_t nsketches = file.size();
  if (!params.levels_v.empty() && nsketches > 0) {
    const Sketch probe = file.open(0, SketchLoad::Buckets);
    const a_range_t probe_ar = admissible_distance_range<double>(
      make_llhf(probe, params.hdist_th), get_nwinmers(params.tau, params.bin_shift), params.lr_th);
    if (!probe_ar.is_ok) {
      error_exit(concat_msg(
        "--levels: no distance is resolvable for a -l window at --lr-th ", params.lr_th, " (raise -l or lower --lr-th)"));
    }
    cerr_msg("--levels: bounds [", probe_ar.t_min, ", ", probe_ar.t_max, "] for a -l window at --lr-th ", params.lr_th);
  }

  set_precision(*output_stream, 5);
  write_provenance(*output_stream);

  QSeq qs(query_path, std::numeric_limits<uint64_t>::max());

  // read_next_batch appends the last batch before returning false.
  while (qs.read_next_batch()) {
  }
  total_qseq = qs.get_batch_v().size();
  const vec<qseq_t>& batch_v = qs.get_batch_v();

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
  // The levels partition what a `-l` window can resolve. The two bounds of that range are always
  // lanes, so a run always reports a near-identical and a near-ceiling claim; the levels fill in
  // between. 0.05 / 0.01 / 0.005 put the loosest low lane at the 0.025 quantile, which covered
  // ~5-7 % of mouse chr7 per reference and scored a higher F1 against the known introgression than
  // the 0.1 ladder, whose loosest lanes were mostly noise. The tightest level needs ~1/(alpha/2)
  // background windows to be estimable, hence --sample-size.
  params.levels_v = {0.05, 0.01, 0.005};
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
                "Two-sided tail probabilities: 1 to 3 (each adds a low and a high threshold); the two "
                "bounds a -l window can resolve are always lanes as well [0.05 0.01 0.005]")
    ->expected(1, max_levels)
    ->check(CLI::Range(1e-8, 0.5));
  sc.add_option("--lr-th",
                params.lr_th,
                "Where the two --levels bounds sit: a -l window must clear this likelihood ratio to be "
                "resolvable from perfect identity and from the sketch ceiling [6.635]")
    ->check(CLI::NonNegativeNumber);
  sc.add_option("--bg-model",
                params.bg_model,
                "Where --levels thresholds come from: 'empirical' quantiles of the background pool, or a "
                "'gamma' fitted to that same pool, which also supplies pvalue and fold [empirical]")
    ->check(CLI::IsMember({"empirical", "gamma"}));
  sc.add_option("--sample-size",
                params.sample_size,
                "Background windows sampled per reference; 0 skips the test, otherwise at least 100 [500]")
    ->check(CLI::NonNegativeNumber);
  sc.add_flag("--iterative-removal",
              params.iterative_removal,
              "Let each threshold claim its intervals and search only the remaining gaps, the two --levels "
              "bounds first; off, every threshold enumerates its own intervals [off]");
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
