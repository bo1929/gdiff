#include "map.hpp"

#include "utils.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "records.hpp"
#include "scan.hpp"
#include "windows.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <exception>
#include <limits>
#include <boost/math/distributions/gamma.hpp>
#include <simde/x86/avx512.h>
#include <sys/mman.h>

namespace {
  vec<double> sorted_unique(vec<double> v)
  {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
  }

  // Threshold values as a SIMD lane vector, padded with d_eps in the unused lanes.
  cmlane_t make_dist_th(const vec<double>& th_v)
  {
    cmlane_t lanes{};
    lanes.fill(d_eps);
    for (size_t i = 0; i < th_v.size(); ++i)
      lanes[i] = th_v[i];
    return lanes;
  }

  // The strand columns exist only for strand-aware sketches; `write_call` keeps this column order.
  inline void write_map_header(std::ostream& sout, bool canonical)
  {
    tsv_row row(sout);
    row("seq")("seq_len")("start")("end");
    if (!canonical) row("strand")("is_rc");
    row("reference")("pi")("threshold")("side")("interval_len")("n_kmers");
    row("d")("d_seq")("d_null")("delta")("d_min")("d_max")("d_floor")("d_ceiling");
    row("n_at_floor")("n_at_ceiling")("p_interval")("q_interval")("p_gamma");
    sout << std::endl;
  }
}

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
    if (fit.p0 > 0.0 && !(p > fit.p0)) return fit.d_floor;
    const boost::math::gamma_distribution<double> gamma(fit.shape, fit.scale);
    return boost::math::quantile(gamma, (p - fit.p0) / (1.0 - fit.p0));
  }

  // Hurdle CDF. Unlike a fit truncated to the admissible range it is a distribution over the whole
  // distance domain, so it is the p-value scale as well as the threshold scale.
  double hurdle_cdf(const gamma_fit_t& fit, double d)
  {
    if (fit.p0 > 0.0 && !(d > fit.d_floor)) return fit.p0;
    const boost::math::gamma_distribution<double> gamma(fit.shape, fit.scale);
    return fit.p0 + ((1.0 - fit.p0) * boost::math::cdf(gamma, d));
  }
} // namespace

// The --levels thresholds with the pi that placed each one: for a level pi, the pi quantile and the
// 1 - pi quantile of the hurdle gamma fitted to the sample. Fails closed: a rejected fit or no
// surviving threshold returns false and the caller skips the reference. A repeated threshold keeps
// the smallest pi (the tightest level). Every threshold that is not placed goes to `report`.
bool level_thresholds(const vec<double>& sample_v,
                      const vec<double>& levels_v,
                      vec<double>& th_v,
                      vec<double>& pi_v,
                      size_t min_samples,
                      gamma_fit_t* fit_out,
                      double d_floor,
                      level_report_t* report)
{
  th_v.clear();
  pi_v.clear();
  level_report_t local;
  level_report_t& rep = report != nullptr ? *report : local;
  rep = {};
  if (levels_v.empty() || sample_v.size() < min_samples) {
    rep.skip = concat_msg("only ", sample_v.size(), " background windows; at least ", min_samples, " needed");
    return false;
  }

  struct point_t
  {
    double t;
    double pi;
  };
  vec<double> levels = levels_v;
  std::sort(levels.begin(), levels.end()); // tightest pi first
  vec<point_t> points_v;

  gamma_fit_t fit;
  // Fail closed: thresholds and p_gamma must come from the same fitted model.
  if (!fit_background_gamma(sample_v, fit, min_samples, d_floor)) {
    const size_t n_above =
      sample_v.size() - static_cast<size_t>(std::llround(fit.p0 * static_cast<double>(sample_v.size())));
    rep.skip = n_above < min_samples ? concat_msg("background fit rejected: ",
                                                  n_above,
                                                  " windows above the floor, at least ",
                                                  min_samples,
                                                  " needed (p0=",
                                                  fit.p0,
                                                  ")")
                                     : str("background fit rejected: degenerate distribution above the floor");
    return false;
  }
  // Gamma tails can extend beyond the distance model's domain. Such thresholds cannot be exceeded
  // by a valid estimate, and evaluating derivatives there can create false calls.
  auto place = [&](double pi, bool low) {
    if (low && fit.p0 > 0.0 && !(pi > fit.p0)) {
      rep.drops_v.push_back({fit.d_floor, pi, low, concat_msg("floor mass p0=", fit.p0, " >= pi")});
      return true;
    }
    double t = nanx();
    try {
      t = hurdle_quantile(fit, low ? pi : 1.0 - pi);
    } catch (const std::exception&) {
      return false;
    }
    if (!std::isfinite(t)) return false;
    if (!is_valid_distance(t)) {
      rep.drops_v.push_back({t, pi, low, "outside the distance domain"});
    } else if (!(t > fit.d_floor * (1.0 + 1e-12))) {
      rep.drops_v.push_back({t, pi, low, "at or below the estimable floor"});
    } else {
      points_v.push_back({t, pi});
    }
    return true;
  };
  for (const double pi : levels) {
    if (!place(pi, true) || !place(pi, false)) {
      rep.skip = "quantile of the background fit is not finite";
      return false;
    }
  }

  std::sort(points_v.begin(), points_v.end(), [](const point_t& a, const point_t& b) { return a.t < b.t; });
  for (const point_t& p : points_v) {
    if (!th_v.empty() && p.t == th_v.back()) continue;
    th_v.push_back(p.t);
    pi_v.push_back(p.pi);
  }
  if (th_v.empty()) {
    rep.skip = "no threshold lies inside the estimable range";
    return false;
  }
  if (fit_out != nullptr) *fit_out = fit;
  return true;
}

void apply_gamma_significance(record_t& r, const gamma_fit_t& fit)
{
  if (!is_valid_distance(r.d) || !(fit.shape > 0.0) || !(fit.scale > 0.0)) return;
  if (!std::isfinite(fit.shape) || !std::isfinite(fit.scale)) return;

  double cdf = nanx();
  try {
    cdf = hurdle_cdf(fit, r.d);
  } catch (const std::exception&) {
    return; // outside the fitted distribution domain: leave the row untested
  }
  if (!std::isfinite(cdf)) return;

  r.p_gamma = std::clamp(2.0 * std::min(cdf, 1.0 - cdf), 0.0, 1.0);
}

bool fit_background_gamma(const vec<double>& sample_v, gamma_fit_t& fit, size_t min_samples, double d_floor)
{
  fit = {};
  if (sample_v.empty()) return false;

  // sample_v is sorted, and windows whose histogram repeats get an identical MLE, so the floor is
  // an exact point mass and `p0` counts it rather than modelling it with a continuous density - the
  // reason a bare gamma fit degenerates on a floored pool. Only the windows above the spike carry
  // the gamma, so `min_samples` counts those.
  fit.d_floor = std::isfinite(d_floor) ? d_floor : sample_v.front();
  const double floor_tol = 1e-12 * std::max(1.0, std::abs(fit.d_floor));
  const size_t n_at =
    static_cast<size_t>(std::upper_bound(sample_v.begin(), sample_v.end(), fit.d_floor + floor_tol) - sample_v.begin());
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

IntMap::IntMap(const map_opts& opts, const Sketch& sketch, const vec<qseq_t>& batch_v)
  : opts(opts)
  , sketch(sketch)
  , batch_v(batch_v)
  , k(sketch.get_k())
  , hpos(sketch.get_h())
  , sketch_canonical(sketch.is_canonical())
{
  scratch_v.assign(hdist_bound + 1, 0);
}

void IntMap::build_null(ThreadPool& pool)
{
  if (opts.sample_size == 0) return;

  BackgroundSampler sampler(sketch, batch_v, opts.tau, opts.bin_shift, opts.hdist_th, spurious_f_v);
  const vec<sample_point_t> points = sampler.sample(opts.sample_size, pool);

  bg.null_v.clear();
  for (const sample_point_t& p : points) {
    if (!is_valid_distance(p.d)) continue;
    bg.null_v.push_back(p.d);
  }
  std::sort(bg.null_v.begin(), bg.null_v.end());
  if (bg.null_v.empty()) return;
  bg.d_null = linear_quantile(bg.null_v, 0.5);
  bg.d_min = bg.null_v.front();
  bg.d_max = bg.null_v.back();
  const LLH<double> probe = make_llhf(sketch, opts.hdist_th, spurious_f_v);
  const uint64_t n_winmers = get_nwinmers(opts.tau, opts.bin_shift);
  bg.d_floor = min_estimable_distance(probe, n_winmers);
  bg.d_ceiling = max_estimable_distance(probe, n_winmers);
}

void IntMap::map_sequences(std::ostream& sout, const str& rname, ThreadPool& pool)
{
  spurious_f_v = sketch.get_spurious_f(opts.hdist_th);
  build_null(pool);
  if (bg.null_v.empty()) {
    report_skip(sout, rname, concat_msg("no valid background windows (sample-size=", opts.sample_size, ")"));
    return;
  }
  if (!place_thresholds(sout, rname)) return;

  const size_t n_active = bg.th_v.size();
  const cmlane_t dist_th = make_dist_th(bg.th_v);
  params_t dparams(dist_th, opts.hdist_th, opts.tau, opts.bin_shift, opts.chisq);
  dparams.n_active = static_cast<uint8_t>(n_active);
  const llh_t llhf(k, hpos, sketch.get_rho(), opts.hdist_th, dist_th, true, spurious_f_v);

  for (size_t bix = 0; bix < batch_v.size(); ++bix)
    scan_sequence(dparams, llhf, bix);

  benjamini_hochberg_intervals(records_v);

  if (verbose) report_null(rname);
  report_calls(sout, rname);
}

// Fills bg.th_v / bg.th_pi_v and the hurdle fit from the null sample, either from --levels or from
// --thresholds, and reports them. False when the reference is skipped.
bool IntMap::place_thresholds(std::ostream& sout, const str& rname)
{
  // The bar is a tenth of the windows that survived scoring, and never below min_fit_windows.
  const size_t min_samples = std::max<size_t>(min_fit_windows, bg.null_v.size() / sample_size_reqf);
  vec<lane_drop_t> drops_v;
  if (!opts.thresholds_v.empty()) {
    // --thresholds: the given distances are the lanes; no detection runs. The null still anchors
    // side/delta/p_interval, and the hurdle gamma is still fitted - but only to score p_gamma, so
    // a rejected fit degrades to NaN instead of skipping the sketch.
    bg.th_v = sorted_unique(opts.thresholds_v);
    bg.th_pi_v.assign(bg.th_v.size(), nanx());
    (void)fit_background_gamma(bg.null_v, bg.gamma_fit, min_samples, bg.d_floor);
  } else {
    level_report_t level_rep;
    const bool placed =
      level_thresholds(bg.null_v, opts.levels_v, bg.th_v, bg.th_pi_v, min_samples, &bg.gamma_fit, bg.d_floor, &level_rep);
    drops_v = std::move(level_rep.drops_v);
    for (const auto& d : drops_v)
      warn_msg(rname, ": pi=", d.pi, ": ", d.low ? "lower" : "upper", " threshold dropped (", d.reason, ")");
    if (!placed) {
      report_skip(sout, rname, level_rep.skip);
      return false;
    }
  }
  warn_unresolvable_lanes(rname);
  if (std::isfinite(bg.gamma_fit.shape)) {
    const size_t n_above =
      bg.null_v.size() - static_cast<size_t>(std::llround(bg.gamma_fit.p0 * static_cast<double>(bg.null_v.size())));
    if (n_above < warn_fit_windows)
      warn_msg(rname,
               ": background fit rests on only ",
               n_above,
               " windows above the floor; thresholds are unstable (raise --sample-size)");
  }
  report_lanes(sout, rname, drops_v);
  const size_t n_active = bg.th_v.size();
  if (n_active == 0 || n_active > rwidth)
    error_exit(concat_msg("the threshold set yields ", n_active, " lanes; at most ", rwidth, " fit the register"));
  return true;
}

void IntMap::scan_sequence(const params_t& params, const llh_t& llhf, size_t bix)
{
  const char* cseq = batch_v[bix].seq.data();
  const uint64_t len = batch_v[bix].seq.size();
  if (len < static_cast<uint64_t>(k)) {
    warn_msg("query ", batch_v[bix].qid, ": skipped, shorter than the k-mer length (len=", len, ", k=", k, ")");
    return;
  }

  const uint64_t enmers = len - k + 1;
  const uint64_t nbins = ceil_bins(enmers, params.bin_shift);
  if (nbins < 2) {
    warn_msg("query ", batch_v[bix].qid, ": skipped, fewer than two bins (len=", len, ", bin size=", params.bin_size, ")");
    return;
  }
  if (params.tau_bin > nbins)
    warn_msg("query ", batch_v[bix].qid, ": -l exceeds the query; using the whole query as the minimum interval");

  // Minimum interval span in bins: extraction keeps a run only when `right - left >= tau_eff`, and
  // the run then covers exactly the bins its score difference covers.
  const uint64_t tau_eff = std::min(params.tau_bin, nbins);
  const size_t first_record = records_v.size();
  const scan_ctx_t ctx = make_scan_ctx(sketch, params.bin_shift, params.hdist_th);

  // strands_v[0] is the forward (or canonical) strand; strand-aware sketches add the reverse one.
  std::deque<ext_t> strands_v;
  strands_v.emplace_back(params, llhf, nbins, enmers);
  if (sketch_canonical) {
    scan_mers_range<true>(ctx, cseq, 0, enmers, intext_agg_t<cmlane_t>{strands_v[0], nullptr});
  } else {
    strands_v.emplace_back(params, llhf, nbins, enmers);
    scan_mers_range<false>(ctx, cseq, 0, enmers, intext_agg_t<cmlane_t>{strands_v[0], &strands_v[1]});
  }

  vec<double> d_seq_v;
  for (ext_t& ext : strands_v) {
    ext.inclusive_scan();
    ext.compute_prefhistsum();
    vec<uint64_t> v_seq_v;
    uint64_t u_seq = 0, t_seq = 0;
    ext.total_histogram(v_seq_v, u_seq, t_seq);
    d_seq_v.push_back(llhf.mle(v_seq_v.data(), u_seq));
    ext.set_seq_distance(d_seq_v.back());
    ext.extrema_scan();
  }

  for (size_t si = 0; si < strands_v.size(); ++si)
    extract_intervals(strands_v[si], llhf, si == 1, tau_eff, bix);

  const double d_diff = sketch_canonical ? nanx() : strand_diff(d_seq_v[0], d_seq_v[1]);
  for (size_t ri = first_record; ri < records_v.size(); ++ri) {
    record_t& r = records_v[ri];
    r.d_seq = d_seq_v[r.is_rc ? 1 : 0];
    r.d_diff = d_diff;
  }
}

void IntMap::extract_intervals(ext_t& ext, const llh_t& llhf, bool is_rc, uint64_t tau_eff, size_t bix)
{
  if (opts.independent_levels) {
    extract_lane_intervals(ext, llhf, is_rc, tau_eff, bix);
    return;
  }
  for (const bool low_side : {true, false})
    extract_iterative_intervals(ext, llhf, is_rc, tau_eff, bix, low_side);
}

// --independent-levels: every lane enumerates its own maximal runs over the whole sequence.
void IntMap::extract_lane_intervals(ext_t& ext, const llh_t& llhf, bool is_rc, uint64_t tau_eff, size_t bix)
{
  const uint64_t nbins = ext.get_nbins();
  // thrank_v lists every active lane, most extreme first (ordered against the sequence's own
  // distance in set_seq_distance); both extraction modes walk the same order.
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

// Default gap-claiming mode: lanes are walked from the most extreme threshold inward, and each
// lane extracts only in the gaps the earlier lanes left, so records do not repeat across
// thresholds.
void IntMap::extract_iterative_intervals(ext_t& ext,
                                            const llh_t& llhf,
                                            bool is_rc,
                                            uint64_t tau_eff,
                                            size_t bix,
                                            bool low_side)
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
    if (ext.is_high_lane(ix) == low_side) continue;
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

  // No fallback records: if no level claims a region, it is simply not reported.
  for (const auto& s : bp_v)
    emit_record(ext, llhf, bix, s.a_bin, s.b_bin, s.ix, is_rc);
}

void IntMap::emit_record(ext_t& ext,
                            const llh_t& llhf,
                            size_t bix,
                            uint64_t a_bin,
                            uint64_t b_bin,
                            size_t th_ix,
                            bool is_rc)
{
  const uint64_t L = batch_v[bix].seq.size();
  const uint64_t enmers = L - k + 1;

  // The scored range [a_bin, b_bin) is also the reported range, so the distance and the p-value
  // describe exactly the interval the extractor selected.
  uint64_t u = 0, t = 0;
  ext.extract_histogram(a_bin, b_bin, scratch_v, u, t);
  if (t == 0) ++nunmapped;

  const double d = t == 0 ? nanx() : llhf.mle(scratch_v.data(), u);
  double p_interval = nanx();
  if (is_valid_distance(d) && is_valid_distance(bg.d_null)) {
    const double nll_bg = llhf.nll(bg.d_null, scratch_v.data(), u);
    const double nll_mle = llhf.nll(d, scratch_v.data(), u);
    p_interval = lr_pvalue(likelihood_ratio_statistic(nll_bg, nll_mle));
  }

  // get_coordinates takes the 1-based inclusive bin range covering the same k-mers.
  const interval_t bin_iv{a_bin + 1, b_bin + 1};
  const interval_t seq_iv = get_coordinates(bin_iv, opts.bin_shift, enmers, k);
  records_v.emplace_back(bix, L, seq_iv, is_rc, d, th_ix);
  record_t& r = records_v.back();
  r.p_interval = p_interval;
  if (std::isfinite(bg.gamma_fit.shape)) apply_gamma_significance(r, bg.gamma_fit);
}

// A skipped reference leaves one comment line in the output beside the stderr warning, so an empty
// result is never mistaken for "no outliers": `# skipped<TAB>reference<TAB>reason`.
void IntMap::report_skip(std::ostream& sout, const str& rname, const str& why) const
{
  warn_msg(rname, ": skipped, ", why);
  sout << "# skipped\t" << rname << '\t' << why << '\n';
}

// The thresholds a reference was scanned with, one comment line each, then the levels that lost a
// threshold: `# threshold<TAB>reference<TAB>side<TAB>distance<TAB>pi` (pi is nan under
// --thresholds) and `# dropped<TAB>reference<TAB>side<TAB>distance<TAB>pi<TAB>reason`.
void IntMap::report_lanes(std::ostream& sout, const str& rname, const vec<lane_drop_t>& drops_v) const
{
  for (size_t i = 0; i < bg.th_v.size(); ++i)
    sout << "# threshold\t" << rname << '\t' << (bg.th_v[i] <= bg.d_null ? "low" : "high") << '\t' << bg.th_v[i] << '\t' << bg.th_pi_v[i]
         << '\n';
  for (const lane_drop_t& d : drops_v)
    sout << "# dropped\t" << rname << '\t' << (d.low ? "low" : "high") << '\t' << d.distance << '\t' << d.pi << '\t'
         << d.reason << '\n';
}

// Warn about kept lanes the estimator cannot reach: at or below the floor, or above the ceiling.
void IntMap::warn_unresolvable_lanes(const str& rname) const
{
  for (const double t : bg.th_v) {
    if (std::isfinite(bg.d_floor) && !(t > bg.d_floor * (1.0 + 1e-12)))
      warn_msg(rname, ": threshold=", t, " is at or below the estimable floor (", bg.d_floor, "); expect no calls");
    else if (std::isfinite(bg.d_ceiling) && t > bg.d_ceiling)
      warn_msg(
        rname, ": threshold=", t, " is above the window ceiling (", bg.d_ceiling, "); only intervals longer than -l can reach it");
  }
}

void IntMap::report_null(const str& rname) const
{
  info_msg("[",
           rname,
           "] null: n=",
           bg.null_v.size(),
           " d_null=",
           bg.d_null,
           " d_min=",
           bg.d_min,
           " d_max=",
           bg.d_max,
           " d_floor=",
           bg.d_floor,
           " d_ceiling=",
           bg.d_ceiling);
  if (std::isfinite(bg.gamma_fit.shape)) {
    info_msg("[",
             rname,
             "] hurdle gamma: shape=",
             bg.gamma_fit.shape,
             " scale=",
             bg.gamma_fit.scale,
             " p0=",
             bg.gamma_fit.p0,
             " (",
             bg.gamma_fit.via_ml ? "ml" : "mom",
             ")");
  } else {
    // Only reachable with --thresholds; the levels path skips the reference when the fit fails.
    info_msg("[", rname, "] hurdle gamma rejected; p_gamma unavailable");
  }
}

void IntMap::report_calls(std::ostream& sout, const str& rname) const
{
  size_t n_calls = 0;
  for (const auto& r : records_v) {
    if (r.th_ix >= bg.th_v.size()) continue;
    write_call(sout, rname, r);
    if (is_valid_distance(r.d)) ++n_calls;
  }
  verbose_msg("[", rname, "] calls: ", n_calls, " / ", records_v.size(), " candidates");
}

// One output row; the column order is the one `write_map_header` declares. An interval with no valid
// distance still reports its coordinates and the sequence-level columns, with NA elsewhere.
void IntMap::write_call(std::ostream& sout, const str& rname, const record_t& r) const
{
  const uint64_t interval_len = r.seq_iv.b - r.seq_iv.a + 1;
  const uint64_t n_kmers = interval_len >= k ? interval_len - k + 1 : 0;

  tsv_row row(sout);
  row(batch_v[r.bix].qid)(r.L)(r.seq_iv.a)(r.seq_iv.b);
  if (!sketch_canonical) row(report_strand(r.is_rc, r.d_diff))(static_cast<uint32_t>(r.is_rc));
  row(rname)(bg.th_pi_v[r.th_ix])(bg.th_v[r.th_ix])(is_high_side(bg.th_v[r.th_ix], r.d_seq) ? "high" : "low");
  row(interval_len)(n_kmers);

  if (!is_valid_distance(r.d)) {
    row("NA")(r.d_seq)(bg.d_null)("NA")(bg.d_min)(bg.d_max)(bg.d_floor)(bg.d_ceiling);
    row(0u)(0u)("NA")("NA")(0.0);
  } else {
    const uint32_t at_floor = std::isfinite(bg.d_floor) && r.d <= bg.d_floor * (1.0 + 1e-9) + 1e-12 ? 1 : 0;
    const uint32_t at_ceiling = std::isfinite(bg.d_ceiling) && r.d >= bg.d_ceiling * (1.0 - 1e-9) ? 1 : 0;
    row(r.d)(r.d_seq)(bg.d_null)(r.d - bg.d_null)(bg.d_min)(bg.d_max)(bg.d_floor)(bg.d_ceiling);
    row(at_floor)(at_ceiling)(r.p_interval)(r.q_interval)(r.p_gamma);
  }
  sout << '\n';
}

// LLH is used by dist.cpp, intext.cpp and roll.cpp as well, but there is no llh.cpp, so its only
// explicit instantiation lives here and every TU links against these two definitions.
template class LLH<double>;
template class LLH<cmlane_t>;

BackgroundSampler::BackgroundSampler(const Sketch& sketch,
                                     const vec<qseq_t>& batch_v,
                                     uint64_t tau,
                                     uint64_t bin_shift,
                                     uint32_t hdist_th,
                                     const vec<double>& spurious_f_v)
  : sketch(sketch)
  , batch_v(batch_v)
  , hdist_th(hdist_th)
  , tau(tau)
  , bin_shift(bin_shift)
  , nwinmers(get_nwinmers(tau, bin_shift))
  , llhf(make_llhf(sketch, hdist_th, spurious_f_v))
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
  bool is_valid = true;
  auto require = [&](bool condition, const auto&... msg) {
    if (condition) return;
    is_valid = false;
    error_msg(msg...);
  };

  is_valid &= check_input_type(query_path, input_type::sequence, "map query-path (first positional)");
  is_valid &= check_input_type(sketch_path, input_type::container, "map sketch-path (second positional)");

  if (!params.thresholds_v.empty()) {
    // --thresholds replaces --levels: the distances are given, so there is no tail budget to
    // check, only the lane count, the distance domain and uniqueness.
    const vec<double>& th_v = params.thresholds_v;
    require(th_v.size() <= max_thresholds, "--thresholds accepts at most ", max_thresholds, " values; got ", th_v.size());
    for (size_t i = 0; i < th_v.size(); ++i)
      require(th_v[i] > LB && th_v[i] < UB, "--thresholds[", i, "] must be in (0, ", UB, "); got ", th_v[i]);
    require(sorted_unique(th_v).size() == th_v.size(), "--thresholds values must be unique");
  } else {
    const vec<double>& levels_v = params.levels_v;
    require(!levels_v.empty(), "--levels is required");
    require(levels_v.size() <= max_levels, "--levels accepts at most ", max_levels, " values; got ", levels_v.size());
    for (size_t i = 0; i < levels_v.size(); ++i)
      require(levels_v[i] > 1e-8 && levels_v[i] < 0.5, "--levels[", i, "] must be in (1e-8, 0.5); got ", levels_v[i]);
    require(sorted_unique(levels_v).size() == levels_v.size(), "--levels values must be unique");
    for (const double pi : levels_v) {
      const double expected = static_cast<double>(params.sample_size) * pi;
      if (expected < 1.0)
        warn_msg("--levels ",
                 pi,
                 ": only ~",
                 expected,
                 " background windows expected in each tail with --sample-size ",
                 params.sample_size,
                 " (raise --sample-size)");
    }
  }

  require(params.sample_size >= min_sample_size, "--sample-size must be at least ", min_sample_size, "; got ", params.sample_size);
  require(params.hdist_th <= hdist_bound, "--hdist-th must be in [0, ", hdist_bound, "]; got ", params.hdist_th);
  is_valid &= validate_binning(params.bin_shift, params.tau);
  const uint64_t bin_size = (params.bin_shift <= 16) ? (uint64_t(1) << params.bin_shift) : 0;
  const uint64_t tau_bin = (bin_size > 0) ? ceil_bins(params.tau, params.bin_shift) : 0;
  if (bin_size > 1 && tau_bin * bin_size != params.tau)
    warn_msg("-l ", params.tau, " is rounded up to ", tau_bin * bin_size, " k-mers: intervals span whole bins of ", bin_size);
  require(tau_bin >= 2, "-l must span at least two bins of ", bin_size, " k-mers; got ", params.tau);
  return is_valid;
}

void MapSC::map()
{
  const Container file(sketch_path);
  const uint32_t nsketches = file.size();

  set_precision(*output_stream, 5);
  write_provenance(*output_stream);

  QSeq qs(query_path, std::numeric_limits<uint64_t>::max());
  while (qs.read_next_batch()) {
  }
  const vec<qseq_t>& batch_v = qs.get_batch_v();

  write_map_header(*output_stream, file.get_config().canonical);

  ThreadPool pool(num_threads);
  info_msg("Mapping ", batch_v.size(), " query sequence(s) against ", nsketches, " sketch(es) with ", pool.size(), " thread(s)");
  init_thread_rng(1);

  uint64_t nunmapped = 0;
  for (uint32_t i = 0; i < nsketches; ++i) {
    const Sketch sketch = file.open(i, SketchLoad::Buckets);
    IntMap intmap(params, sketch, batch_v);
    intmap.map_sequences(*output_stream, sketch.get_rname(), pool);
    nunmapped += intmap.get_nunmapped();

    const scentry& e = file.get_entry(i);
    file.advise(e.offset, e.len, MADV_DONTNEED);
    progress_continue("Mapped sketch", i + 1, nsketches);
  }
  progress_done();

  if (nunmapped > 0) info_msg("Unmapped intervals (no k-mer hits): ", nunmapped, " (distance reported as NA)");
}

MapSC::MapSC(CLI::App& sc)
{
  sc.validate_optional_arguments();
  params.levels_v = {0.01, 0.001, 0.0001};

  sc.add_option("query-path", query_path, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible)")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("sketch-path", sketch_path, "Reference container <path>")->required()->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_option("-l", params.tau, "Window and minimum interval length in k-mers; also the background window length")
    ->required()
    ->check(CLI::PositiveNumber);
  sc.add_option("--hdist-th", params.hdist_th, "Maximum Hamming distance for a k-mer to match [3]");
  sc.add_option("--chisq", params.chisq, "Chi-square threshold for merging adjacent intervals [33.00051]")
    ->check(CLI::NonNegativeNumber);
  // The Number checks on --levels and --thresholds end each list at the first non-number, so the
  // positionals after it are not swallowed; the ranges are validated in validate_configuration.
  auto* levels_opt =
    sc.add_option("--levels",
                  params.levels_v,
                  "One to three tail probabilities pi; each places a lower threshold at the pi quantile and "
                  "an upper one at the 1-pi quantile of the sampled hurdle null [0.01 0.001 0.0001]")
      ->expected(1, max_levels)
      ->check(CLI::Number);
  sc.add_option("--thresholds",
                params.thresholds_v,
                "Distance thresholds in (0, 0.5) used as-is, one lane each; skips the background "
                "threshold detection --levels performs. Mutually exclusive with --levels")
    ->expected(1, max_thresholds)
    ->check(CLI::Number)
    ->excludes(levels_opt);
  sc.add_option("--sample-size", params.sample_size, "Background windows sampled per reference; at least 100 [1000]")
    ->check(CLI::NonNegativeNumber);
  sc.add_option("-b,--bin-shift", params.bin_shift, "Bin width as a power of two (0 = one k-mer per bin)");
  sc.add_flag("--independent-levels",
              params.independent_levels,
              "Let every level enumerate its own intervals; off, levels claim in order of extremity and "
              "later levels fill only the gaps [off]");
  sc.callback([&]() {
    if (!validate_configuration()) exit_failure();
    open_output(output_file, output_path, output_stream);
  });
}
