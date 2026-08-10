#include "gamma.hpp"
#include "dim.hpp"
#include "msg.hpp"

bool test_significance(record_t& r,
                       const vec<sample_t>& bg_samples,
                       const uint64_t sample_size,
                       const str& qid,
                       gamma_fit_t* fit)
{
  if (!is_valid_distance(r.d)) return false;
  const double d_obs = r.d;

  vec<sample_t> filtered;
  const bool excluded = filter_background_samples(bg_samples, r, sample_size, filtered);

  if (filtered.size() < GammaModel::min_nsamples) {
    warn_pmsg(qid, "not enough background samples; skipping significance test");
    return false;
  }

  GammaModel::params_t gp{1.0, 1.0};
  double median = nanx();
  bool ok = false;

  if (fit && !excluded && r.bix == fit->bix && r.nbins == fit->nwin_bins && sample_size == fit->sample_size) {
    gp = fit->params;
    median = fit->median;
    ok = fit->ok;
  } else {
    vec<double> d_raw;
    d_raw.reserve(filtered.size());
    for (const auto& est : filtered)
      d_raw.push_back(est.d);
    const auto prepared = GammaModel::prepare_samples(d_raw, d_eps);
    gp = GammaModel::fit_from_samples(prepared.x);
    ok = GammaModel::validate_params(gp);
    median = ok ? GammaModel::median_from_params(gp, d_eps, d_ub - d_eps) : nanx();
    if (fit && !excluded) {
      fit->bix = r.bix;
      fit->nwin_bins = r.nbins;
      fit->sample_size = sample_size;
      fit->params = gp;
      fit->median = median;
      fit->ok = ok;
    }
  }

  if (!ok) {
    warn_pmsg(qid, "gamma fit failed; skipping significance test");
    return false;
  }

  const double cdf = GammaModel::cdf(d_obs, gp.shape, gp.scale);
  if (!std::isfinite(cdf)) return false;
  const double prob = std::clamp(cdf, 0.0, 1.0);
  const bool two_sided = !std::isnan(r.d_diff) && (r.is_rc == (r.d_diff > 0.0));
  r.percentile = std::clamp(two_sided ? 2.0 * std::min(prob, 1.0 - prob) : prob, 0.0, 1.0);
  if (median > eps) r.fold = d_obs / median;
  return true;
}
