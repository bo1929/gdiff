#include "diststat.hpp"
#include "msg.hpp"
#include <algorithm>
#include <cmath>

namespace {
  // Applies a fitted null to one record: the cdf at the observed distance
  // (two-sided for the reference strand), and the fold change vs the latent
  // median. Returns false when the cdf is not finite.
  bool apply_null_score(record_t& r, double d_obs, const GammaModel::params_t& gp, double median)
  {
    const double prob = GammaModel::cdf(d_obs, gp.shape, gp.scale);
    if (!std::isfinite(prob)) return false;
    const bool two_sided = !std::isnan(r.d_diff) && (r.is_rc == (r.d_diff > 0.0));
    r.percentile = two_sided ? (2.0 * std::min(prob, 1.0 - prob)) : prob;
    if (median > eps) r.fold = r.d / median;
    return true;
  }
} // namespace

template<typename T>
DistanceStat<T>::DistanceStat(const params_t<T>& params, const llh_sptr_t<T>& llhf)
  : params(params)
  , llhf(llhf)
{
}

template<typename T>
void DistanceStat<T>::clear_samples()
{
  samples_v.clear();
}

template<typename T>
void DistanceStat<T>::sample_null_pool(const DIM<T>& dim, uint64_t tau_eff, uint64_t bix)
{
  // samples_v is not cleared here: null windows accumulate across queries in the batch.
  const uint64_t nbins_q = dim.get_nbins();
  const uint64_t win_len = tau_eff + 1; // minimum window length in bins
  if (nbins_q < win_len) return;

  const uint64_t max_nwindows = nbins_q / win_len;
  const uint64_t target = std::min<uint64_t>(params.sample_size, max_nwindows);

  std::uniform_int_distribution<uint64_t> rvstart(0, nbins_q - win_len);
  vec<uint64_t> v(hdist_bound + 1);

  // Fixed-size random windows, same procedure as dist/detect. Windows with no
  // k-mer hits (t == 0) are unmapped: skipped and counted, never fitted. To
  // keep the pool at the target size despite unmapped windows, draws continue
  // past the target, up to a bounded budget; if the hit rate is too low to
  // build a reliable null within that budget, the pool stays short and
  // test_significance reports the fit as unreliable.
  const uint64_t max_attempts = std::max<uint64_t>(target * 8, 64);
  uint64_t n_usable = 0;
  for (uint64_t attempt = 0; attempt < max_attempts && n_usable < target; ++attempt) {
    const uint64_t x = rvstart(gen);
    const uint64_t a_bin = x + 1;
    const uint64_t b_bin = x + win_len + 1;

    uint64_t u, t;
    dim.extract_histogram(a_bin - 1, b_bin - 1, v, u, t);
    if (t == 0) {
      ++n_unmapped;
      continue;
    }

    const double d = llhf->mle(v.data(), u);
    if (!std::isfinite(d)) continue;
    const double I = llhf->compute_fisher_info(v.data(), u, d);
    if (!std::isfinite(I) || I <= 0.0) continue;
    samples_v.push_back({d, I, bix, {a_bin, b_bin}});
    ++n_usable;
  }
}

template<typename T>
bool DistanceStat<T>::filter_sample(const record_t& r, vec<p_t>& p_v, uint64_t sample_size) const
{
  bool excluded = false;
  p_v.clear();
  p_v.reserve(samples_v.size());
  for (const auto& s : samples_v) {
    if (s.bix == r.bix && overlaps_half_open(s.bin_iv, r.bin_iv)) {
      excluded = true;
      continue;
    }
    if (!std::isfinite(s.I)) continue;
    const double d = validate_distance(s.d);
    if (!std::isfinite(d)) continue;
    p_v.push_back({d, s.I});
  }

  if (p_v.size() <= sample_size) return excluded;

  size_t w = 0;
  for (size_t i = 0; i < p_v.size(); ++i) {
    if (w < sample_size) {
      p_v[w++] = p_v[i];
    } else {
      const size_t j = std::uniform_int_distribution<size_t>(0, i)(gen);
      if (j < sample_size) p_v[j] = p_v[i];
    }
  }
  p_v.resize(static_cast<size_t>(sample_size));
  return excluded;
}

template<typename T>
bool DistanceStat<T>::test_significance(record_t& r, uint64_t sample_size, const str& qid)
{
  const double d_obs = validate_distance(r.d);
  if (!std::isfinite(d_obs)) return false; // unmapped interval (no k-mer hits): no significance

  vec<p_t> p_v;
  const bool excluded = filter_sample(r, p_v, sample_size);

  if (p_v.size() < GammaModel::min_nsamples) {
    warn_pmsg(qid, "not enough null samples; skipping significance test");
    return false;
  }

  // The gamma fit and latent median depend only on the filtered null set.
  // Records of one query are contiguous and share that set unless an overlap
  // exclusion applied, so the Nelder-Mead fit is computed once per query.
  if (excluded || r.bix != fit_cache_bix) {
    vec<double> d_v;
    d_v.reserve(p_v.size());
    for (const auto& est : p_v) {
      if (const double d = validate_distance(est.d); std::isfinite(d)) {
        d_v.push_back(d);
      }
    }
    const GammaModel::params_t gp = GammaModel::fit_from_samples(d_v);
    const bool ok = GammaModel::validate_params(gp);
    const double median = ok ? GammaModel::median_from_params(gp, d_eps, d_ub - d_eps) : nanx();
    if (!excluded) {
      fit_cache_bix = r.bix;
      fit_cache_params = gp;
      fit_cache_ok = ok;
      fit_cache_median = median;
    }
    if (!ok) {
      warn_pmsg(qid, "gamma fit failed; skipping significance test");
      return false;
    }
    return apply_null_score(r, d_obs, gp, median);
  }

  if (!fit_cache_ok) {
    warn_pmsg(qid, "gamma fit failed; skipping significance test");
    return false;
  }
  return apply_null_score(r, d_obs, fit_cache_params, fit_cache_median);
}

template<typename T>
void DistanceStat<T>::benjamini_hochberg_correction(vec<record_t>& records_v)
{
  for (bool is_rc : {false, true}) {
    vec<size_t> idx;
    for (size_t i = 0; i < records_v.size(); ++i)
      if (records_v[i].is_rc == is_rc && std::isfinite(records_v[i].percentile)) idx.push_back(i);
    if (idx.empty()) continue;

    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return records_v[a].percentile < records_v[b].percentile; });

    const double m = static_cast<double>(idx.size());
    double q_min = 1.0;
    for (size_t rank = idx.size(); rank >= 1; --rank) {
      record_t& r = records_v[idx[rank - 1]];
      q_min = std::min(q_min, std::min(1.0, r.percentile * m / static_cast<double>(rank)));
      r.qvalue = q_min;
    }
  }
}

template class DistanceStat<double>;
template class DistanceStat<cm512_t>;
