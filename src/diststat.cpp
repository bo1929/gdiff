#include "diststat.hpp"
#include "msg.hpp"
#include <algorithm>
#include <cmath>

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

  const uint64_t max_nconcat = 10;
  const uint64_t max_len = std::min<uint64_t>(max_nconcat * win_len, nbins_q);
  const uint64_t max_nwindows = nbins_q / win_len;
  const uint64_t nsamples = std::min<uint64_t>(params.sample_size, max_nwindows);

  std::uniform_int_distribution<uint64_t> rvstart(0, nbins_q - win_len);
  vec<uint64_t> v(hdist_bound + 1);

  for (uint64_t s = 0; s < nsamples; ++s) {
    const uint64_t x = rvstart(gen);
    uint64_t a_bin = x + 1;
    uint64_t b_bin = x + win_len + 1;

    uint64_t u, t;
    dim.extract_histogram(a_bin - 1, b_bin - 1, v, u, t);

    bool right = true;
    while (t == 0 && (b_bin - a_bin) < max_len) {
      const uint64_t exlen = std::min(win_len, max_len - (b_bin - a_bin));
      if (exlen == 0) break;
      bool concat = false;
      if (right && b_bin + exlen <= nbins_q + 1) {
        b_bin += exlen; // extend right
        concat = true;
      } else if (a_bin > exlen) {
        a_bin -= exlen; // else extend left
        concat = true;
      } else if (b_bin + exlen <= nbins_q + 1) {
        b_bin += exlen; // if blocked, try right regardless
        concat = true;
      }
      if (!concat) break; // hit a boundary

      dim.extract_histogram(a_bin - 1, b_bin - 1, v, u, t);
      right = !right; // alternate sides so growth remains symmetric
    }
    if (t == 0) continue;

    double d = llhf->mle(v.data(), u);
    double I = llhf->compute_fisher_info(d);
    d = validate_distance(d);
    if (std::isfinite(d) && std::isfinite(I)) samples_v.push_back({d, I, bix, {a_bin, b_bin}});
  }
}

template<typename T>
double DistanceStat<T>::sample_box_muller(std::mt19937& rng)
{
  std::uniform_real_distribution<double> U(0.0, 1.0);
  const double u1 = U(rng);
  const double u2 = U(rng);
  const double r = std::sqrt(-2.0 * std::log(u1));
  const double phi = 2.0 * M_PI * u2;
  return r * std::cos(phi);
}

template<typename T>
bool DistanceStat<T>::sample_metropolis_hastings(const p_t& init, vec<p_t>& p_v)
{
  // Metropolis-Hastings posterior draws kept (S) and burn-in iterations (B) per null window.
  bool is_valid = true;
  auto f = [&](const double& D) { return (*llhf)(D); };
  double d = init.d;
  const double step = init.I;
  double nll = f(d);

  std::uniform_real_distribution<double> ruv(0, 1);

  p_v.reserve(p_v.size() + S);

  for (uint64_t iter = 0; iter < B + S; ++iter) {
    double d_i = d + step * sample_box_muller(gen);
    if (d_i <= 0.0) d_i = -d_i;
    if (d_i >= d_ub) d_i = 2.0 * d_ub - d_i;
    d_i = std::clamp(d_i, eps, d_ub - eps);

    const double nll_i = f(d_i);
    const double log_alpha = nll - nll_i;
    if (std::log(ruv(gen)) < log_alpha) {
      d = d_i;
      nll = nll_i;
    }
    if (iter >= B) {
      const double I = llhf->compute_fisher_info(d);
      if ((std::isfinite(I) && std::isfinite(d)) && (d > 0.0 && I > 0.0)) {
        p_v.push_back({d, I});
      } else {
        is_valid = false;
      }
    }
  }
  return is_valid;
}

template<typename T>
void DistanceStat<T>::filter_sample(const record_t& r, vec<p_t>& p_v, uint64_t sample_size) const
{
  p_v.clear();
  p_v.reserve(samples_v.size());
  for (const auto& s : samples_v) {
    if (s.bix == r.bix && overlaps_half_open(s.bin_iv, r.bin_iv)) continue;
    if (!std::isfinite(s.I)) continue;
    const double d = validate_distance(s.d);
    if (!std::isfinite(d)) continue;
    p_v.push_back({d, s.I});
  }

  if (p_v.size() <= sample_size) return;

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
}

template<typename T>
bool DistanceStat<T>::test_significance(record_t& r, uint64_t sample_size, const str& qid)
{
  vec<p_t> p_v;
  filter_sample(r, p_v, sample_size);

  if (p_v.size() < GammaModel::min_nsamples) {
    warn_pmsg(qid, "not enough null samples; skipping significance test");
    return false;
  }

  const double d_obs = validate_distance(r.d);
  vec<double> d_v;
  d_v.reserve(p_v.size());
  for (const auto& est : p_v) {
    if (const double d = validate_distance(est.d); std::isfinite(d)) {
      d_v.push_back(d);
    }
  }
  auto [prob, median] = GammaModel::score_from_samples(d_obs, d_v, d_eps, d_ub - d_eps);
  if (std::isfinite(median)) median = validate_distance(median);

  if (!std::isfinite(prob) || !std::isfinite(median)) {
    warn_pmsg(qid, "gamma fit failed; skipping significance test");
    return false;
  }
  const bool two_sided = !std::isnan(r.d_diff) && (r.is_rc == (r.d_diff > 0.0));
  r.percentile = two_sided ? (2.0 * std::min(prob, 1.0 - prob)) : prob;
  if (std::isfinite(r.d) && median > eps) r.fold = r.d / median;
  return true;
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
