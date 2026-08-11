#ifndef _LLH_HPP
#define _LLH_HPP

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>
#include <boost/math/tools/minima.hpp>
#include "stils.hpp"

template<typename T>
class LLH
{
  static constexpr size_t WIDTH = std::is_same_v<T, double> ? 1 : RWIDTH;
  static_assert(std::is_same_v<T, double> || std::is_same_v<T, cm512_t>, "LLH supports only double or cm512_t");

public:
  const uint32_t k;
  const uint32_t h;
  const double rho;
  const uint32_t hdist_th;
  const T extrema;
  std::vector<uint64_t> binom_coef_k;
  std::vector<uint64_t> binom_coef_hnk;
  std::vector<double> binom_coef_k_d;
  std::vector<double> binom_coef_hnk_d;

  LLH(uint32_t k, uint32_t h, double rho, uint32_t hdist_th, T extrema, bool compute_derivatives = true)
    : k(k)
    , h(h)
    , rho(rho)
    , hdist_th(hdist_th)
    , extrema(extrema)
    , binom_coef_k(k + 1)
    , binom_coef_hnk(hdist_th + 1)
    , fdc_v(compute_derivatives ? hdist_th + 1 : 0)
    , sdc_v(compute_derivatives ? hdist_th + 1 : 0)
  {
    // Binomial coefficients for the likelihood model.
    const uint32_t nh = k - h;
    binom_coef_k[0] = 1;
    for (uint32_t d = 0; d < k; ++d)
      binom_coef_k[d + 1] = (binom_coef_k[d] * (k - d)) / (d + 1);
    binom_coef_hnk[0] = 0;
    uint64_t vc = 1;
    for (uint32_t d = 1; d <= hdist_th; ++d) {
      vc = (vc * (nh - d + 1)) / d;
      binom_coef_hnk[d] = binom_coef_k[d] - vc;
    }
    binom_coef_k_d.assign(binom_coef_k.begin(), binom_coef_k.end());
    binom_coef_hnk_d.assign(binom_coef_hnk.begin(), binom_coef_hnk.end());

    if (!compute_derivatives) return;

    if constexpr (std::is_same_v<T, double>) {
      sign = extrema < 0 ? -1.0 : 1.0;
      const double axtrema = extrema * sign;
      for (uint32_t d = 0; d <= hdist_th; ++d) {
        fdc_v[d] = sign * compute_fdc_v(axtrema, d);
        sdc_v[d] = compute_sdc_v(axtrema, d);
      }
      fdc_u = sign * compute_fdc_u(axtrema);
      sdc_u = compute_sdc_u(axtrema);
    } else if constexpr (std::is_same_v<T, cm512_t>) {
      alignas(64) double axtrema[WIDTH];
      for (uint32_t i = 0; i < WIDTH; ++i) {
        sign[i] = extrema[i] < 0 ? -1.0 : 1.0;
        axtrema[i] = extrema[i] * sign[i];
      }
      for (uint32_t d = 0; d <= hdist_th; ++d) {
        for (uint32_t i = 0; i < WIDTH; ++i) {
          fdc_v[d][i] = sign[i] * compute_fdc_v(axtrema[i], d);
          sdc_v[d][i] = compute_sdc_v(axtrema[i], d);
        }
      }
      for (uint32_t i = 0; i < WIDTH; ++i) {
        fdc_u[i] = sign[i] * compute_fdc_u(axtrema[i]);
        sdc_u[i] = compute_sdc_u(axtrema[i]);
      }
    } else {
      static_assert(std::is_same_v<T, double> || std::is_same_v<T, cm512_t>, "LLH supports only double or cm512_t");
    }
  }

  // Passing a reference (&) fails due to SIMD for these
  const T get_sdc(uint32_t d) const { return sdc_v[d]; }

  const T get_fdc(uint32_t d) const { return fdc_v[d]; }

  const T get_sdc() const { return sdc_u; }

  const T get_fdc() const { return fdc_u; }

  T get_sign() const { return sign; }

  T get_extrema() const { return extrema; }

  std::pair<uint8_t, uint8_t> get_sign_bv() const
  {
    uint8_t pos_bv = 0, neg_bv = 0;
    if constexpr (std::is_same_v<T, double>) {
      if (sign > 0)
        pos_bv = 1;
      else
        neg_bv = 1;
    } else {
      for (size_t ti = 0; ti < WIDTH; ++ti) {
        if (sign[ti] > 0)
          pos_bv |= static_cast<uint8_t>(1u << ti);
        else
          neg_bv |= static_cast<uint8_t>(1u << ti);
      }
    }
    return {pos_bv, neg_bv};
  }

  double prob_elude(uint32_t d) const
  {
    return static_cast<double>(binom_coef_hnk[d]) / static_cast<double>(binom_coef_k[d]);
  }

  double prob_collide(uint32_t d) const
  {
    return 1.0 - (static_cast<double>(binom_coef_hnk[d]) / static_cast<double>(binom_coef_k[d]));
  }

  double prob_mutate(double D, uint32_t d) const { return std::pow(1.0 - D, k - d) * std::pow(D, d) * binom_coef_k[d]; }

  double prob_miss(double D) const
  {
    double p = 0;
    for (uint32_t d = 0; d <= hdist_th; ++d) {
      p += prob_elude(d) * prob_mutate(D, d);
    }
    for (uint32_t d = hdist_th + 1; d <= k; ++d) {
      p += prob_mutate(D, d);
    }
    return (rho * p) + 1.0 - rho;
  }

  double prob_hit(double D, uint32_t d) const { return rho * prob_collide(d) * prob_mutate(D, d); }

  double nll(const double& D, const uint64_t* vv, const uint64_t uu) const
  {
    if (D == 0.0) {
      for (uint32_t d = 1; d <= hdist_th; ++d)
        if (vv[d] > 0) return pinf();
    }

    double lsum = 0.0;
    double lv_m = 0.0;
    double powdc = std::pow(1.0 - D, k);
    const double logdn = k * std::log(1.0 - D);
    const double logdp = D > 0.0 ? std::log(D) - std::log(1.0 - D) : 0.0;
    const double ratioD = D / (1.0 - D);

    for (uint32_t d = 0; d <= k; ++d) {
      if (d <= hdist_th) {
        lsum -= (logdn + (d * logdp)) * vv[d];
        lv_m += binom_coef_hnk_d[d] * powdc;
      } else {
        lv_m += powdc * binom_coef_k_d[d];
      }
      powdc *= ratioD;
    }

    return lsum - (std::log((rho * lv_m) + 1.0 - rho) * uu);
  }

  // Analytic observed Fisher information:
  // I(D) = -d^2/dD^2 log L(D) (the negative log-likelihood evaluated at the given D)
  double compute_fisher_info(const uint64_t* v_r, uint64_t u_r, double D) const
  {
    double ll_dd = 0.0;
    for (uint32_t d = 0; d <= hdist_th; ++d) {
      ll_dd += static_cast<double>(v_r[d]) * compute_sdc_v(D, d);
    }
    ll_dd += static_cast<double>(u_r) * compute_sdc_u(D);
    return -ll_dd;
  }

  // Returns NaN when there are no k-mer hits.
  // The distance is undefined, instead of relying on a plateau MLE.
  double mle(const uint64_t* v_r, uint64_t u_r, double* nll_min = nullptr) const
  {
    uint64_t t = 0;
    for (uint32_t d = 0; d <= hdist_th; ++d)
      t += v_r[d];
    if (t == 0) {
      if (nll_min) *nll_min = nanx();
      return nanx();
    }
    auto f = [&](const double& D) { return nll(D, v_r, u_r); };
    xy_t result = boost::math::tools::brent_find_minima(f, LB, UB, 24);
    if (nll_min) *nll_min = result.second;
    return validate_distance(result.first);
  }

private:
  double compute_fdc_v(const double D, const uint32_t d) const
  {
    const double S = 1.0 - D;
    return (d - (k * D)) / (D * S);
  }

  double compute_sdc_v(const double D, const uint32_t d) const
  {
    const double S = 1.0 - D;
    const double numerator = (d * ((2 * D) - 1.0)) - (k * D * D);
    const double denominator = D * D * S * S;
    return numerator / denominator;
  }

  double compute_fdc_u(const double D) const
  {
    const double S = 1.0 - D;
    double gd = 0;
    double fd = 0;
    for (uint32_t d = 0; d <= k; ++d) {
      const double pd = (d - (k * D)) / (D * S);
      const double pe = std::pow(S, k - d) * std::pow(D, d);
      double wt = pe;
      if (d <= hdist_th) {
        wt *= binom_coef_hnk[d];
      } else {
        wt *= binom_coef_k[d];
      }
      gd += wt * pd;
      fd += wt;
    }
    return rho * gd / (1.0 - rho + (rho * fd));
  }

  double compute_sdc_u(const double D) const
  {
    const double S = 1.0 - D;
    double gd = 0;
    double fd = 0;
    double fpd = 0;
    double gpd = 0;
    const double D_sq = D * D;
    const double S_sq = S * S;
    const double denom = D_sq * S_sq;
    for (uint32_t d = 0; d <= k; ++d) {
      const double pd = (d - (k * D)) / (D * S);
      const double pe = std::pow(S, k - d) * std::pow(D, d);
      const double vy = ((d * d) + ((k - 1) * k * D_sq) - (d * (1 + ((k - 1) * 2 * D)))) / denom;
      double wt = pe;
      if (d <= hdist_th) {
        wt *= binom_coef_hnk[d];
      } else {
        wt *= binom_coef_k[d];
      }
      gd += wt * pd;
      fd += wt;
      fpd += wt * pd;
      gpd += wt * vy;
    }
    gd *= rho;
    fpd *= rho;
    gpd *= rho;
    fd = 1.0 - rho + (rho * fd);
    return ((fd * gpd) - (gd * fpd)) / (fd * fd);
  }

  T sign;
  T fdc_u;
  T sdc_u;
  std::vector<T> fdc_v;
  std::vector<T> sdc_v;
};

// -2 log[L(D_ref) / L(d_mle)] from negative log-likelihoods.
inline double likelihood_ratio_statistic(const double nll_ref, const double nll_mle) noexcept
{
  if (!std::isfinite(nll_ref) || !std::isfinite(nll_mle)) return nanx();
  return 2.0 * std::max(0.0, nll_ref - nll_mle);
}

// Weakest detectable match pattern: one hit at hdist_th, all other observed
// k-mers are misses. Its MLE is the largest distance the model can estimate for
// a sketch given n_total observed k-mers.
template<typename T>
inline double max_estimable_distance(const LLH<T>& llhf, uint64_t n_total, double* nll_min = nullptr)
{
  if (n_total == 0) {
    if (nll_min) *nll_min = nanx();
    return nanx();
  }
  arr<uint64_t, hdist_bound + 1> v{};
  v[llhf.hdist_th] = 1;
  return llhf.mle(v.data(), n_total - 1, nll_min);
}

// Likelihood-ratio of the window MLE distance vs the sketch's max estimable
// distance, evaluated on the extreme match counts that define that ceiling.
// Large values mean d is well below the detection limit; near 0 means d is
// indistinguishable from the weakest detectable homology.
template<typename T>
inline double compute_lr_ub(const LLH<T>& llhf, double d, uint64_t n_total)
{
  if (!is_valid_distance(d) || n_total == 0) return nanx();
  arr<uint64_t, hdist_bound + 1> v{};
  v[llhf.hdist_th] = 1;
  const uint64_t u = n_total - 1;
  double nll_ub = nanx();
  const double d_ub_est = llhf.mle(v.data(), u, &nll_ub);
  if (!is_valid_distance(d_ub_est)) return nanx();
  return likelihood_ratio_statistic(llhf.nll(d, v.data(), u), nll_ub);
}

struct likelihood_estimate_t
{
  double d = nanx();
  double I = nanx();
  double lr_bg = nanx();
  double lr_ub = nanx();
  bool has_hits = false;
};

template<typename T>
inline likelihood_estimate_t
compute_likelihood_estimate(const LLH<T>& llhf, const uint64_t* v, uint64_t u, uint64_t t, double d_bg)
{
  likelihood_estimate_t est;
  if (t == 0) return est;
  est.has_hits = true;

  double nll = nanx();
  est.d = llhf.mle(v, u, &nll);
  if (!is_valid_distance(est.d)) return est;

  const double I = llhf.compute_fisher_info(v, u, est.d);
  est.I = (std::isfinite(I) && I > 0.0) ? I : nanx();
  if (is_valid_distance(d_bg)) est.lr_bg = likelihood_ratio_statistic(llhf.nll(d_bg, v, u), nll);
  est.lr_ub = compute_lr_ub(llhf, est.d, t + u);
  return est;
}

template<typename T>
using llh_sptr_t = std::shared_ptr<LLH<T>>;

#endif
