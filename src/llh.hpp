#ifndef _LLH_HPP
#define _LLH_HPP

#include <cmath>
#include <cstdint>
#include <boost/math/tools/minima.hpp>
#include "distance.hpp"

template<typename T>
class LLH
{
  static constexpr size_t WIDTH = std::is_same_v<T, double> ? 1 : rwidth;
  static_assert(std::is_same_v<T, double> || std::is_same_v<T, cmlane_t>, "LLH supports only double or cmlane_t");

public:
  const uint32_t k;
  const uint32_t h;
  const double rho;
  const uint32_t hdist_th;
  const T extrema;
  vec<uint64_t> binom_coef_k;
  vec<uint64_t> binom_coef_hnk;
  vec<double> dcoef_k;
  vec<double> dcoef_hnk;

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
    for (uint32_t x = 0; x < k; ++x)
      binom_coef_k[x + 1] = (binom_coef_k[x] * (k - x)) / (x + 1);
    binom_coef_hnk[0] = 0;
    uint64_t vc = 1;
    for (uint32_t x = 1; x <= hdist_th; ++x) {
      vc = (vc * (nh - x + 1)) / x;
      binom_coef_hnk[x] = binom_coef_k[x] - vc;
    }
    dcoef_k.assign(binom_coef_k.begin(), binom_coef_k.end());
    dcoef_hnk.assign(binom_coef_hnk.begin(), binom_coef_hnk.end());

    if (!compute_derivatives) return;

    if constexpr (std::is_same_v<T, double>) {
      sign = extrema < 0 ? -1.0 : 1.0;
      const double axtrema = extrema * sign;
      for (uint32_t x = 0; x <= hdist_th; ++x) {
        fdc_v[x] = sign * compute_fdc_v(axtrema, x);
        sdc_v[x] = compute_sdc_v(axtrema, x);
      }
      fdc_u = sign * compute_fdc_u(axtrema);
      sdc_u = compute_sdc_u(axtrema);
    } else if constexpr (std::is_same_v<T, cmlane_t>) {
      alignas(64) double axtrema[WIDTH];
      for (uint32_t i = 0; i < WIDTH; ++i) {
        sign[i] = extrema[i] < 0 ? -1.0 : 1.0;
        axtrema[i] = extrema[i] * sign[i];
      }
      for (uint32_t x = 0; x <= hdist_th; ++x) {
        for (uint32_t i = 0; i < WIDTH; ++i) {
          fdc_v[x][i] = sign[i] * compute_fdc_v(axtrema[i], x);
          sdc_v[x][i] = compute_sdc_v(axtrema[i], x);
        }
      }
      for (uint32_t i = 0; i < WIDTH; ++i) {
        fdc_u[i] = sign[i] * compute_fdc_u(axtrema[i]);
        sdc_u[i] = compute_sdc_u(axtrema[i]);
      }
    } else {
      static_assert(std::is_same_v<T, double> || std::is_same_v<T, cmlane_t>, "LLH supports only double or cmlane_t");
    }
  }

  // Passing a reference (&) fails due to SIMD for these
  const T get_sdc(uint32_t x) const { return sdc_v[x]; }

  const T get_fdc(uint32_t x) const { return fdc_v[x]; }

  const T get_sdc() const { return sdc_u; }

  const T get_fdc() const { return fdc_u; }

  T get_extrema() const { return extrema; }

  std::pair<uint8_t, uint8_t> get_sign_bv() const
  {
    uint8_t positive_bv = 0, negative_bv = 0;
    if constexpr (std::is_same_v<T, double>) {
      if (sign > 0)
        positive_bv = 1;
      else
        negative_bv = 1;
    } else {
      for (size_t ti = 0; ti < WIDTH; ++ti) {
        if (sign[ti] > 0)
          positive_bv |= static_cast<uint8_t>(1u << ti);
        else
          negative_bv |= static_cast<uint8_t>(1u << ti);
      }
    }
    return {positive_bv, negative_bv};
  }

  double prob_elude(uint32_t x) const
  {
    return static_cast<double>(binom_coef_hnk[x]) / static_cast<double>(binom_coef_k[x]);
  }

  double prob_collide(uint32_t x) const
  {
    return 1.0 - (static_cast<double>(binom_coef_hnk[x]) / static_cast<double>(binom_coef_k[x]));
  }

  double prob_mutate(double D, uint32_t x) const { return std::pow(1.0 - D, k - x) * std::pow(D, x) * binom_coef_k[x]; }

  double prob_miss(double D) const
  {
    double p = 0;
    for (uint32_t x = 0; x <= hdist_th; ++x) {
      p += prob_elude(x) * prob_mutate(D, x);
    }
    for (uint32_t x = hdist_th + 1; x <= k; ++x) {
      p += prob_mutate(D, x);
    }
    return (rho * p) + 1.0 - rho;
  }

  double prob_hit(double D, uint32_t x) const { return rho * prob_collide(x) * prob_mutate(D, x); }

  double nll(const double& D, const uint64_t* vv, uint64_t uu) const
  {
    if (D == 0.0) {
      for (uint32_t x = 1; x <= hdist_th; ++x)
        if (vv[x] > 0) return pinf();
    }

    double lsum = 0.0;
    double lv_m = 0.0;
    double powdc = std::pow(1.0 - D, k);
    const double logdn = k * std::log(1.0 - D);
    const double logdp = D > 0.0 ? std::log(D) - std::log(1.0 - D) : 0.0;
    const double ratioD = D / (1.0 - D);

    for (uint32_t x = 0; x <= k; ++x) {
      if (x <= hdist_th) {
        lsum -= (logdn + (x * logdp)) * vv[x];
        lv_m += dcoef_hnk[x] * powdc;
      } else {
        lv_m += powdc * dcoef_k[x];
      }
      powdc *= ratioD;
    }

    return lsum - (std::log((rho * lv_m) + 1.0 - rho) * uu);
  }

  // Observed Fisher information I(D) = -d^2/dD^2 log L(D).
  double compute_fisher_info(const uint64_t* v_r, uint64_t u_r, double D) const
  {
    double ll_dd = 0.0;
    for (uint32_t x = 0; x <= hdist_th; ++x) {
      ll_dd += static_cast<double>(v_r[x]) * compute_sdc_v(D, x);
    }
    ll_dd += static_cast<double>(u_r) * compute_sdc_u(D);
    return -ll_dd;
  }

  // NaN when there are no k-mer hits (distance undefined).
  double mle(const uint64_t* v_r, uint64_t u_r, double* nll_min = nullptr) const
  {
    uint64_t t = 0;
    for (uint32_t x = 0; x <= hdist_th; ++x)
      t += v_r[x];
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
  double compute_fdc_v(double D, uint32_t x) const
  {
    const double S = 1.0 - D;
    return (x - (k * D)) / (D * S);
  }

  double compute_sdc_v(double D, uint32_t x) const
  {
    const double S = 1.0 - D;
    const double numerator = (x * ((2 * D) - 1.0)) - (k * D * D);
    const double denominator = D * D * S * S;
    return numerator / denominator;
  }

  double compute_fdc_u(double D) const
  {
    const double S = 1.0 - D;
    double gd = 0;
    double fd = 0;
    for (uint32_t x = 0; x <= k; ++x) {
      const double pd = (x - (k * D)) / (D * S);
      const double pe = std::pow(S, k - x) * std::pow(D, x);
      double wt = pe;
      if (x <= hdist_th) {
        wt *= binom_coef_hnk[x];
      } else {
        wt *= binom_coef_k[x];
      }
      gd += wt * pd;
      fd += wt;
    }
    return rho * gd / (1.0 - rho + (rho * fd));
  }

  double compute_sdc_u(double D) const
  {
    const double S = 1.0 - D;
    double gd = 0;
    double fd = 0;
    double fpd = 0;
    double gpd = 0;
    const double D_sq = D * D;
    const double S_sq = S * S;
    const double denom = D_sq * S_sq;
    for (uint32_t x = 0; x <= k; ++x) {
      const double pd = (x - (k * D)) / (D * S);
      const double pe = std::pow(S, k - x) * std::pow(D, x);
      const double vy = ((x * x) + ((k - 1) * k * D_sq) - (x * (1 + ((k - 1) * 2 * D)))) / denom;
      double wt = pe;
      if (x <= hdist_th) {
        wt *= binom_coef_hnk[x];
      } else {
        wt *= binom_coef_k[x];
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
  // fdc/sdc: first/second d/dD log-likelihood contributions; *_v per Hamming-distance hit, *_u for misses.
  T fdc_u;
  T sdc_u;
  vec<T> fdc_v;
  vec<T> sdc_v;
};

// -2 log[L(D_ref) / L(d_mle)] from negative log-likelihoods.
inline double likelihood_ratio_statistic(double nll_ref, double nll_mle) noexcept
{
  if (!std::isfinite(nll_ref) || !std::isfinite(nll_mle)) return nanx();
  return 2.0 * std::max(0.0, nll_ref - nll_mle);
}

// Weakest detectable match: one hit at hdist_th, rest misses. Largest estimable d.
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

// LR of d vs the sketch's max estimable distance on the ceiling match counts.
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

#endif
