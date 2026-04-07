#ifndef _LLH_HPP
#define _LLH_HPP

#include <cmath>
#include <cstdint>
#include <sys/types.h>
#include <vector>
#include "types.hpp"

class PLLH
{
public:
  PLLH(uint32_t h, uint32_t k, double rho, uint32_t hdist_th)
    : h(h)
    , k(k)
    , rho(rho)
    , hdist_th(hdist_th)
    , binom_coef_k(k + 1)
    , binom_coef_hnk(hdist_th + 1)
  {
    const uint32_t nh = k - h;

    binom_coef_k[0] = 1;
    for (uint32_t d = 0; d < k; ++d) {
      binom_coef_k[d + 1] = (binom_coef_k[d] * (k - d)) / (d + 1);
    }
    binom_coef_hnk[0] = 0;
    uint64_t vc = 1;
    for (uint32_t d = 1; d <= hdist_th; ++d) {
      vc = (vc * (nh - d + 1)) / d;
      binom_coef_hnk[d] = binom_coef_k[d] - vc;
    }
  }

  const double rho;
  const uint32_t h;
  const uint32_t k;
  const uint32_t hdist_th;
  std::vector<uint64_t> binom_coef_k;
  std::vector<uint64_t> binom_coef_hnk;
};

template<typename T>
class LLH : public PLLH
{
  static constexpr size_t WIDTH = std::is_same_v<T, double> ? 1 : RWIDTH;
  static_assert(std::is_same_v<T, double> || std::is_same_v<T, cm512_t>, "LLH supports only double or cm512_t");

public:
  const T extrema;

  LLH(uint32_t k, uint32_t h, double rho, uint32_t hdist_th, T extrema)
    : PLLH(h, k, rho, hdist_th)
    , extrema(extrema)
    , fdc_v(hdist_th + 1)
    , sdc_v(hdist_th + 1)
  {
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

  void set_counts(const uint64_t* v_r, uint64_t u_r)
  {
    v = v_r;
    u = u_r;
  }

  // Passing a reference (&) fails due to SIMD for these...
  const T get_sdc(uint32_t d) const { return sdc_v[d]; }

  const T get_fdc(uint32_t d) const { return fdc_v[d]; }

  const T get_sdc() const { return sdc_u; }

  const T get_fdc() const { return fdc_u; }

  T get_sign() const { return sign; }

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
    return 1.0 - static_cast<double>(binom_coef_hnk[d]) / static_cast<double>(binom_coef_k[d]);
  }

  double prob_collide(uint32_t d) const
  {
    return static_cast<double>(binom_coef_hnk[d]) / static_cast<double>(binom_coef_k[d]);
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
    return rho * p + 1.0 - rho;
  }

  double prob_hit(double D, uint32_t d) const { return rho * prob_collide(d) * prob_mutate(D, d); }

  double operator()(const double& D) const
  {
    double lsum = 0.0;
    double lv_m = 0.0;
    double powdc = std::pow(1.0 - D, k);
    const double logdn = k * std::log(1.0 - D);
    const double logdp = std::log(D) - std::log(1.0 - D);
    const double ratioD = D / (1.0 - D);

    for (uint32_t d = 0; d <= k; ++d) {
      if (d <= hdist_th) {
        lsum -= (logdn + d * logdp) * v[d];
        lv_m += binom_coef_hnk[d] * powdc;
      } else {
        lv_m += powdc * binom_coef_k[d];
      }
      powdc *= ratioD;
    }

    return lsum - std::log(rho * lv_m + 1.0 - rho) * u;
  }

private:
  double compute_fdc_v(const double D, const uint32_t d) const
  {
    const double S = 1.0 - D;
    return (d - k * D) / (D * S);
  }

  double compute_sdc_v(const double D, const uint32_t d) const
  {
    const double S = 1.0 - D;
    const double numerator = d * (2 * D - 1.0) - (k * D * D);
    const double denominator = D * D * S * S;
    return numerator / denominator;
  }

  double compute_fdc_u(const double D) const
  {
    const double S = 1.0 - D;
    double gd = 0;
    double fd = 0;
    for (uint32_t d = 0; d <= k; ++d) {
      const double pd = (d - k * D) / (D * S);
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
    return rho * gd / (1.0 - rho + rho * fd);
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
      const double pd = (d - k * D) / (D * S);
      const double pe = std::pow(S, k - d) * std::pow(D, d);
      const double vy = (d * d + (k - 1) * k * D_sq - d * (1 + (k - 1) * 2 * D)) / denom;
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
    fd = 1.0 - rho + rho * fd;
    return (fd * gpd - gd * fpd) / (fd * fd);
  }

  const uint64_t* v = nullptr;
  uint64_t u = 0;
  T sign;
  T fdc_u;
  T sdc_u;
  std::vector<T> fdc_v;
  std::vector<T> sdc_v;
};

#endif
