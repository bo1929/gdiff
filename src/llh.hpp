#ifndef _LLH_H
#define _LLH_H

#include <cmath>
#include <vector>
#include <cstdint>
#include <sys/types.h>

class LLH
{
public:
  const uint32_t h;
  const uint32_t k;
  const uint32_t hdist_th;
  const double extrema;
  const double ixtrema;
  const bool opposite;
  const double rho;

  LLH(uint32_t h, uint32_t k, uint32_t hdist_th, double extrema, double rho)
    : h(h)
    , k(k)
    , hdist_th(hdist_th)
    , extrema(extrema < 0 ? -extrema : extrema)
    , ixtrema(1.0 - (extrema < 0 ? -extrema : extrema))
    , opposite(extrema < 0)
    , rho(rho)
    , binom_coef_k(k + 1)
    , binom_coef_hnk(hdist_th + 1)
    , fdc_v(hdist_th + 1)
    , sdc_v(hdist_th + 1)
  {
    const uint32_t nh = k - h;

    binom_coef_k[0] = 1;
    binom_coef_hnk[0] = 0;
    for (uint32_t d = 0; d < k; ++d) {
      binom_coef_k[d + 1] = (binom_coef_k[d] * (k - d)) / (d + 1);
    }
    uint64_t vc = 1;
    for (uint32_t d = 1; d <= hdist_th; ++d) {
      vc = (vc * (nh - d + 1)) / d;
      binom_coef_hnk[d] = binom_coef_k[d] - vc;
    }

    for (uint32_t d = 0; d <= hdist_th; ++d) {
      const double fdc = compute_fdc_v(d);
      const double sdc = compute_sdc_v(d);
      fdc_v[d] = opposite ? -fdc : fdc;
      sdc_v[d] = sdc;
    }
    const double fdc = compute_fdc_u();
    const double sdc = compute_sdc_u();
    fdc_u = opposite ? -fdc : fdc;
    sdc_u = sdc;
  }

  void set_counts(uint64_t* v_r, uint64_t u_r)
  {
    v = v_r;
    u = u_r;
  }

  double get_sdc(uint32_t d) const { return sdc_v[d]; }

  double get_fdc(uint32_t d) const { return fdc_v[d]; }

  double get_sdc() const { return sdc_u; }

  double get_fdc() const { return fdc_u; }

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
  double compute_fdc_v(uint32_t d) const { return (d - k * extrema) / (extrema * ixtrema); }

  double compute_sdc_v(uint32_t d) const
  {
    const double numerator = d * (2 * extrema - 1.0) - (k * extrema * extrema);
    const double denominator = extrema * extrema * ixtrema * ixtrema;
    return numerator / denominator;
  }

  double compute_fdc_u() const
  {
    double vnp = 0;
    double vdp = 0;
    for (uint32_t d = 0; d <= k; ++d) {
      const double pd = (d - k * extrema) / (extrema * ixtrema);
      const double pe = std::pow(ixtrema, k - d) * std::pow(extrema, d);
      const double pc = (d > hdist_th) ? 1.0 : (binom_coef_hnk[d] / static_cast<double>(binom_coef_k[d]));
      const double wt = binom_coef_k[d] * pe * pc;
      vnp += wt * pd;
      vdp += wt;
    }
    return rho * vnp / (1.0 - rho + rho * vdp);
  }

  double compute_sdc_u() const
  {
    double ggd = 0;
    double ffd = 0;
    double fpd = 0;
    double gpd = 0;
    const double extrema_sq = extrema * extrema;
    const double ixtrema_sq = ixtrema * ixtrema;
    const double denom = extrema_sq * ixtrema_sq;
    for (uint32_t d = 0; d <= k; ++d) {
      const double pd = (d - k * extrema) / (extrema * ixtrema);
      const double pe = std::pow(ixtrema, k - d) * std::pow(extrema, d);
      const double vy = (d * d + (k - 1) * k * extrema_sq - d * (1 + (k - 1) * 2 * extrema)) / denom;
      const double wt = binom_coef_k[d] * pe * (d > hdist_th) +
                        binom_coef_k[d] * (d <= hdist_th) * ((pe * binom_coef_hnk[d]) / binom_coef_k[d]);
      ggd += wt * pd;
      ffd += wt;
      fpd += wt * pd;
      gpd += wt * vy;
    }
    ggd *= rho;
    fpd *= rho;
    gpd *= rho;
    const double fd = 1.0 - rho + rho * ffd;
    return (fd * gpd - ggd * fpd) / (fd * fd);
  }

  uint64_t* v = nullptr;
  uint64_t u = 0;
  double fdc_u;
  double sdc_u;
  std::vector<double> fdc_v;
  std::vector<double> sdc_v;
  std::vector<uint64_t> binom_coef_k;
  std::vector<uint64_t> binom_coef_hnk;
};

#endif
