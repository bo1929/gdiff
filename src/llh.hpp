#ifndef _LLH_H
#define _LLH_H

#include <cmath>
#include <sys/types.h>

class LLH
{
private:
  double extrema;
  double ixtrema;
  uint64_t* v;
  uint64_t u;
  std::vector<uint64_t> binom_coef_k;
  std::vector<uint64_t> binom_coef_hnk;
  std::vector<double> fdc_v;
  std::vector<double> sdc_v;
  double fdc_u;
  double sdc_u;

public:
  const uint32_t k;
  const uint32_t h;
  const double rho;
  const uint32_t hdist_th;
  const bool opposite;

  LLH(uint32_t h, uint32_t k, uint32_t hdist_th, double extrema, double rho)
    : h(h)
    , k(k)
    , hdist_th(hdist_th)
    , extrema(extrema)
    , opposite(extrema < 0)
    , rho(rho)
  {
    uint64_t vc = 1;
    uint32_t nh = k - h;

    binom_coef_k.resize(k + 1);
    binom_coef_hnk.resize(hdist_th + 1);
    binom_coef_k[0] = 1;
    binom_coef_hnk[0] = 0;
    for (int32_t d = 0; d < k; ++d) {
      binom_coef_k[d + 1] = (binom_coef_k[d] * (k - d)) / (d + 1);
    }
    for (uint32_t d = 1; d <= hdist_th; ++d) {
      vc = (vc * (nh - d + 1)) / d;
      binom_coef_hnk[d] = binom_coef_k[d] - vc;
    }

    if (opposite) {
      extrema = -extrema;
    }
    ixtrema = 1.0 - extrema;
    fdc_v.resize(hdist_th + 1);
    sdc_v.resize(hdist_th + 1);
    for (uint32_t d = 0; d <= hdist_th; ++d) {
      fdc_v[d] = opposite ? -compute_fdc_v(d) : compute_fdc_v(d);
      sdc_v[d] = opposite ? -compute_sdc_v(d) : compute_sdc_v(d);
      ;
    }
    fdc_u = opposite ? -compute_fdc_u() : compute_fdc_u();
    sdc_u = opposite ? -compute_sdc_u() : compute_sdc_u();
  }

  void set_counts(uint64_t* v_r, uint64_t u_r)
  {
    v = v_r;
    u = u_r;
  }

  double compute_fdc_v(uint32_t d) { return (d - k * extrema) / (extrema * ixtrema); }

  double compute_sdc_v(uint32_t d)
  {
    return (static_cast<double>(d) * (2 * extrema - 1) - k * extrema * extrema) / (extrema * extrema * ixtrema * ixtrema);
  }

  double compute_fdc_u()
  {
    double vnp = 0;
    double vdp = 0;
    for (uint32_t d = 0; d <= k; ++d) {
      const double pd = (d - (k * extrema)) / (extrema * (1.0 - extrema));
      const double pe = pow((1.0 - extrema), k - d) * pow(extrema, d);
      const double pc = (d > hdist_th) + (d <= hdist_th) * (binom_coef_hnk[d] / static_cast<double>(binom_coef_k[d]));
      vnp += pc * binom_coef_k[d] * pe * pd;
      vdp += pe * binom_coef_k[d] * pc;
    }
    return rho * (vnp / (1.0 - rho + rho * vdp));
  }

  double compute_sdc_u()
  {
    double gd = 0;
    double ffd = 0;
    double fpd = 0;
    double gpd = 0;
    for (uint32_t d = 0; d <= k; ++d) {
      const double pd = (d - (k * extrema)) / (extrema * (1.0 - extrema));
      const double pe = pow((1.0 - extrema), k - d) * pow(extrema, d);
      const double pc = (d > hdist_th) + (d <= hdist_th) * (binom_coef_hnk[d] / static_cast<double>(binom_coef_k[d]));
      const double vy = ((d * d + (k - 1) * k * extrema * extrema) - d * (1 + (k - 1) * 2 * extrema)) /
                        ((extrema * ixtrema) * (extrema * ixtrema));

      gd += pc * binom_coef_k[d] * pe * pd;
      gpd += pc * binom_coef_k[d] * pe * vy;
      ffd += (pe * binom_coef_k[d] * pc);
      fpd += (pe * binom_coef_k[d] * pd * pc);
    }
    gd = gd * rho;
    fpd = fpd * rho;
    gpd = gpd * rho;
    const double fd = (1.0 - rho + rho * ffd);
    return (((fd * gpd) - (gd * fpd)) / (fd * fd));
  }

  double get_sdc(uint32_t d) { return sdc_v[d]; }

  double get_fdc(uint32_t d) { return fdc_v[d]; }

  double get_sdc() { return sdc_u; }

  double get_fdc() { return fdc_u; }

  double prob_elude(const uint32_t d)
  {
    return 1.0 - static_cast<double>(binom_coef_hnk[d]) / static_cast<double>(binom_coef_k[d]);
  }

  double prob_collide(const uint32_t d)
  {
    return static_cast<double>(binom_coef_hnk[d]) / static_cast<double>(binom_coef_k[d]);
  }

  double prob_mutate(const double D, const uint32_t d) { return pow((1.0 - D), (k - d)) * pow(D, d) * binom_coef_k[d]; }

  double prob_miss(const double D)
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

  double prob_hit(const double D, const uint32_t d) { return rho * prob_collide(d) * prob_mutate(D, d); }

  double operator()(double const& D)
  {
    double sum = 0.0, lv_m = 0.0;
    double powdc = pow((1.0 - D), k);
    double logdn = log(1.0 - D);
    double logdp = log(D) - logdn;
    logdn *= k;
    double ratioD = D / (1.0 - D);
    for (uint32_t d = 0; d <= k; ++d) {
      if (d <= hdist_th) {
        sum -= (logdn + d * logdp) * (*(v + d));
        lv_m += binom_coef_hnk[d] * powdc;
      } else {
        lv_m += powdc * binom_coef_k[d];
      }
      powdc *= ratioD;
    }
    return sum - log(rho * lv_m + 1.0 - rho) * u;
  }
};

#endif
