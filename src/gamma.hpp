#ifndef _GAMMA_HPP
#define _GAMMA_HPP

#include <algorithm>
#include <boost/math/distributions/gamma.hpp>
#include "types.hpp"

struct GammaModel
{
  struct params_t
  {
    double alpha;
    double beta;
  };
  static constexpr double eps = 1e-10;
  static constexpr double tol = 1e-6;
  static constexpr int max_niter = 250;
  static constexpr uint64_t min_nsamples = 50;
  static constexpr std::array<double, 3> fit_probs = {0.25, 0.50, 0.75};

  // Fit Gamma to a sorted sample vector using quantile matching + 2D Nelder-Mead
  static inline params_t fit(const vec<double>& sorted_samples)
  {
    const size_t n = sorted_samples.size();
    if (n < 3) return {1.0, 1.0};

    // Compute empirical quantiles
    std::array<double, fit_probs.size()> emp_q;
    for (size_t i = 0; i < fit_probs.size(); ++i) {
      const double idx = fit_probs[i] * (n - 1);
      const size_t low = static_cast<size_t>(idx);
      const size_t high = std::min(low + 1, n - 1);
      const double frac = idx - low;
      emp_q[i] = (sorted_samples[low] * (1.0 - frac)) + (sorted_samples[high] * frac);
    }

    if (emp_q.back() < eps) return {1.0, eps};

    // Mean distance for initialization
    double d_mean = 0.0;
    for (double x : sorted_samples)
      d_mean += x;
    d_mean /= n;
    if (d_mean < eps) d_mean = eps;

    // Our objective over both shape (alpha) and scale (beta)
    auto sse = [&](double x, double y) -> double {
      if (x < eps || y < eps) return 1e30;
      boost::math::gamma_distribution<double> g(x, y);
      double s = 0.0;
      for (size_t i = 0; i < fit_probs.size(); ++i) {
        const double diff = quantile(g, fit_probs[i]) - emp_q[i];
        s += diff * diff;
      }
      return s;
    };

    // Nelder-Mead on (alpha, beta), initialized at (1, d_mean)
    using pt = params_t;
    std::array<pt, 3> S = {pt{1.0, d_mean}, pt{2.0, d_mean}, pt{1.0, 2.0 * d_mean}};
    std::array<double, 3> F;
    for (int i = 0; i < 3; ++i)
      F[i] = sse(S[i].alpha, S[i].beta);

    for (int iter = 0; iter < max_niter; ++iter) {
      // Sort vertices by objective
      if (F[0] > F[1]) {
        std::swap(S[0], S[1]);
        std::swap(F[0], F[1]);
      }
      if (F[1] > F[2]) {
        std::swap(S[1], S[2]);
        std::swap(F[1], F[2]);
      }
      if (F[0] > F[1]) {
        std::swap(S[0], S[1]);
        std::swap(F[0], F[1]);
      }

      if (F[2] - F[0] < tol && iter > 10) break;

      pt c = {0.5 * (S[0].alpha + S[1].alpha), 0.5 * (S[0].beta + S[1].beta)};
      pt r = {(2.0 * c.alpha) - S[2].alpha, (2.0 * c.beta) - S[2].beta};
      double fr = sse(r.alpha, r.beta);

      if (fr < F[0]) {
        pt e = {c.alpha + (2.0 * (r.alpha - c.alpha)), c.beta + (2.0 * (r.beta - c.beta))};
        double fe = sse(e.alpha, e.beta);
        if (fe < fr) {
          S[2] = e;
          F[2] = fe;
        } else {
          S[2] = r;
          F[2] = fr;
        }
      } else if (fr < F[1]) {
        S[2] = r;
        F[2] = fr;
      } else {
        pt w = (fr < F[2]) ? r : S[2];
        double fw = std::min(fr, F[2]);
        pt cc = {0.5 * (c.alpha + w.alpha), 0.5 * (c.beta + w.beta)};
        double fc = sse(cc.alpha, cc.beta);
        if (fc <= fw) {
          S[2] = cc;
          F[2] = fc;
        } else {
          for (int j = 1; j < 3; ++j) {
            S[j].alpha = 0.5 * (S[0].alpha + S[j].alpha);
            S[j].beta = 0.5 * (S[0].beta + S[j].beta);
            F[j] = sse(S[j].alpha, S[j].beta);
          }
        }
      }
    }

    // Return best vertex
    int bidx = 0;
    if (F[1] < F[bidx]) bidx = 1;
    if (F[2] < F[bidx]) bidx = 2;
    return {std::max(S[bidx].alpha, eps), std::max(S[bidx].beta, eps)};
  }
};

#endif
