#ifndef _GAMMA_HPP
#define _GAMMA_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <boost/math/distributions/gamma.hpp>
#include "types.hpp"

typedef boost::math::policies::policy<boost::math::policies::max_series_iterations<1000000>> hpolicy;

// Latent-Gamma model with additive Gaussian noise:
//   d_i = D_i + N(0, sigma2_i),  D_i ~ Gamma(shape, scale)
//
// Marginals are integrated via 15-point Gauss-Hermite quadrature.
// Fitting uses 2D Nelder-Mead over (log shape, log scale).
struct GammaModel
{
  // Public types

  struct params_t
  {
    double shape; // > 0
    double scale; // > 0
  };

  [[nodiscard]] static bool validate_params(const params_t& gp) noexcept
  {
    return gp.shape > 0.0 && gp.scale > 0.0 && std::isfinite(gp.shape) && std::isfinite(gp.scale);
  }

  enum class Method
  {
    QM,  // noise-aware quantile matching at cfg.quantile_probs (default);
         // fewer degrees of information but more robust to outliers
    MEVL // marginal evidence: maximize the full marginal log-likelihood
  };

  // All runtime tuning knobs. Defaults reproduce the original behaviour.
  struct Config
  {
    // Fitting
    Method method = Method::QM;
    std::array<double, 3> quantile_probs = {0.25, 0.50, 0.75};

    // Nelder-Mead
    int max_niter = 200;   // iteration cap
    double tol = 1e-6;     // convergence threshold (objective spread and simplex diameter)
    double nm_step = 0.5;  // initial simplex offset in log-space
    double nm_alpha = 1.0; // reflection coefficient
    double nm_gamma = 2.0; // expansion coefficient
    double nm_rho = 0.5;   // contraction coefficient
    double nm_sigma = 0.5; // shrink coefficient
  };

  // Constants

  static constexpr double eps = 1e-10;
  static constexpr double sqrt2 = 1.4142135623730951;
  static constexpr double sqrt_pi = 1.7724538509055159;
  static constexpr size_t min_nsamples = 8;

  // 15-point Gauss-Hermite nodes and weights for int e^{-x^2} f(x) dx.
  // More points extend tail accuracy, which matters when marginal_cdf is
  // evaluated far from the null median.
  static constexpr std::array<double, 15> gh_nodes = {-4.49999070730939,
                                                      -3.66995037340445,
                                                      -2.96716692790560,
                                                      -2.32573248617386,
                                                      -1.71999257518649,
                                                      -1.13611558521092,
                                                      -0.565069583255576,
                                                      0.0,
                                                      0.565069583255576,
                                                      1.13611558521092,
                                                      1.71999257518649,
                                                      2.32573248617386,
                                                      2.96716692790560,
                                                      3.66995037340445,
                                                      4.49999070730939};
  static constexpr std::array<double, 15> gh_weights = {1.52247580425352e-9,
                                                        1.05911554771107e-6,
                                                        1.00004441232499e-4,
                                                        2.77806884291280e-3,
                                                        3.07800338725462e-2,
                                                        1.58488915795936e-1,
                                                        4.12028687498899e-1,
                                                        5.64100308726418e-1,
                                                        4.12028687498899e-1,
                                                        1.58488915795936e-1,
                                                        3.07800338725462e-2,
                                                        2.77806884291280e-3,
                                                        1.00004441232499e-4,
                                                        1.05911554771107e-6,
                                                        1.52247580425352e-9};
  static_assert(gh_nodes.size() == gh_weights.size(), "GH node and weight arrays must have the same length.");

  // Gamma distribution primitives

  // Returns 0 for x <= 0 or invalid parameters.
  [[nodiscard]] static double gamma_pdf(double x, double shape, double scale)
  {
    if (!(x > 0.0) || !(shape > 0.0) || !(scale > 0.0)) return 0.0;
    return boost::math::pdf(boost::math::gamma_distribution<double, hpolicy>(shape, scale), x);
  }

  // Returns 0 for x <= 0, NaN for invalid parameters.
  [[nodiscard]] static double gamma_cdf(double x, double shape, double scale)
  {
    if (!(shape > 0.0) || !(scale > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    if (!(x > 0.0)) return 0.0;
    return boost::math::cdf(boost::math::gamma_distribution<double, hpolicy>(shape, scale), x);
  }

  // Marginal distribution (Gamma convolved with Gaussian noise)

  // Marginal pdf of d = D + N(0, sigma^2), D ~ Gamma(shape, scale).
  // Falls back to gamma_pdf when sigma is not positive-finite.
  [[nodiscard]] static double marginal_pdf(double d, double sigma, double shape, double scale)
  {
    if (!(sigma > 0.0) || !std::isfinite(sigma)) return gamma_pdf(d, shape, scale);
    const double s = sigma * sqrt2;
    return gh_integrate([&](double xi) { return gamma_pdf(d - s * xi, shape, scale); });
  }

  // Marginal cdf P(D + N(0, sigma^2) <= t).
  // Falls back to gamma_cdf when sigma is not positive-finite.
  [[nodiscard]] static double marginal_cdf(double t, double sigma, double shape, double scale)
  {
    if (!(sigma > 0.0) || !std::isfinite(sigma)) return gamma_cdf(t, shape, scale);
    const double s = sigma * sqrt2;
    return std::clamp(gh_integrate([&](double xi) { return gamma_cdf(t - s * xi, shape, scale); }), 0.0, 1.0);
  }

  // Objective functions

  // Sum_i -log p(d_i | shape, scale, sigma_i). Returns +inf for invalid params.
  [[nodiscard]] static double neg_log_likelihood(const vec<double>& d_v, const vec<double>& s2_v, double shape, double scale)
  {
    return nll_from_sigmas(d_v, make_sigmas(d_v.size(), s2_v), shape, scale);
  }

  // Fitting

  // Fit (shape, scale) from noisy observations d_v.
  // We have per-observation noise variances s2_v.
  // Returns {1, 1} for empty input.
  [[nodiscard]] static params_t fit(const vec<double>& d_v, const vec<double>& s2_v)
  {
    Config cfg{};
    return fit(d_v, s2_v, cfg);
  }

  [[nodiscard]] static params_t fit(const vec<double>& d_v, const vec<double>& s2_v, const Config& cfg)
  {
    if (d_v.empty()) return {1.0, 1.0};
    const params_t p0 = moment_init(d_v, s2_v);
    const vec<double> sigmas = make_sigmas(d_v.size(), s2_v);

    if (cfg.method == Method::MEVL) {
      return nelder_mead_bivariate(
        p0, [&](double shape, double scale) { return nll_from_sigmas(d_v, sigmas, shape, scale); }, cfg);
    }

    // QM: minimise sum_k (mean_i F(q_k; sigma_i, shape, scale) - p_k)^2
    const auto emp_q = compute_quantiles(d_v, cfg.quantile_probs);
    return nelder_mead_bivariate(
      p0,
      [&](double shape, double scale) { return qm_loss_from_sigmas(emp_q, sigmas, cfg.quantile_probs, shape, scale); },
      cfg);
  }

  // Fit Gamma directly to posterior draws (no noise model).
  // Uses plain quantile matching -- faster than MEVL (no quadrature).
  [[nodiscard]] static params_t fit_from_samples(const vec<double>& d_v)
  {
    Config cfg{};
    return fit_from_samples(d_v, cfg);
  }

  [[nodiscard]] static params_t fit_from_samples(const vec<double>& d_v, const Config& cfg)
  {
    if (d_v.empty()) return {1.0, 1.0};
    const params_t p0 = moment_init(d_v, {});
    const auto emp_q = compute_quantiles(d_v, cfg.quantile_probs);
    return nelder_mead_bivariate(
      p0,
      [&](double shape, double scale) {
        if (!(shape > 0.0) || !(scale > 0.0)) return std::numeric_limits<double>::infinity();
        double L = 0.0;
        for (size_t k = 0; k < emp_q.size(); ++k) {
          const double diff = gamma_cdf(emp_q[k], shape, scale) - cfg.quantile_probs[k];
          L += diff * diff;
        }
        return L;
      },
      cfg);
  }

  // Score an observed distance against null samples (plain Gamma, no noise convolution).
  // Returns {gamma_cdf(d_obs), latent median} or NaN components on failure.
  [[nodiscard]] static xy_t
  score_from_samples(const double d_obs, const vec<double>& d_v, const double d_lo, const double d_hi)
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (!std::isfinite(d_obs)) return {nan, nan};
    if (d_v.size() < min_nsamples) return {nan, nan};

    const params_t gp = fit_from_samples(d_v);
    if (!validate_params(gp)) return {nan, nan};

    const double prob = gamma_cdf(d_obs, gp.shape, gp.scale);
    if (!std::isfinite(prob)) return {nan, nan};

    auto F = [&](double x) { return gamma_cdf(x, gp.shape, gp.scale); };
    if (!(d_hi > d_lo) || F(d_hi) < 0.5) return {prob, nan};

    // Median of the latent Gamma distribution by binary search between distance lower and upper bounds.
    double lo = d_lo;
    double hi = d_hi;
    for (int it = 0; it < 40; ++it) {
      const double mid = 0.5 * (lo + hi);
      if (F(mid) < 0.5)
        lo = mid;
      else
        hi = mid;
    }
    const double median = 0.5 * (lo + hi);
    if (!(median > eps) || !std::isfinite(median)) return {prob, nan};
    return {prob, median};
  }

private:
  // Internal types

  // 2D point in log-space, used by Nelder-Mead.
  struct Point2
  {
    double x, y;
  };

  // Private helpers

  // Gauss-Hermite quadrature: (1/sqrt(pi)) * int e^{-x^2} f(x) dx.
  template<typename F>
  [[nodiscard]] static double gh_integrate(F&& f)
  {
    double acc = 0.0;
    for (size_t j = 0; j < gh_nodes.size(); ++j)
      acc += gh_weights[j] * f(gh_nodes[j]);
    return acc / sqrt_pi;
  }

  // Resolve sigma from a variance vector. Returns 0 if the entry is absent or invalid.
  [[nodiscard]] static double resolve_sigma(const vec<double>& s2_v, size_t i)
  {
    return (i < s2_v.size() && s2_v[i] > 0.0 && std::isfinite(s2_v[i])) ? std::sqrt(s2_v[i]) : 0.0;
  }

  // Pre-compute a dense sigma vector so hot-path objectives avoid repeated
  // sqrt() and bounds-checks. Called once per fit(), not per NM step.
  [[nodiscard]] static vec<double> make_sigmas(size_t n, const vec<double>& s2_v)
  {
    vec<double> sigmas(n);
    for (size_t i = 0; i < n; ++i)
      sigmas[i] = resolve_sigma(s2_v, i);
    return sigmas;
  }

  // NLL using a pre-computed sigma vector (hot path inside NM).
  [[nodiscard]] static double nll_from_sigmas(const vec<double>& d_v, const vec<double>& sigmas, double shape, double scale)
  {
    if (!(shape > 0.0) || !(scale > 0.0)) return std::numeric_limits<double>::infinity();
    double nll = 0.0;
    for (size_t i = 0; i < d_v.size(); ++i) {
      const double p = marginal_pdf(d_v[i], sigmas[i], shape, scale);
      nll -= std::log(std::max(p, 1e-300));
    }
    return nll;
  }

  // QM loss using a pre-computed sigma vector (hot path inside NM).
  // Loss = sum_k (mean_i F(q_k | sigma_i, shape, scale) - p_k)^2
  [[nodiscard]] static double qm_loss_from_sigmas(const std::array<double, 3>& emp_q,
                                                  const vec<double>& sigmas,
                                                  const std::array<double, 3>& probs,
                                                  double shape,
                                                  double scale)
  {
    if (!(shape > 0.0) || !(scale > 0.0)) return std::numeric_limits<double>::infinity();
    const double inv_n = 1.0 / static_cast<double>(sigmas.size());
    double L = 0.0;
    for (size_t k = 0; k < probs.size(); ++k) {
      double Fk = 0.0;
      for (size_t i = 0; i < sigmas.size(); ++i)
        Fk += marginal_cdf(emp_q[k], sigmas[i], shape, scale);
      const double diff = Fk * inv_n - probs[k];
      L += diff * diff;
    }
    return L;
  }

  // Moment-based initializer. Subtracts mean noise variance from sample variance
  // to estimate latent Gamma moments before optimization begins.
  [[nodiscard]] static params_t moment_init(const vec<double>& d_v, const vec<double>& s2_v)
  {
    const double N = static_cast<double>(d_v.size());
    double m = 0.0, v = 0.0, mean_s2 = 0.0;
    for (double x : d_v)
      m += x;
    m /= N;
    for (double x : d_v)
      v += (x - m) * (x - m);
    v /= N;
    for (double s2 : s2_v)
      if (s2 > 0.0 && std::isfinite(s2)) mean_s2 += s2;
    if (!s2_v.empty()) mean_s2 /= static_cast<double>(s2_v.size());

    const double v_lat = std::max(v - mean_s2, 1e-6);
    const double shape0 = (m > eps) ? (m * m) / v_lat : 1.0;
    const double scale0 = (m > eps) ? v_lat / m : std::max(m, eps);
    return {std::max(shape0, 1e-2), std::max(scale0, 1e-6)};
  }

  // Compute interpolated empirical quantiles of d_v at the given probabilities.
  [[nodiscard]] static std::array<double, 3> compute_quantiles(const vec<double>& d_v, const std::array<double, 3>& probs)
  {
    vec<double> sorted = d_v;
    std::sort(sorted.begin(), sorted.end());
    const double last = static_cast<double>(sorted.size() - 1);
    std::array<double, 3> q{};
    for (size_t k = 0; k < probs.size(); ++k) {
      const double ix = probs[k] * last;
      const size_t lo = static_cast<size_t>(ix);
      const size_t hi = std::min(lo + 1, sorted.size() - 1);
      const double frac = ix - static_cast<double>(lo);
      q[k] = sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
    }
    return q;
  }

  // Sort a 3-element simplex (S, F) so that F[0] <= F[1] <= F[2].
  // Optimal 3-element sort network (3 conditional swaps, no branches on average).
  static void sort3(std::array<Point2, 3>& S, std::array<double, 3>& F)
  {
    if (F[0] > F[1]) {
      std::swap(F[0], F[1]);
      std::swap(S[0], S[1]);
    }
    if (F[1] > F[2]) {
      std::swap(F[1], F[2]);
      std::swap(S[1], S[2]);
    }
    if (F[0] > F[1]) {
      std::swap(F[0], F[1]);
      std::swap(S[0], S[1]);
    }
  }

  // 2D Nelder-Mead in (log shape, log scale). obj receives linear-space params.
  // Stops when both objective spread and simplex diameter fall below cfg.tol;
  // the diameter check catches flat-likelihood drift that objective-only misses.
  template<typename Obj>
  [[nodiscard]] static params_t nelder_mead_bivariate(params_t p0, Obj&& obj, const Config& cfg)
  {
    auto eval = [&](const Point2& p) { return obj(std::exp(p.x), std::exp(p.y)); };

    std::array<Point2, 3> S = {{{std::log(p0.shape), std::log(p0.scale)},
                                {std::log(p0.shape) + cfg.nm_step, std::log(p0.scale)},
                                {std::log(p0.shape), std::log(p0.scale) + cfg.nm_step}}};
    std::array<double, 3> F = {eval(S[0]), eval(S[1]), eval(S[2])};
    sort3(S, F);

    for (int iter = 0; iter < cfg.max_niter; ++iter) {
      const bool converged =
        (F[2] - F[0]) < cfg.tol && std::abs(S[2].x - S[0].x) < cfg.tol && std::abs(S[2].y - S[0].y) < cfg.tol;
      if (converged && iter > 5) break;

      const Point2 c = {0.5 * (S[0].x + S[1].x), 0.5 * (S[0].y + S[1].y)};
      const Point2 r = {c.x + cfg.nm_alpha * (c.x - S[2].x), c.y + cfg.nm_alpha * (c.y - S[2].y)};
      const double fr = eval(r);

      if (fr < F[0]) {
        // Expansion
        const Point2 e = {c.x + cfg.nm_gamma * (r.x - c.x), c.y + cfg.nm_gamma * (r.y - c.y)};
        const double fe = eval(e);
        S[2] = (fe < fr) ? e : r;
        F[2] = std::min(fe, fr);
      } else if (fr < F[1]) {
        // Accept reflection
        S[2] = r;
        F[2] = fr;
      } else {
        // Contraction toward the better of r and S[2]
        const bool use_r = fr < F[2];
        const Point2 w = use_r ? r : S[2];
        const double fw = use_r ? fr : F[2];
        const Point2 cc = {c.x + cfg.nm_rho * (w.x - c.x), c.y + cfg.nm_rho * (w.y - c.y)};
        const double fc = eval(cc);
        if (fc <= fw) {
          S[2] = cc;
          F[2] = fc;
        } else {
          // Shrink: pull all non-best vertices toward S[0]
          for (int j = 1; j < 3; ++j) {
            S[j] = {S[0].x + cfg.nm_sigma * (S[j].x - S[0].x), S[0].y + cfg.nm_sigma * (S[j].y - S[0].y)};
            F[j] = eval(S[j]);
          }
        }
      }
      sort3(S, F);
    }

    // S[0] is always the best vertex after sort3.
    return {std::exp(S[0].x), std::exp(S[0].y)};
  }
};

#endif // _GAMMA_HPP
