#ifndef _GAMMA_HPP
#define _GAMMA_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>
#include <boost/math/distributions/gamma.hpp>

// Gamma distribution model.
// Fitting uses 2D Nelder-Mead over (log shape, log scale).
class GammaModel
{
public:
  struct params_t
  {
    double shape; // > 0
    double scale; // > 0
  };

  [[nodiscard]] static bool validate_params(const params_t& gp) noexcept
  {
    return gp.shape > 0.0 && gp.scale > 0.0 && std::isfinite(gp.shape) && std::isfinite(gp.scale);
  }

  struct Config
  {
    std::array<double, 3> quantile_probs = {0.25, 0.50, 0.75};

    // Nelder-Mead
    int max_niter = 200;    // iteration cap
    double tol = 1e-6;      // convergence threshold (objective spread and simplex diameter)
    double nm_step = 0.5;   // initial simplex offset in log-space
    double nm_alpha = 1.0;  // reflection coefficient
    double nm_expand = 2.0; // expansion coefficient
    double nm_rho = 0.5;    // contraction coefficient
    double nm_sigma = 0.5;  // shrink coefficient
  };

  // Constants

  static constexpr double eps = 1e-10;
  static constexpr size_t min_nsamples = 8;
  static constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
  static constexpr double INF = std::numeric_limits<double>::infinity();

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
    if (!(shape > 0.0) || !(scale > 0.0)) return NaN;
    if (!(x > 0.0)) return 0.0;
    return boost::math::cdf(boost::math::gamma_distribution<double, hpolicy>(shape, scale), x);
  }

  // Fit Gamma directly to draws via quantile matching.
  [[nodiscard]] static params_t fit_from_samples(const std::vector<double>& x_v)
  {
    Config cfg{};
    return fit_from_samples(x_v, cfg);
  }

  [[nodiscard]] static params_t fit_from_samples(const std::vector<double>& x_v, const Config& cfg)
  {
    if (x_v.empty()) return {1.0, 1.0};
    const params_t p0 = init_from_moments(x_v);
    const auto emp_q = compute_quantiles(x_v, cfg.quantile_probs);
    return nelder_mead_bivariate(
      p0,
      [&](double shape, double scale) {
        if (!(shape > 0.0) || !(scale > 0.0)) return INF;
        double L = 0.0;
        for (size_t k = 0; k < emp_q.size(); ++k) {
          const double diff = gamma_cdf(emp_q[k], shape, scale) - cfg.quantile_probs[k];
          L += diff * diff;
        }
        return L;
      },
      cfg);
  }

  [[nodiscard]] static std::pair<double, double>
  score_from_samples(double x, const std::vector<double>& samples, double lower, double upper)
  {
    if (!std::isfinite(x)) return {NaN, NaN};
    if (samples.size() < min_nsamples) return {NaN, NaN};

    const params_t gp = fit_from_samples(samples);
    if (!validate_params(gp)) return {NaN, NaN};

    const double prob = gamma_cdf(x, gp.shape, gp.scale);
    if (!std::isfinite(prob)) return {NaN, NaN};

    auto F = [&](double t) { return gamma_cdf(t, gp.shape, gp.scale); };
    if (!(upper > lower) || F(upper) < 0.5) return {prob, NaN};

    // Median of the latent Gamma distribution by binary search.
    double lo = lower;
    double hi = upper;
    for (int it = 0; it < 40; ++it) {
      const double mid = 0.5 * (lo + hi);
      if (F(mid) < 0.5)
        lo = mid;
      else
        hi = mid;
    }
    const double median = 0.5 * (lo + hi);
    if (!(median > eps) || !std::isfinite(median)) return {prob, NaN};
    return {prob, median};
    // Returns {gamma_cdf(x), latent median} or NaN components on failure.
  }

private:
  // 2D coordinate in log-space, used by Nelder-Mead.
  struct coord_t
  {
    double x, y;
  };

  using hpolicy = boost::math::policies::policy<boost::math::policies::max_series_iterations<1000000>>;

  // Moment-based initializer from sample mean and variance.
  [[nodiscard]] static params_t init_from_moments(const std::vector<double>& x_v)
  {
    const double N = static_cast<double>(x_v.size());
    double m = 0.0, v = 0.0;
    for (double x : x_v)
      m += x;
    m /= N;
    for (double x : x_v)
      v += (x - m) * (x - m);
    v /= N;

    const double v_l = std::max(v, 1e-6);
    const double shape0 = (m > eps) ? (m * m) / v_l : 1.0;
    const double scale0 = (m > eps) ? v_l / m : std::max(m, eps);
    return {std::max(shape0, 1e-2), std::max(scale0, 1e-6)};
  }

  // Compute interpolated empirical quantiles.
  [[nodiscard]] static std::array<double, 3>
  compute_quantiles(const std::vector<double>& x_v, const std::array<double, 3>& probs)
  {
    std::vector<double> sorted = x_v;
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
  static void sort3(std::array<coord_t, 3>& S, std::array<double, 3>& F)
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

  // 2D Nelder-Mead in (log shape, log scale), takes linear-space parameters.
  // Stops when both objective spread and simplex diameter fall below tolerance.
  // The diameter check catches flat-likelihood drift that the objective only misses.
  template<typename Obj>
  [[nodiscard]] static params_t nelder_mead_bivariate(params_t p0, Obj&& obj, const Config& cfg)
  {
    auto eval = [&](const coord_t& p) { return obj(std::exp(p.x), std::exp(p.y)); };

    std::array<coord_t, 3> S = {{{std::log(p0.shape), std::log(p0.scale)},
                                 {std::log(p0.shape) + cfg.nm_step, std::log(p0.scale)},
                                 {std::log(p0.shape), std::log(p0.scale) + cfg.nm_step}}};
    std::array<double, 3> F = {eval(S[0]), eval(S[1]), eval(S[2])};
    sort3(S, F);

    for (int iter = 0; iter < cfg.max_niter; ++iter) {
      const bool converged =
        (F[2] - F[0]) < cfg.tol && std::abs(S[2].x - S[0].x) < cfg.tol && std::abs(S[2].y - S[0].y) < cfg.tol;
      if (converged && iter > 5) break;

      const coord_t c = {0.5 * (S[0].x + S[1].x), 0.5 * (S[0].y + S[1].y)};
      const coord_t r = {c.x + cfg.nm_alpha * (c.x - S[2].x), c.y + cfg.nm_alpha * (c.y - S[2].y)};
      const double fr = eval(r);

      if (fr < F[0]) {
        // Expansion
        const coord_t e = {c.x + cfg.nm_expand * (r.x - c.x), c.y + cfg.nm_expand * (r.y - c.y)};
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
        const coord_t w = use_r ? r : S[2];
        const double fw = use_r ? fr : F[2];
        const coord_t cc = {c.x + cfg.nm_rho * (w.x - c.x), c.y + cfg.nm_rho * (w.y - c.y)};
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
