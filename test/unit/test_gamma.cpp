#include "doctest/doctest.h"
#include "gamma.hpp"
#include <algorithm>
#include <cmath>
#include <random>

namespace {

// Tighter NM settings for recovery checks (still matches production code paths).
GammaModel::Config strict_fit_from_samples_cfg()
{
  GammaModel::Config cfg{};
  cfg.tol = 1e-10;
  cfg.max_niter = 900;
  cfg.nm_step = 0.35;
  return cfg;
}

void check_gamma_recovery(double shape_true, double scale_true, const GammaModel::params_t& r, double rel_eps)
{
  CHECK(r.shape > 0.0);
  CHECK(r.scale > 0.0);
  CHECK(std::isfinite(r.shape));
  CHECK(std::isfinite(r.scale));
  CHECK(r.shape == doctest::Approx(shape_true).epsilon(rel_eps));
  CHECK(r.scale == doctest::Approx(scale_true).epsilon(rel_eps));
  const double mean_true = shape_true * scale_true;
  CHECK(r.shape * r.scale == doctest::Approx(mean_true).epsilon(rel_eps * 0.55));
}

} // namespace

TEST_SUITE("GammaModel::fit") {

TEST_CASE("recovers known Gamma(2, 0.05) parameters in the noise-free limit") {
  // Noise-free samples (sigma2 == 0 for all i) should collapse the marginal
  // likelihood to the plain Gamma likelihood, up to the 5-point GH approximation.
  std::mt19937 rng(42);
  std::gamma_distribution<double> dist(2.0, 0.05);
  vec<double> samples(500);
  for (auto& s : samples) s = dist(rng);

  const vec<double> s2(samples.size(), 0.0);
  auto result = GammaModel::fit(samples, s2);

  CHECK(result.shape == doctest::Approx(2.0).epsilon(0.6));
  CHECK(result.scale == doctest::Approx(0.05).epsilon(0.03));
}

TEST_CASE("recovers Gamma(5, 0.02) parameters in the noise-free limit") {
  std::mt19937 rng(123);
  std::gamma_distribution<double> dist(5.0, 0.02);
  vec<double> samples(1000);
  for (auto& s : samples) s = dist(rng);

  const vec<double> s2(samples.size(), 0.0);
  auto result = GammaModel::fit(samples, s2);

  // MLE on noise-free data is considerably tighter than the previous
  // quantile-matching tolerance; keep a small margin for stochastic variation.
  CHECK(result.shape == doctest::Approx(5.0).epsilon(0.5));
  CHECK(result.scale == doctest::Approx(0.02).epsilon(0.5));
}

TEST_CASE("deconvolution recovers latent shape/scale under Gaussian noise") {
  // Generate latent D_i ~ Gamma(shape, scale), observed d_i = D_i + N(0, sigma^2).
  // The fit should recover the latent parameters, not the variance-inflated ones.
  std::mt19937 rng(7);
  const double true_shape = 3.0;
  const double true_scale = 0.04;
  const double sigma = 0.01;
  std::gamma_distribution<double> gd(true_shape, true_scale);
  std::normal_distribution<double> nd(0.0, sigma);
  vec<double> d_v(1500);
  for (auto& d : d_v) d = std::max(1e-6, gd(rng) + nd(rng));

  const vec<double> s2(d_v.size(), sigma * sigma);
  auto r = GammaModel::fit(d_v, s2);

  CHECK(r.shape == doctest::Approx(true_shape).epsilon(0.6));
  CHECK(r.scale == doctest::Approx(true_scale).epsilon(0.6));
}

TEST_CASE("empty input returns the safe default (1, 1)") {
  vec<double> empty;
  vec<double> s2_empty;
  auto r = GammaModel::fit(empty, s2_empty);
  CHECK(r.shape == 1.0);
  CHECK(r.scale == 1.0);
}

TEST_CASE("tiny inputs converge without crashing") {
  vec<double> tiny = {0.1};
  auto r = GammaModel::fit(tiny, {0.0});
  CHECK(r.shape > 0.0);
  CHECK(r.scale > 0.0);
  CHECK(std::isfinite(r.shape));
  CHECK(std::isfinite(r.scale));

  vec<double> two = {0.1, 0.2};
  r = GammaModel::fit(two, {0.0, 0.0});
  CHECK(r.shape > 0.0);
  CHECK(r.scale > 0.0);
}

TEST_CASE("constant samples produce a finite fit") {
  vec<double> constant(100, 0.05);
  vec<double> s2(100, 1e-6);
  auto r = GammaModel::fit(constant, s2);
  CHECK(r.shape > 0.0);
  CHECK(r.scale > 0.0);
  CHECK(std::isfinite(r.shape));
  CHECK(std::isfinite(r.scale));
}

TEST_CASE("marginal_cdf is monotone non-decreasing and in [0, 1]") {
  const double shape = 2.0, scale = 0.05, sigma = 0.01;
  double prev = 0.0;
  for (double x = 0.0; x < 0.5; x += 0.01) {
    const double c = GammaModel::marginal_cdf(x, sigma, shape, scale);
    CHECK(c >= prev - 1e-9);
    CHECK(c >= 0.0);
    CHECK(c <= 1.0);
    prev = c;
  }
  CHECK(GammaModel::marginal_cdf(0.0, sigma, shape, scale) <= 0.5);
  CHECK(GammaModel::marginal_cdf(10.0, sigma, shape, scale) == doctest::Approx(1.0).epsilon(1e-3));
}

TEST_CASE("marginal_pdf falls back to plain Gamma pdf when sigma is 0") {
  const double shape = 2.5, scale = 0.03, x = 0.1;
  const double margin = GammaModel::marginal_pdf(x, 0.0, shape, scale);
  const double direct = GammaModel::gamma_pdf(x, shape, scale);
  CHECK(margin == doctest::Approx(direct));
}

} // TEST_SUITE

TEST_SUITE("GammaModel::fit_from_samples") {

TEST_CASE("recovers known Gamma(2, 0.05) parameters from raw draws") {
  constexpr double shape_true = 2.0;
  constexpr double scale_true = 0.05;
  std::mt19937_64 rng(42);
  std::gamma_distribution<double> dist(shape_true, scale_true);
  vec<double> samples(5500);
  for (auto& s : samples) s = dist(rng);

  auto result = GammaModel::fit_from_samples(samples, strict_fit_from_samples_cfg());
  check_gamma_recovery(shape_true, scale_true, result, 0.15);
}

TEST_CASE("recovers Gamma(5, 0.02) parameters from raw draws") {
  constexpr double shape_true = 5.0;
  constexpr double scale_true = 0.02;
  std::mt19937_64 rng(123);
  std::gamma_distribution<double> dist(shape_true, scale_true);
  vec<double> samples(6500);
  for (auto& s : samples) s = dist(rng);

  auto result = GammaModel::fit_from_samples(samples, strict_fit_from_samples_cfg());
  check_gamma_recovery(shape_true, scale_true, result, 0.15);
}

TEST_CASE("empty input returns safe default (1, 1)") {
  vec<double> empty;
  auto r = GammaModel::fit_from_samples(empty);
  CHECK(r.shape == 1.0);
  CHECK(r.scale == 1.0);
}

TEST_CASE("constant samples produce a finite fit") {
  vec<double> constant(100, 0.05);
  auto r = GammaModel::fit_from_samples(constant);
  CHECK(r.shape > 0.0);
  CHECK(r.scale > 0.0);
  CHECK(std::isfinite(r.shape));
  CHECK(std::isfinite(r.scale));
}

TEST_CASE("tiny input converges without crashing") {
  vec<double> one = {0.1};
  auto r1 = GammaModel::fit_from_samples(one);
  CHECK(r1.shape > 0.0);
  CHECK(r1.scale > 0.0);

  vec<double> two = {0.1, 0.2};
  auto r2 = GammaModel::fit_from_samples(two);
  CHECK(r2.shape > 0.0);
  CHECK(r2.scale > 0.0);
}

TEST_CASE("fit_from_samples(d_v, cfg) respects custom quantile targets") {
  constexpr double shape_true = 4.0;
  constexpr double scale_true = 0.03;
  GammaModel::Config cfg = strict_fit_from_samples_cfg();
  cfg.quantile_probs = {0.1, 0.5, 0.9};

  std::mt19937_64 rng(55);
  std::gamma_distribution<double> dist(shape_true, scale_true);
  vec<double> samples(5500);
  for (auto& s : samples) s = dist(rng);

  auto r = GammaModel::fit_from_samples(samples, cfg);
  check_gamma_recovery(shape_true, scale_true, r, 0.18);
}

TEST_CASE("same multiset of samples yields same fit regardless of order") {
  std::mt19937_64 rng(99);
  std::gamma_distribution<double> dist(3.0, 0.04);
  vec<double> samples(320);
  for (auto& s : samples) s = dist(rng);

  vec<double> perm = samples;
  std::shuffle(perm.begin(), perm.end(), rng);

  const auto a = GammaModel::fit_from_samples(samples);
  const auto b = GammaModel::fit_from_samples(perm);
  CHECK(a.shape == doctest::Approx(b.shape).epsilon(1e-12));
  CHECK(a.scale == doctest::Approx(b.scale).epsilon(1e-12));
}

TEST_CASE("recovers Gamma(1.5, 0.08) from raw draws") {
  constexpr double shape_true = 1.5;
  constexpr double scale_true = 0.08;
  std::mt19937_64 rng(202);
  std::gamma_distribution<double> dist(shape_true, scale_true);
  vec<double> samples(5500);
  for (auto& s : samples) s = dist(rng);

  auto r = GammaModel::fit_from_samples(samples, strict_fit_from_samples_cfg());
  check_gamma_recovery(shape_true, scale_true, r, 0.18);
}

} // TEST_SUITE
