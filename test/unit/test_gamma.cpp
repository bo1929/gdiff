#include "doctest/doctest.h"
#include "gamma.hpp"
#include <algorithm>
#include <cmath>
#include <random>

TEST_SUITE("GammaModel::fit") {

TEST_CASE("recovers known Gamma(2, 0.05) parameters") {
  // Generate sorted samples from Gamma(alpha=2, beta=0.05)
  std::mt19937 rng(42);
  std::gamma_distribution<double> dist(2.0, 0.05);
  vec<double> samples(500);
  for (auto& s : samples) s = dist(rng);
  std::sort(samples.begin(), samples.end());

  auto result = GammaModel::fit(samples);

  // Fitted parameters should be within ~30% of true values
  CHECK(result.alpha == doctest::Approx(2.0).epsilon(0.6));
  CHECK(result.beta  == doctest::Approx(0.05).epsilon(0.03));
}

TEST_CASE("recovers Gamma(5, 0.02) parameters") {
  std::mt19937 rng(123);
  std::gamma_distribution<double> dist(5.0, 0.02);
  vec<double> samples(1000);
  for (auto& s : samples) s = dist(rng);
  std::sort(samples.begin(), samples.end());

  auto result = GammaModel::fit(samples);
  CHECK(result.alpha == doctest::Approx(5.0).epsilon(2.0));
  CHECK(result.beta  == doctest::Approx(0.02).epsilon(0.01));
}

TEST_CASE("minimum sample size: returns default for n < 3") {
  vec<double> tiny = {0.1};
  auto r = GammaModel::fit(tiny);
  CHECK(r.alpha == 1.0);
  CHECK(r.beta == 1.0);

  vec<double> two = {0.1, 0.2};
  r = GammaModel::fit(two);
  CHECK(r.alpha == 1.0);
  CHECK(r.beta == 1.0);
}

TEST_CASE("handles near-zero samples gracefully") {
  vec<double> zeros(100, 0.0);
  auto r = GammaModel::fit(zeros);
  CHECK(r.alpha > 0);
  CHECK(r.beta > 0);
  // Should return fallback: alpha=1, beta=eps
  CHECK(r.alpha == doctest::Approx(1.0));
}

TEST_CASE("handles constant samples") {
  vec<double> constant(100, 0.05);
  auto r = GammaModel::fit(constant);
  CHECK(r.alpha > 0);
  CHECK(r.beta > 0);
  // With all samples identical, the fit should still converge
  CHECK(std::isfinite(r.alpha));
  CHECK(std::isfinite(r.beta));
}

TEST_CASE("fitted CDF at empirical quantiles is close to target") {
  std::mt19937 rng(77);
  std::gamma_distribution<double> dist(3.0, 0.03);
  vec<double> samples(500);
  for (auto& s : samples) s = dist(rng);
  std::sort(samples.begin(), samples.end());

  auto result = GammaModel::fit(samples);
  boost::math::gamma_distribution<double> fitted(result.alpha, result.beta);

  // Check CDF at empirical quartiles
  size_t n = samples.size();
  double q25 = samples[n / 4];
  double q50 = samples[n / 2];
  double q75 = samples[3 * n / 4];

  double cdf25 = boost::math::cdf(fitted, q25);
  double cdf50 = boost::math::cdf(fitted, q50);
  double cdf75 = boost::math::cdf(fitted, q75);

  CHECK(cdf25 == doctest::Approx(0.25).epsilon(0.15));
  CHECK(cdf50 == doctest::Approx(0.50).epsilon(0.15));
  CHECK(cdf75 == doctest::Approx(0.75).epsilon(0.15));
}

TEST_CASE("min_nsamples constant is reasonable") {
  CHECK(GammaModel::min_nsamples >= 3);
  CHECK(GammaModel::min_nsamples <= 100);
}

} // TEST_SUITE
