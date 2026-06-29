#include "doctest/doctest.h"
#include "gamma.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

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

TEST_SUITE("GammaModel::fit_from_samples") {

TEST_CASE("recovers known Gamma(2, 0.05) parameters from raw draws") {
  constexpr double shape_true = 2.0;
  constexpr double scale_true = 0.05;
  std::mt19937_64 rng(42);
  std::gamma_distribution<double> dist(shape_true, scale_true);
  std::vector<double> samples(5500);
  for (auto& s : samples) s = dist(rng);

  auto result = GammaModel::fit_from_samples(samples, strict_fit_from_samples_cfg());
  check_gamma_recovery(shape_true, scale_true, result, 0.15);
}

TEST_CASE("recovers Gamma(5, 0.02) parameters from raw draws") {
  constexpr double shape_true = 5.0;
  constexpr double scale_true = 0.02;
  std::mt19937_64 rng(123);
  std::gamma_distribution<double> dist(shape_true, scale_true);
  std::vector<double> samples(6500);
  for (auto& s : samples) s = dist(rng);

  auto result = GammaModel::fit_from_samples(samples, strict_fit_from_samples_cfg());
  check_gamma_recovery(shape_true, scale_true, result, 0.15);
}

TEST_CASE("empty input returns safe default (1, 1)") {
  std::vector<double> empty;
  auto r = GammaModel::fit_from_samples(empty);
  CHECK(r.shape == 1.0);
  CHECK(r.scale == 1.0);
}

TEST_CASE("constant samples produce a finite fit") {
  std::vector<double> constant(100, 0.05);
  auto r = GammaModel::fit_from_samples(constant);
  CHECK(r.shape > 0.0);
  CHECK(r.scale > 0.0);
  CHECK(std::isfinite(r.shape));
  CHECK(std::isfinite(r.scale));
}

TEST_CASE("tiny input converges without crashing") {
  std::vector<double> one = {0.1};
  auto r1 = GammaModel::fit_from_samples(one);
  CHECK(r1.shape > 0.0);
  CHECK(r1.scale > 0.0);

  std::vector<double> two = {0.1, 0.2};
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
  std::vector<double> samples(5500);
  for (auto& s : samples) s = dist(rng);

  auto r = GammaModel::fit_from_samples(samples, cfg);
  check_gamma_recovery(shape_true, scale_true, r, 0.18);
}

TEST_CASE("same multiset of samples yields same fit regardless of order") {
  std::mt19937_64 rng(99);
  std::gamma_distribution<double> dist(3.0, 0.04);
  std::vector<double> samples(320);
  for (auto& s : samples) s = dist(rng);

  std::vector<double> perm = samples;
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
  std::vector<double> samples(5500);
  for (auto& s : samples) s = dist(rng);

  auto r = GammaModel::fit_from_samples(samples, strict_fit_from_samples_cfg());
  check_gamma_recovery(shape_true, scale_true, r, 0.18);
}

} // TEST_SUITE

TEST_SUITE("GammaModel::validate_params") {

TEST_CASE("accepts positive finite shape and scale") {
  CHECK(GammaModel::validate_params({2.0, 0.05}));
}

TEST_CASE("rejects non-positive or non-finite parameters") {
  CHECK_FALSE(GammaModel::validate_params({0.0, 0.05}));
  CHECK_FALSE(GammaModel::validate_params({2.0, 0.0}));
  CHECK_FALSE(GammaModel::validate_params({std::numeric_limits<double>::quiet_NaN(), 0.05}));
  CHECK_FALSE(GammaModel::validate_params({2.0, std::numeric_limits<double>::infinity()}));
}

} // TEST_SUITE

TEST_SUITE("GammaModel::score_from_samples") {

TEST_CASE("rejects too few samples and invalid observations") {
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> few(4, 0.05);
  const auto r1 = GammaModel::score_from_samples(0.1, few, 1e-5, 0.99);
  CHECK(std::isnan(r1.first));
  CHECK(std::isnan(r1.second));

  std::vector<double> many(20, 0.05);
  const auto r2 = GammaModel::score_from_samples(nan, many, 1e-5, 0.99);
  CHECK(std::isnan(r2.first));
  CHECK(std::isnan(r2.second));
}

TEST_CASE("CDF increases with observed distance") {
  std::mt19937_64 rng(11);
  std::gamma_distribution<double> dist(2.0, 0.05);
  std::vector<double> nulls(256);
  for (auto& d : nulls) d = dist(rng);

  const auto [prob_low, median] = GammaModel::score_from_samples(0.01, nulls, 1e-5, 0.99);
  const auto [prob_mid, _] = GammaModel::score_from_samples(median, nulls, 1e-5, 0.99);
  CHECK(std::isfinite(prob_low));
  CHECK(std::isfinite(median));
  CHECK(prob_low < prob_mid);
  CHECK(median > 0.03);
  CHECK(median < 0.2);
}

TEST_CASE("null draw near median gets moderate probability") {
  std::mt19937_64 rng(22);
  std::gamma_distribution<double> dist(3.0, 0.04);
  std::vector<double> nulls(64);
  for (auto& d : nulls) d = dist(rng);

  const double target = dist(rng);
  const auto [prob, median] = GammaModel::score_from_samples(target, nulls, 1e-5, 0.99);
  CHECK(std::isfinite(prob));
  CHECK(std::isfinite(median));
  CHECK(prob > 0.01);
  CHECK(prob < 0.99);
}

} // TEST_SUITE

TEST_SUITE("GammaModel::pdf and cdf") {

TEST_CASE("pdf and cdf are positive at interior points") {
  const double shape = 2.0;
  const double scale = 0.05;
  const double x = 0.1;
  CHECK(GammaModel::pdf(x, shape, scale) > 0.0);
  const double p = GammaModel::cdf(x, shape, scale);
  CHECK(p > 0.0);
  CHECK(p < 1.0);
}

TEST_CASE("pdf returns 0 for non-positive x") {
  CHECK(GammaModel::pdf(0.0, 2.0, 0.05) == 0.0);
  CHECK(GammaModel::pdf(-1.0, 2.0, 0.05) == 0.0);
}

TEST_CASE("cdf returns 0 for non-positive x and NaN for invalid parameters") {
  CHECK(GammaModel::cdf(0.0, 2.0, 0.05) == 0.0);
  CHECK(std::isnan(GammaModel::cdf(0.1, 0.0, 0.05)));
  CHECK(std::isnan(GammaModel::cdf(0.1, 2.0, 0.0)));
}

} // TEST_SUITE
