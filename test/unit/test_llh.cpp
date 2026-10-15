#include "doctest/doctest.h"
#include "llh.hpp"
#include <cmath>
#include <limits>
#include <numeric>

static constexpr double tol = 1e-6;

TEST_SUITE("LLH binomial coefficients") {

TEST_CASE("binom_coef_k for k=27") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  CHECK(llh.binom_coef_k[0] == 1);
  CHECK(llh.binom_coef_k[1] == 27);
  CHECK(llh.binom_coef_k[2] == 351);
  CHECK(llh.binom_coef_k[27] == 1);
  for (uint32_t d = 0; d <= 27; ++d)
    CHECK(llh.binom_coef_k[d] == llh.binom_coef_k[27 - d]);
}

TEST_CASE("binom_coef_hnk[0] is always 0") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  CHECK(llh.binom_coef_hnk[0] == 0);
}

TEST_CASE("binom_coef_hnk relation: C(k,d) - C(k-h,d)") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  CHECK(llh.binom_coef_hnk[1] == 11);
  CHECK(llh.binom_coef_hnk[2] == 231);
}

} // TEST_SUITE

TEST_SUITE("LLH<double>") {

TEST_CASE("prob_hit at zero distance equals rho for d=0") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  CHECK(llh.prob_hit(0.0, 0) == doctest::Approx(0.5));
}

TEST_CASE("prob_miss at zero distance equals 1-rho") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  CHECK(llh.prob_miss(0.0) == doctest::Approx(0.5));
}

TEST_CASE("prob_miss is between 0 and 1 and varies with distance") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  for (double D : {0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5}) {
    double pm = llh.prob_miss(D);
    CHECK(pm >= 0.0);
    CHECK(pm <= 1.0);
  }
  CHECK(std::abs(llh.prob_miss(0.001) - llh.prob_miss(0.5)) > 1e-6);
}

TEST_CASE("prob_elude + prob_collide = 1") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  for (uint32_t d = 0; d <= 4; ++d) {
    double sum = llh.prob_elude(d) + llh.prob_collide(d);
    CHECK(sum == doctest::Approx(1.0).epsilon(tol));
  }
}

TEST_CASE("hit and miss probabilities sum to one") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  for (const double D : {0.0, 0.001, 0.05, 0.2, 0.5}) {
    double p = llh.prob_miss(D);
    for (uint32_t d = 0; d <= 4; ++d)
      p += llh.prob_hit(D, d);
    CHECK(p == doctest::Approx(1.0).epsilon(tol));
  }
}

TEST_CASE("negative log-likelihood is minimized near true distance") {
  const double D_true = 0.05;
  const uint32_t k = 27, h = 11, hdist_th = 4;
  const double rho = 0.5;
  LLH<double> llh(k, h, rho, hdist_th, D_true);

  const uint64_t N = 100000;
  std::vector<uint64_t> counts(hdist_th + 1, 0);
  for (uint32_t d = 0; d <= hdist_th; ++d) {
    counts[d] = static_cast<uint64_t>(N * llh.prob_hit(D_true, d));
  }
  uint64_t misses = static_cast<uint64_t>(N * llh.prob_miss(D_true));

  double nll_true = llh.nll(D_true, counts.data(), misses);
  double nll_low  = llh.nll(0.001, counts.data(), misses);
  double nll_high = llh.nll(0.3, counts.data(), misses);

  CHECK(nll_true < nll_low);
  CHECK(nll_true < nll_high);
  CHECK(llh.mle(counts.data(), misses) == doctest::Approx(D_true).epsilon(0.02));
}

TEST_CASE("negative log-likelihood differences match explicit probabilities") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  const std::vector<uint64_t> counts{50, 20, 10, 4, 1};
  const uint64_t misses = 100;
  auto explicit_nll = [&](const double D) {
    double value = -static_cast<double>(misses) * std::log(llh.prob_miss(D));
    for (uint32_t d = 0; d < counts.size(); ++d)
      value -= static_cast<double>(counts[d]) * std::log(llh.prob_hit(D, d));
    return value;
  };
  const double d0 = 0.04, d1 = 0.16;
  const double actual = llh.nll(d0, counts.data(), misses) - llh.nll(d1, counts.data(), misses);
  const double expected = explicit_nll(d0) - explicit_nll(d1);
  CHECK(actual == doctest::Approx(expected).epsilon(tol));
}

TEST_CASE("negative log-likelihood handles zero distance") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  std::vector<uint64_t> counts(5, 0);
  counts[0] = 10;
  CHECK(llh.nll(0.0, counts.data(), 5) == doctest::Approx(-5.0 * std::log(0.5)));
  counts[1] = 1;
  CHECK(std::isinf(llh.nll(0.0, counts.data(), 5)));
}

TEST_CASE("likelihood-ratio statistic rejects invalid likelihoods") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  CHECK(likelihood_ratio_statistic(3.0, 1.0) == doctest::Approx(4.0));
  CHECK(likelihood_ratio_statistic(1.0, 1.0) == doctest::Approx(0.0));
  CHECK(likelihood_ratio_statistic(1.0, 1.0 + 1e-12) == doctest::Approx(0.0));
  CHECK(std::isnan(likelihood_ratio_statistic(nan, 1.0)));
  CHECK(std::isnan(likelihood_ratio_statistic(1.0, nan)));
  CHECK(std::isnan(likelihood_ratio_statistic(inf, 1.0)));
  CHECK(std::isnan(likelihood_ratio_statistic(1.0, inf)));
}

TEST_CASE("max estimable distance is MLE of one hit at hdist_th") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1, false);
  const uint64_t n_total = 1000;
  arr<uint64_t, hdist_bound + 1> v{};
  v[4] = 1;
  const double d_ref = llh.mle(v.data(), n_total - 1);
  const double d_max = max_estimable_distance(llh, n_total);
  CHECK(d_max == doctest::Approx(d_ref).epsilon(tol));
  CHECK(is_valid_distance(d_max));
  CHECK(std::isnan(max_estimable_distance(llh, 0)));
}

TEST_CASE("lr_pvalue calibrates a likelihood-ratio statistic on chi-square, one degree of freedom") {
  CHECK(lr_pvalue(0.0) == doctest::Approx(1.0));
  CHECK(lr_pvalue(2.0) == doctest::Approx(0.1572992).epsilon(1e-6)); // erfc(1)
  CHECK(lr_pvalue(6.634897) == doctest::Approx(0.01).epsilon(1e-5)); // the --lr-th default
  CHECK(lr_pvalue(50.0) < 1e-11);
  CHECK(lr_pvalue(std::numeric_limits<double>::infinity()) == 0.0);
  CHECK(lr_pvalue(-1.0) == 1.0); // no evidence, or round-off below zero
  CHECK(std::isnan(lr_pvalue(nanx())));
}

TEST_CASE("compute_lr_ub is near zero at the max estimable distance") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1, false);
  const uint64_t n_total = 1000;
  const double d_max = max_estimable_distance(llh, n_total);
  CHECK(compute_lr_ub(llh, d_max, n_total) == doctest::Approx(0.0).epsilon(1e-8));
  const double lr_close = compute_lr_ub(llh, 0.01, n_total);
  CHECK(std::isfinite(lr_close));
  CHECK(lr_close > 1.0);
  CHECK(std::isnan(compute_lr_ub(llh, d_max, 0)));
  CHECK(std::isnan(compute_lr_ub(llh, nanx(), n_total)));
}

TEST_CASE("get_fdc sign correctness") {
  // For positive extrema, fdc should reflect the direction of the first derivative
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  // fdc_v[0] at D=0.1 is sign * (-30).
  double fdc0 = llh.get_fdc(0);
  CHECK(fdc0 < 0); // Negative because 0 < k*D for d=0
}

TEST_CASE("get_sdc values are finite") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  for (uint32_t d = 0; d <= 4; ++d) {
    CHECK(std::isfinite(llh.get_sdc(d)));
    CHECK(std::isfinite(llh.get_fdc(d)));
  }
  CHECK(std::isfinite(llh.get_sdc()));
  CHECK(std::isfinite(llh.get_fdc()));
}

} // TEST_SUITE

TEST_SUITE("LLH<cmlane_t>") {

TEST_CASE("SIMD path produces same probabilities as scalar") {
  const uint32_t k = 27, h = 11, hdist_th = 4;
  const double rho = 0.5;

  // Create scalar LLH at D=0.1
  LLH<double> llh_scalar(k, h, rho, hdist_th, 0.1);

  // Create SIMD LLH with all 8 thresholds = 0.1
  cmlane_t dths{};
  for (int i = 0; i < 8; ++i) dths[i] = 0.1;
  LLH<cmlane_t> llh_simd(k, h, rho, hdist_th, dths);

  // All 8 lanes should match the scalar fdc/sdc
  for (uint32_t d = 0; d <= hdist_th; ++d) {
    cmlane_t fdc = llh_simd.get_fdc(d);
    cmlane_t sdc = llh_simd.get_sdc(d);
    for (int i = 0; i < 8; ++i) {
      CHECK(fdc[i] == doctest::Approx(llh_scalar.get_fdc(d)).epsilon(tol));
      CHECK(sdc[i] == doctest::Approx(llh_scalar.get_sdc(d)).epsilon(tol));
    }
  }
}

TEST_CASE("SIMD path with varying thresholds") {
  const uint32_t k = 27, h = 11, hdist_th = 4;
  const double rho = 0.5;

  cmlane_t dths{};
  for (int i = 0; i < 8; ++i) dths[i] = 0.05 * (i + 1);
  LLH<cmlane_t> llh_simd(k, h, rho, hdist_th, dths);

  // Each lane should match the corresponding scalar LLH
  for (int lane = 0; lane < 8; ++lane) {
    LLH<double> llh_scalar(k, h, rho, hdist_th, dths[lane]);
    for (uint32_t d = 0; d <= hdist_th; ++d) {
      cmlane_t fdc = llh_simd.get_fdc(d);
      CHECK(fdc[lane] == doctest::Approx(llh_scalar.get_fdc(d)).epsilon(tol));
    }
  }
}

TEST_CASE("sign bitvector for mixed positive/negative thresholds") {
  const uint32_t k = 27, h = 11, hdist_th = 4;
  const double rho = 0.5;
  cmlane_t dths{};
  dths[0] = 0.1;
  dths[1] = -0.1;
  dths[2] = 0.2;
  dths[3] = -0.2;
  dths[4] = 0.3;
  dths[5] = -0.3;
  dths[6] = 0.4;
  dths[7] = -0.4;

  LLH<cmlane_t> llh(k, h, rho, hdist_th, dths);
  auto [pos_bv, neg_bv] = llh.get_sign_bv();

  // Positive: lanes 0,2,4,6 -> bits 0,2,4,6 -> 0b01010101 = 0x55
  CHECK(pos_bv == 0x55);
  // Negative: lanes 1,3,5,7 -> bits 1,3,5,7 -> 0b10101010 = 0xAA
  CHECK(neg_bv == 0xAA);
}

} // TEST_SUITE

// Reference values from an independent Python port of nll/mle for k=23, h=9, hdist=3, rho=0.5,
// n=500: max_estimable_distance=0.358841, LR roots at lr_th=6.635 give t_min=2.8844e-4,
// t_max=0.276462.
TEST_SUITE("admissible distance range") {

TEST_CASE("compute_lr_lb is zero at the floor and grows with distance") {
  LLH<double> llh(23, 9, 0.5, 3, 0.0, false);
  CHECK(compute_lr_lb(llh, d_eps, 500) < 1.0); // still essentially at the floor
  CHECK(compute_lr_lb(llh, 1e-4, 500) == doctest::Approx(2.30).epsilon(0.02));
  CHECK(compute_lr_lb(llh, 1e-3, 500) == doctest::Approx(23.0).epsilon(0.02));
  CHECK(compute_lr_lb(llh, 0.05, 500) > compute_lr_lb(llh, 0.01, 500));
  CHECK(std::isnan(compute_lr_lb(llh, 0.1, 0)));
}

TEST_CASE("compute_lr_ub falls to zero at the estimable ceiling") {
  LLH<double> llh(23, 9, 0.5, 3, 0.0, false);
  const double d_max_est = max_estimable_distance(llh, 500);
  CHECK(d_max_est == doctest::Approx(0.358841).epsilon(1e-3));
  CHECK(compute_lr_ub(llh, 0.25, 500) == doctest::Approx(13.79).epsilon(0.02));
  CHECK(compute_lr_ub(llh, 0.3, 500) < 6.635);
  CHECK(std::abs(compute_lr_ub(llh, d_max_est, 500)) < 1e-6);
}

TEST_CASE("admissible_distance_range brackets the resolvable distances") {
  LLH<double> llh(23, 9, 0.5, 3, 0.0, false);
  const a_range_t r = admissible_distance_range(llh, 500, 6.635);
  REQUIRE(r.is_ok);
  CHECK(r.t_min == doctest::Approx(2.8844e-4).epsilon(0.02));
  CHECK(r.t_max == doctest::Approx(0.276462).epsilon(1e-3));
  CHECK(r.t_min < r.t_max);
  CHECK(r.t_max < max_estimable_distance(llh, 500));

  // Inside the range both LRs clear the cut; outside it, the violated end does not.
  CHECK(compute_lr_lb(llh, r.t_min, 500) == doctest::Approx(6.635).epsilon(0.01));
  CHECK(compute_lr_ub(llh, r.t_max, 500) == doctest::Approx(6.635).epsilon(0.01));
  CHECK(compute_lr_lb(llh, 1e-5, 500) < 6.635);
  CHECK(compute_lr_ub(llh, 0.35, 500) < 6.635);
}

TEST_CASE("the p-value columns calibrate the statistics the admissible range solves") {
  // `map` reports p_lb / p_ub as these two statistics calibrated on chi-square, and the resolvable
  // bounds are exactly where they cross --lr-th, so the bounds must land on p = 0.01.
  LLH<double> llh(23, 9, 0.5, 3, 0.0, false);
  constexpr uint64_t n_win = 500;
  const a_range_t r = admissible_distance_range(llh, n_win, 6.635);
  REQUIRE(r.is_ok);
  CHECK(lr_pvalue(compute_lr_lb(llh, r.t_min, n_win)) == doctest::Approx(0.01).epsilon(1e-3));
  CHECK(lr_pvalue(compute_lr_ub(llh, r.t_max, n_win)) == doctest::Approx(0.01).epsilon(1e-3));
  // Inside the band the distance is resolved from both ends: each contrast rejects its own extreme.
  CHECK(lr_pvalue(compute_lr_lb(llh, 1.5 * r.t_min, n_win)) < 0.01);
  CHECK(lr_pvalue(compute_lr_ub(llh, 0.5 * r.t_max, n_win)) < 0.01);
  // Below the floor the lower contrast stops rejecting identity, above the ceiling the upper stops
  // rejecting the weakest detectable match, so the two values read the band from either side.
  CHECK(lr_pvalue(compute_lr_lb(llh, 0.5 * r.t_min, n_win)) > 0.01);
  CHECK(lr_pvalue(compute_lr_ub(llh, 1.5 * r.t_max, n_win)) > 0.01);
}

TEST_CASE("admissible_distance_range rejects a window nothing can be resolved in") {
  LLH<double> llh(23, 9, 0.5, 3, 0.0, false);
  CHECK_FALSE(admissible_distance_range(llh, 0, 6.635).is_ok);
  CHECK_FALSE(admissible_distance_range(llh, 500, 1e9).is_ok); // no distance can clear the cut
}

} // TEST_SUITE
