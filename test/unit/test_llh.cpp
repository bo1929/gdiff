#include "doctest/doctest.h"
#include "llh.hpp"
#include <cmath>
#include <numeric>

static constexpr double tol = 1e-6;

TEST_SUITE("PLLH") {

TEST_CASE("binom_coef_k for k=27") {
  PLLH p(11, 27, 0.5, 4);
  // C(27,0)=1, C(27,1)=27, C(27,2)=351
  CHECK(p.binom_coef_k[0] == 1);
  CHECK(p.binom_coef_k[1] == 27);
  CHECK(p.binom_coef_k[2] == 351);
  CHECK(p.binom_coef_k[27] == 1);
  // Symmetry: C(27,k) == C(27, 27-k)
  for (uint32_t d = 0; d <= 27; ++d)
    CHECK(p.binom_coef_k[d] == p.binom_coef_k[27 - d]);
}

TEST_CASE("binom_coef_hnk[0] is always 0") {
  PLLH p(11, 27, 0.5, 4);
  CHECK(p.binom_coef_hnk[0] == 0);
}

TEST_CASE("binom_coef_hnk relation: C(k,d) - C(k-h,d)") {
  // For h=11, k=27, nh=16
  // binom_coef_hnk[d] = C(27,d) - C(16,d)
  PLLH p(11, 27, 0.5, 4);
  // C(16,1)=16, so binom_coef_hnk[1] = 27 - 16 = 11
  CHECK(p.binom_coef_hnk[1] == 11);
  // C(16,2)=120, so binom_coef_hnk[2] = 351 - 120 = 231
  CHECK(p.binom_coef_hnk[2] == 231);
}

} // TEST_SUITE

TEST_SUITE("LLH<double>") {

TEST_CASE("prob_hit at D=0 should be rho for d=0") {
  // At D=0: prob_mutate(0, 0) = (1-0)^k * 0^0 * C(k,0) = 1
  // prob_collide(0) = binom_coef_hnk[0]/binom_coef_k[0] = 0
  // So prob_hit(0, 0) = rho * 0 * 1 = 0
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  CHECK(llh.prob_hit(1e-10, 0) == doctest::Approx(0.0).epsilon(1e-8));
}

TEST_CASE("prob_miss at D=0 should be close to 1-rho") {
  // At D~0: most k-mers are exact matches that elude (d=0, prob_elude(0)=1)
  // prob_miss ~ rho * 1 + (1-rho) = 1
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  double pm = llh.prob_miss(1e-10);
  CHECK(pm == doctest::Approx(1.0).epsilon(1e-4));
}

TEST_CASE("prob_miss is between 0 and 1 and varies with distance") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  for (double D : {0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5}) {
    double pm = llh.prob_miss(D);
    CHECK(pm >= 0.0);
    CHECK(pm <= 1.0);
  }
  // At very low distance, miss probability should be high (most k-mers elude the sketch)
  // At medium distance, the behavior depends on rho, k, h, hdist_th
  double pm_low = llh.prob_miss(0.001);
  double pm_high = llh.prob_miss(0.5);
  // Both should be valid probabilities
  CHECK(pm_low > 0);
  CHECK(pm_high > 0);
}

TEST_CASE("prob_elude + prob_collide = 1") {
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  for (uint32_t d = 0; d <= 4; ++d) {
    double sum = llh.prob_elude(d) + llh.prob_collide(d);
    CHECK(sum == doctest::Approx(1.0).epsilon(tol));
  }
}

TEST_CASE("negative log-likelihood is minimized near true distance") {
  // Create a synthetic histogram at distance D=0.05
  const double D_true = 0.05;
  const uint32_t k = 27, h = 11, hdist_th = 4;
  const double rho = 0.5;
  LLH<double> llh(k, h, rho, hdist_th, D_true);

  // Build counts that are consistent with D_true
  // Use probabilities to generate expected counts for a large number of k-mers
  const uint64_t N = 100000;
  std::vector<uint64_t> counts(hdist_th + 1, 0);
  for (uint32_t d = 0; d <= hdist_th; ++d) {
    counts[d] = static_cast<uint64_t>(N * llh.prob_hit(D_true, d));
  }
  uint64_t total_hits = 0;
  for (auto c : counts) total_hits += c;
  uint64_t misses = static_cast<uint64_t>(N * llh.prob_miss(D_true));

  llh.set_counts(counts.data(), misses);

  // NLL at true distance should be less than at significantly different distances
  double nll_true = llh(D_true);
  double nll_low  = llh(0.001);
  double nll_high = llh(0.3);

  CHECK(nll_true < nll_low);
  CHECK(nll_true < nll_high);
}

TEST_CASE("get_fdc sign correctness") {
  // For positive extrema, fdc should reflect the direction of the first derivative
  LLH<double> llh(27, 11, 0.5, 4, 0.1);
  // fdc_v[d] = sign * (d - k*D) / (D*S) where D=0.1, S=0.9
  // For d=0: sign * (0 - 27*0.1)/(0.1*0.9) = sign * (-2.7)/(0.09) = sign * (-30)
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

TEST_SUITE("LLH<cm512_t>") {

TEST_CASE("SIMD path produces same probabilities as scalar") {
  const uint32_t k = 27, h = 11, hdist_th = 4;
  const double rho = 0.5;

  // Create scalar LLH at D=0.1
  LLH<double> llh_scalar(k, h, rho, hdist_th, 0.1);

  // Create SIMD LLH with all 8 thresholds = 0.1
  cm512_t dths{};
  for (int i = 0; i < 8; ++i) dths[i] = 0.1;
  LLH<cm512_t> llh_simd(k, h, rho, hdist_th, dths);

  // All 8 lanes should match the scalar fdc/sdc
  for (uint32_t d = 0; d <= hdist_th; ++d) {
    cm512_t fdc = llh_simd.get_fdc(d);
    cm512_t sdc = llh_simd.get_sdc(d);
    for (int i = 0; i < 8; ++i) {
      CHECK(fdc[i] == doctest::Approx(llh_scalar.get_fdc(d)).epsilon(tol));
      CHECK(sdc[i] == doctest::Approx(llh_scalar.get_sdc(d)).epsilon(tol));
    }
  }
}

TEST_CASE("SIMD path with varying thresholds") {
  const uint32_t k = 27, h = 11, hdist_th = 4;
  const double rho = 0.5;

  cm512_t dths{};
  for (int i = 0; i < 8; ++i) dths[i] = 0.05 * (i + 1);
  LLH<cm512_t> llh_simd(k, h, rho, hdist_th, dths);

  // Each lane should match the corresponding scalar LLH
  for (int lane = 0; lane < 8; ++lane) {
    LLH<double> llh_scalar(k, h, rho, hdist_th, dths[lane]);
    for (uint32_t d = 0; d <= hdist_th; ++d) {
      cm512_t fdc = llh_simd.get_fdc(d);
      CHECK(fdc[lane] == doctest::Approx(llh_scalar.get_fdc(d)).epsilon(tol));
    }
  }
}

TEST_CASE("sign bitvector for mixed positive/negative thresholds") {
  const uint32_t k = 27, h = 11, hdist_th = 4;
  const double rho = 0.5;
  cm512_t dths{};
  dths[0] = 0.1;
  dths[1] = -0.1;
  dths[2] = 0.2;
  dths[3] = -0.2;
  dths[4] = 0.3;
  dths[5] = -0.3;
  dths[6] = 0.4;
  dths[7] = -0.4;

  LLH<cm512_t> llh(k, h, rho, hdist_th, dths);
  auto [pos_bv, neg_bv] = llh.get_sign_bv();

  // Positive: lanes 0,2,4,6 -> bits 0,2,4,6 -> 0b01010101 = 0x55
  CHECK(pos_bv == 0x55);
  // Negative: lanes 1,3,5,7 -> bits 1,3,5,7 -> 0b10101010 = 0xAA
  CHECK(neg_bv == 0xAA);
}

} // TEST_SUITE
