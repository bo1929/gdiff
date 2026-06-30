#include "doctest/doctest.h"
#include "diststat.hpp"
#include "map.hpp"
#include "random.hpp"
#include <cmath>

namespace {

static std::pair<llh_sptr_t<double>, params_t<double>> make_diststat_params(uint64_t sample_size = 32)
{
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 0, sample_size, true);
  params.canonical = true;
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  return {llhf, params};
}

static void inject_hits_all_bins(DIM<double>& dim, uint64_t nbins, uint32_t hdist = 0, int reps = 4)
{
  for (uint64_t i = 0; i < nbins; ++i) {
    for (int r = 0; r < reps; ++r) {
      dim.aggregate_mer(hdist, i);
    }
  }
}

} // namespace

TEST_SUITE("DistanceStat<double>") {

TEST_CASE("clear_samples resets null pool") {
  auto [llhf, params] = make_diststat_params();
  DistanceStat<double> diststat(params, llhf);

  DIM<double> dim(params, llhf, 16, 16);
  inject_hits_all_bins(dim, 16);
  dim.compute_prefhistsum();

  diststat.sample_null_pool(dim, 2, 0);
  CHECK(diststat.samples().size() > 0);
  diststat.clear_samples();
  CHECK(diststat.samples().empty());
}

TEST_CASE("sample_null_pool collects windows from hit-rich query") {
  seed = 7;
  init_thread_rng(0);

  auto [llhf, params] = make_diststat_params(20);
  DistanceStat<double> diststat(params, llhf);

  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();

  diststat.sample_null_pool(dim, 2, 0);
  CHECK(diststat.samples().size() >= GammaModel::min_nsamples);
  for (const auto& s : diststat.samples()) {
    CHECK(std::isfinite(s.d));
    CHECK(std::isfinite(s.I));
    CHECK(s.bix == 0);
    CHECK(s.bin_iv.a < s.bin_iv.b);
  }
}

TEST_CASE("test_significance scores a single record") {
  seed = 11;
  init_thread_rng(0);

  auto [llhf, params] = make_diststat_params(20);
  DistanceStat<double> diststat(params, llhf);

  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();
  diststat.sample_null_pool(dim, 2, 0);
  REQUIRE(diststat.samples().size() >= GammaModel::min_nsamples);

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(4, 8, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(d_obs);

  record_t r(0, 32, interval_t{1, 20}, interval_t{5, 9}, false, d_obs, I_obs, 0);
  r.d_diff = -0.05;
  r.d_q = 0.2;

  CHECK(diststat.test_significance(r, params.sample_size, "q0"));
  CHECK(std::isfinite(r.percentile));
  CHECK(r.percentile >= 0.0);
  CHECK(r.percentile <= 1.0);

  const bool two_sided = !std::isnan(r.d_diff) && (r.is_rc == (r.d_diff > 0.0));
  CHECK(two_sided);
}

TEST_CASE("benjamini_hochberg_correction assigns qvalues per strand") {
  seed = 11;
  init_thread_rng(0);

  auto [llhf, params] = make_diststat_params(20);
  DistanceStat<double> diststat(params, llhf);

  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();
  diststat.sample_null_pool(dim, 2, 0);
  REQUIRE(diststat.samples().size() >= GammaModel::min_nsamples);

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(4, 8, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(d_obs);

  vec<record_t> records;
  records.emplace_back(0, 32, interval_t{1, 20}, interval_t{5, 9}, false, d_obs, I_obs, 0);
  records.emplace_back(0, 32, interval_t{10, 25}, interval_t{12, 16}, false, d_obs, I_obs, 0);
  records.emplace_back(0, 32, interval_t{1, 20}, interval_t{6, 10}, true, d_obs, I_obs, 0);
  for (auto& r : records) {
    r.d_diff = -0.05;
    r.d_q = 0.2;
    CHECK(diststat.test_significance(r, params.sample_size, "q0"));
  }

  diststat.benjamini_hochberg_correction(records);

  for (const auto& r : records) {
    CHECK(std::isfinite(r.qvalue));
    CHECK(r.qvalue >= r.percentile);
    CHECK(r.qvalue <= 1.0);
  }
}

TEST_CASE("overlapping null windows on same query are excluded from scoring") {
  seed = 19;
  init_thread_rng(0);

  auto [llhf, params] = make_diststat_params(20);
  DistanceStat<double> diststat(params, llhf);

  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();
  diststat.sample_null_pool(dim, 2, 0);
  REQUIRE(diststat.samples().size() >= GammaModel::min_nsamples);

  const sample_t& anchor = diststat.samples().front();
  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(anchor.bin_iv.a - 1, anchor.bin_iv.b - 1, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(d_obs);

  record_t r(anchor.bix, 20, interval_t{1, 20}, anchor.bin_iv, false, d_obs, I_obs, 0);
  r.d_diff = 0.1;
  r.d_q = 0.15;

  CHECK(diststat.test_significance(r, params.sample_size, "q0"));
  CHECK(std::isfinite(r.percentile));
}

TEST_CASE("sample_null_pool skips when window exceeds query bins") {
  auto [llhf, params] = make_diststat_params(20);
  DistanceStat<double> diststat(params, llhf);

  DIM<double> dim(params, llhf, 2, 2);
  inject_hits_all_bins(dim, 2, 0, 2);
  dim.compute_prefhistsum();

  diststat.sample_null_pool(dim, 2, 0); // win_len = tau_eff + 1 = 3 > nbins
  CHECK(diststat.samples().empty());
}

TEST_CASE("test_significance fails when null pool is too small") {
  auto [llhf, params] = make_diststat_params(20);
  DistanceStat<double> diststat(params, llhf);

  record_t r(0, 32, interval_t{1, 20}, interval_t{5, 9}, false, 0.1, 10.0, 0);
  r.d_diff = -0.05;
  r.d_q = 0.2;

  CHECK_FALSE(diststat.test_significance(r, params.sample_size, "q0"));
  CHECK(std::isnan(r.percentile));
  CHECK(std::isnan(r.fold));
}

TEST_CASE("canonical records use one-sided test when d_diff is NaN") {
  seed = 11;
  init_thread_rng(0);

  auto [llhf, params] = make_diststat_params(20);
  DistanceStat<double> diststat(params, llhf);

  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();
  diststat.sample_null_pool(dim, 2, 0);
  REQUIRE(diststat.samples().size() >= GammaModel::min_nsamples);

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(4, 8, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(d_obs);

  record_t r(0, 32, interval_t{1, 20}, interval_t{5, 9}, false, d_obs, I_obs, 0);
  r.d_diff = nanx();

  CHECK(diststat.test_significance(r, params.sample_size, "q0"));
  CHECK(std::isfinite(r.percentile));
  CHECK(r.percentile >= 0.0);
  CHECK(r.percentile <= 1.0);
  const bool two_sided = !std::isnan(r.d_diff) && (r.is_rc == (r.d_diff > 0.0));
  CHECK_FALSE(two_sided);
}

TEST_CASE("reference strand d_diff zero uses two-sided test") {
  seed = 11;
  init_thread_rng(0);

  auto [llhf, params] = make_diststat_params(20);
  DistanceStat<double> diststat(params, llhf);

  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();
  diststat.sample_null_pool(dim, 2, 0);
  REQUIRE(diststat.samples().size() >= GammaModel::min_nsamples);

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(4, 8, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(d_obs);

  record_t r(0, 32, interval_t{1, 20}, interval_t{5, 9}, false, d_obs, I_obs, 0);
  r.d_diff = 0.0;

  CHECK(diststat.test_significance(r, params.sample_size, "q0"));
  CHECK(std::isfinite(r.percentile));
  CHECK((!std::isnan(r.d_diff) && (r.is_rc == (r.d_diff > 0.0))));
}

TEST_CASE("query strand uses one-sided test") {
  seed = 11;
  init_thread_rng(0);

  auto [llhf, params] = make_diststat_params(20);
  DistanceStat<double> diststat(params, llhf);

  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();
  diststat.sample_null_pool(dim, 2, 0);
  REQUIRE(diststat.samples().size() >= GammaModel::min_nsamples);

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(4, 8, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(d_obs);

  record_t r(0, 32, interval_t{1, 20}, interval_t{5, 9}, true, d_obs, I_obs, 0);
  r.d_diff = -0.1;
  r.d_q = 0.15;

  CHECK(diststat.test_significance(r, params.sample_size, "q0"));
  CHECK(std::isfinite(r.percentile));
  const bool two_sided = !std::isnan(r.d_diff) && (r.is_rc == (r.d_diff > 0.0));
  CHECK_FALSE(two_sided);
}

TEST_CASE("benjamini_hochberg_correction with canonical-only records") {
  vec<record_t> records;
  records.emplace_back(0, 100, interval_t{1, 50}, interval_t{1, 5}, false, 0.1, 10.0, 0);
  records.emplace_back(0, 100, interval_t{51, 100}, interval_t{6, 10}, false, 0.12, 10.0, 0);
  records[0].percentile = 0.05;
  records[1].percentile = 0.10;

  DistanceStat<double> diststat(params_t<double>(1, 0.1, 4, 2, 33.0, 0, 20, true),
                                std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1));
  diststat.benjamini_hochberg_correction(records);

  CHECK(records[0].qvalue == doctest::Approx(0.10));
  CHECK(records[1].qvalue == doctest::Approx(0.10));
}

TEST_CASE("benjamini_hochberg_correction leaves NaN qvalues untouched") {
  vec<record_t> records;
  records.emplace_back(0, 100, interval_t{1, 50}, interval_t{1, 5}, false, 0.1, 10.0, 0);
  records.emplace_back(0, 100, interval_t{51, 100}, interval_t{6, 10}, false, 0.12, 10.0, 0);
  records[0].percentile = 0.05;
  records[1].percentile = nanx();

  DistanceStat<double> diststat(params_t<double>(1, 0.1, 4, 2, 33.0, 0, 20, true),
                                std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1));
  diststat.benjamini_hochberg_correction(records);

  CHECK(records[0].qvalue == doctest::Approx(0.05));
  CHECK(std::isnan(records[1].qvalue));
}

} // TEST_SUITE
