#include "doctest/doctest.h"
#include "dist.hpp"
#include "map.hpp"
#include <algorithm>
#include <boost/math/distributions/gamma.hpp>
#include <cmath>
#include <random>

TEST_SUITE("gamma significance") {

TEST_CASE("apply_gamma_significance reads the p-value off the fitted gamma") {
  gamma_fit_t fit;
  fit.shape = 4.0;
  fit.scale = 0.02;
  const boost::math::gamma_distribution<double> gamma(fit.shape, fit.scale);

  const double median = boost::math::median(gamma);
  record_t med(0, 32, interval_t{1, 20}, false, median, 0);
  apply_gamma_significance(med, fit);
  CHECK(med.p_gamma == doctest::Approx(1.0)); // the median sits in both tails

  // Both tails land on the same scale: an unusually similar record and an unusually divergent
  // one are scored alike.
  record_t low(0, 32, interval_t{1, 20}, false, boost::math::quantile(gamma, 1e-4), 0);
  apply_gamma_significance(low, fit);
  CHECK(low.p_gamma == doctest::Approx(2e-4).epsilon(0.02));
  record_t high(0, 32, interval_t{1, 20}, false, boost::math::quantile(gamma, 1.0 - 1e-4), 0);
  apply_gamma_significance(high, fit);
  CHECK(high.p_gamma == doctest::Approx(2e-4).epsilon(0.02));

  // No usable fit, or no distance: the row stays untested.
  const gamma_fit_t no_fit;
  record_t untested(0, 32, interval_t{1, 20}, false, 0.1, 0);
  apply_gamma_significance(untested, no_fit);
  CHECK(std::isnan(untested.p_gamma));
  record_t unmapped(0, 32, interval_t{1, 20}, false, nanx(), 0);
  apply_gamma_significance(unmapped, fit);
  CHECK(std::isnan(unmapped.p_gamma));
}

} // TEST_SUITE

TEST_SUITE("make_window_plan") {

TEST_CASE("window span rounds up to whole bins") {
  CHECK(get_nwinmers(500, 0) == 500);
  CHECK(get_nwinmers(1, 0) == 1);
  CHECK(get_nwinmers(10, 1) == 10);
  CHECK(get_nwinmers(9, 1) == 10); // ceil(9 / 2) * 2
  CHECK(get_nwinmers(5, 2) == 8);  // ceil(5 / 4) * 4
  CHECK(get_nwinmers(0, 2) == 4);  // never shorter than one bin
  CHECK(get_nwinmers(50, 2) == 52); // the bin_shift=2 case below, directly
}

TEST_CASE("one global sample is split across sources, skipping short ones") {
  constexpr uint64_t k = 5;
  constexpr uint64_t tau = 10;
  const vec<uint64_t> source_lens_v{100, 10, 200}; // 10 bases cannot hold k + tau - 1 = 14
  std::mt19937 rng(7);

  const window_plan_t wp = make_window_plan(source_lens_v, k, tau, 0, 20, rng);

  CHECK(wp.nwinmers == tau);
  CHECK(wp.nwins() == 20);
  REQUIRE(wp.sources_v.size() == 2);
  CHECK(wp.sources_v[0].bix == 0);
  CHECK(wp.sources_v[1].bix == 2);
  CHECK(wp.sources_v[0].enmers == 100 - k + 1);
  CHECK(wp.sources_v[1].enmers == 200 - k + 1);

  vec<uint64_t> all_v;
  for (const window_plan_t::source_t& src : wp.sources_v) {
    CHECK(std::is_sorted(src.starts_v.begin(), src.starts_v.end()));
    for (const uint64_t start : src.starts_v) {
      CHECK(start + wp.nwinmers <= src.enmers);
      all_v.push_back(start);
    }
  }
  std::sort(all_v.begin(), all_v.end());
  CHECK(std::adjacent_find(all_v.begin(), all_v.end()) == all_v.end()); // no window drawn twice
}

TEST_CASE("bin_shift quantises the draw to whole bins") {
  constexpr uint64_t k = 5;
  constexpr uint64_t tau = 9;
  std::mt19937 rng(5);
  // bin_size 2 rounds tau=9 up to a 10-mer window, so 44 whole-bin positions fit in 96.
  const window_plan_t wp = make_window_plan(vec<uint64_t>{100}, k, tau, 1, 20, rng);

  CHECK(wp.nwinmers == 10);
  REQUIRE(wp.sources_v.size() == 1);
  CHECK(wp.nwins() == 20);
  for (const uint64_t start : wp.sources_v[0].starts_v)
    CHECK(start <= 43);
}

TEST_CASE("sources without room yield an empty plan") {
  std::mt19937 rng(3);
  const window_plan_t wp = make_window_plan(vec<uint64_t>{4, 13}, 5, 10, 0, 20, rng);
  CHECK(wp.nwinmers == 10);
  CHECK(wp.nwins() == 0);
  CHECK(wp.sources_v.empty());
}

} // TEST_SUITE

TEST_SUITE("dist helpers") {

TEST_CASE("select_strand_distance prefers finite lower distance") {
  {
    const auto [d, s] = select_strand_distance(0.1, 0.2);
    CHECK(d == doctest::Approx(0.1));
    CHECK(s == '+');
  }
  {
    const auto [d, s] = select_strand_distance(0.3, 0.1);
    CHECK(d == doctest::Approx(0.1));
    CHECK(s == '-');
  }
  {
    const auto [d, s] = select_strand_distance(0.2, 0.2);
    CHECK(d == doctest::Approx(0.2));
    CHECK(s == '+');
  }
  {
    const auto [d, s] = select_strand_distance(nanx(), 0.2);
    CHECK(d == doctest::Approx(0.2));
    CHECK(s == '-');
  }
  {
    const auto [d, s] = select_strand_distance(0.2, nanx());
    CHECK(d == doctest::Approx(0.2));
    CHECK(s == '+');
  }
  {
    const auto [d, s] = select_strand_distance(nanx(), nanx());
    CHECK(std::isnan(d));
    CHECK(s == '.');
  }
}

TEST_CASE("linear_quantile on sorted values") {
  CHECK(std::isnan(linear_quantile({}, 0.5)));
  CHECK(linear_quantile({0.4}, 0.5) == doctest::Approx(0.4));
  const vec<double> v{0.0, 0.5, 1.0};
  CHECK(linear_quantile(v, 0.0) == doctest::Approx(0.0));
  CHECK(linear_quantile(v, 1.0) == doctest::Approx(1.0));
  CHECK(linear_quantile(v, 0.5) == doctest::Approx(0.5));
  CHECK(linear_quantile({0.0, 10.0}, 0.25) == doctest::Approx(2.5));
}

TEST_CASE("validate_binning rejects oversized bins") {
  CHECK(validate_binning(0, 10));
  CHECK(validate_binning(3, 8));
  CHECK_FALSE(validate_binning(4, 8)); // bin_size=16 > tau
  CHECK_FALSE(validate_binning(17, 100));
}

} // TEST_SUITE

TEST_SUITE("gamma background thresholds") {

TEST_CASE("method of moments recovers a known gamma") {
  const double shape = 4.0;
  const double scale = 0.02;
  const boost::math::gamma_distribution<double> gamma(shape, scale);
  const size_t n = 20000;
  vec<double> sample;
  sample.reserve(n);
  for (size_t i = 0; i < n; ++i)
    sample.push_back(boost::math::quantile(gamma, (static_cast<double>(i) + 0.5) / static_cast<double>(n)));

  gamma_fit_t fit;
  REQUIRE(fit_gamma_mom(sample, fit, 8));
  CHECK(fit.shape == doctest::Approx(shape).epsilon(0.02));
  CHECK(fit.scale == doctest::Approx(scale).epsilon(0.02));
}

TEST_CASE("a short sample fails") {
  const vec<double> sample{0.1, 0.2, 0.3};
  gamma_fit_t fit;
  CHECK_FALSE(fit_background_gamma(sample, fit, 8));

  // The caller's minimum - a tenth of the pool - gates the gamma fit as well: twelve values
  // fit when ten are asked for, and do not when more are.
  vec<double> twelve;
  for (int i = 0; i < 12; ++i)
    twelve.push_back(0.02 + 0.02 * static_cast<double>(i));
  CHECK(fit_background_gamma(twelve, fit, 10));
  CHECK(std::isfinite(fit.shape));
  CHECK_FALSE(fit_background_gamma(twelve, fit, 20));
}

TEST_CASE("the fit carries the floor spike separately and keeps the gamma on the rest") {
  // A pool with a fifth of its windows pinned at the estimable floor: windows with the same
  // histogram get the same MLE, so that is an exact point mass, and a gamma fitted straight to the
  // pool collapses onto it (shape 1.18 here instead of 2). The hurdle counts the spike and fits the
  // windows above it, which recovers both the share and the shape.
  const double shape = 2.0;
  const double scale = 0.165;
  const double floor_d = 2e-7;
  const boost::math::gamma_distribution<double> gamma(shape, scale);
  const size_t n = 500;
  const size_t n_floor = n / 5;

  vec<double> pool(n_floor, floor_d);
  for (size_t i = 0; i < n - n_floor; ++i)
    pool.push_back(boost::math::quantile(gamma, (static_cast<double>(i) + 0.5) / static_cast<double>(n - n_floor)));
  std::sort(pool.begin(), pool.end());

  gamma_fit_t bare;
  REQUIRE(fit_gamma_mom(pool, bare, 8));
  CHECK(bare.shape < 1.5); // the spike drags the moment fit down

  gamma_fit_t fit;
  REQUIRE(fit_background_gamma(pool, fit, 8));
  CHECK(fit.p0 == doctest::Approx(0.2).epsilon(0.02));
  CHECK(fit.d_floor == doctest::Approx(floor_d));
  CHECK(fit.shape == doctest::Approx(shape).epsilon(0.05));
  CHECK(fit.via_ml);

  // The hurdle is a distribution over the whole domain, so it also scores records: one at the floor
  // sits in the spike's mass.
  record_t r(0, 100, interval_t{1, 50}, false, floor_d, 0);
  apply_gamma_significance(r, fit);
  CHECK(r.p_gamma == doctest::Approx(2.0 * fit.p0).epsilon(0.02));
}

TEST_CASE("the hurdle anchors at the physical floor when the sample minimum is above it")
{
  const double floor_d = 2e-7;
  vec<double> sample;
  for (int i = 0; i < 500; ++i)
    sample.push_back(floor_d + 0.001 + 0.0001 * static_cast<double>(i));
  gamma_fit_t fit;
  REQUIRE(fit_background_gamma(sample, fit, 10, floor_d));
  CHECK(fit.d_floor == doctest::Approx(floor_d));
  CHECK(fit.p0 == doctest::Approx(0.0));
  CHECK(fit.shape > 0.0);

  vec<double> th, pi;
  REQUIRE(level_thresholds(sample, {0.1}, th, pi, 10, nullptr, floor_d));
  REQUIRE(!th.empty());
  CHECK(th.front() > floor_d);
}

TEST_CASE("map skips gamma thresholds outside the distance model domain")
{
  vec<double> pool;
  for (size_t i = 0; i < 1000; ++i)
    pool.push_back(0.01 + 0.0008 * static_cast<double>(i));

  gamma_fit_t fit;
  vec<double> thresholds, pi;
  const vec<double> levels{0.01, 0.001, 0.0001};
  REQUIRE(level_thresholds(pool, levels, thresholds, pi, 10, &fit, 1e-9));
  const boost::math::gamma_distribution<double> model(fit.shape, fit.scale);
  REQUIRE(boost::math::quantile(model, 1.0 - levels.back()) > UB);
  REQUIRE(!thresholds.empty());
  REQUIRE(thresholds.size() == pi.size());
  for (size_t i = 0; i < thresholds.size(); ++i) {
    CHECK(is_valid_distance(thresholds[i]));
    // Only the lower tails lie inside the model domain for this pool.
    CHECK(thresholds[i] == doctest::Approx(boost::math::quantile(model, pi[i])));
  }
}

TEST_CASE("a level inside the floor mass loses only its lower threshold, and the report says so")
{
  // A fifth of the pool sits at the floor, so p0 = 0.2: pi = 0.1 asks for a lower quantile inside the
  // spike, pi = 0.3 does not.
  const double floor_d = 2e-7;
  const boost::math::gamma_distribution<double> gamma(2.0, 0.05);
  const size_t n = 500;
  vec<double> pool(n / 5, floor_d);
  for (size_t i = 0; i < n - (n / 5); ++i)
    pool.push_back(boost::math::quantile(gamma, (static_cast<double>(i) + 0.5) / static_cast<double>(n - (n / 5))));

  vec<double> th, pi;
  level_report_t report;
  REQUIRE(level_thresholds(pool, {0.1, 0.3}, th, pi, 10, nullptr, floor_d, &report));
  REQUIRE(report.drops_v.size() == 1);
  CHECK(report.drops_v[0].low);
  CHECK(report.drops_v[0].pi == doctest::Approx(0.1));
  CHECK(report.drops_v[0].reason.find("floor mass") != std::string::npos);
  CHECK(report.skip.empty());
  // Both upper thresholds and the 0.3 lower one are placed: pi 0.3, 0.1 above and 0.3 below.
  REQUIRE(th.size() == 3);
  CHECK(pi.front() == doctest::Approx(0.3));
  CHECK(std::is_sorted(th.begin(), th.end()));
}

TEST_CASE("a level is the tail probability on each side, undivided")
{
  const boost::math::gamma_distribution<double> gamma(3.0, 0.02);
  vec<double> pool;
  for (size_t i = 0; i < 1000; ++i)
    pool.push_back(boost::math::quantile(gamma, (static_cast<double>(i) + 0.5) / 1000.0));
  vec<double> th, pi;
  gamma_fit_t fit;
  REQUIRE(level_thresholds(pool, {0.05}, th, pi, 10, &fit, 1e-9));
  REQUIRE(th.size() == 2);
  const boost::math::gamma_distribution<double> model(fit.shape, fit.scale);
  CHECK(th[0] == doctest::Approx(boost::math::quantile(model, 0.05)));
  CHECK(th[1] == doctest::Approx(boost::math::quantile(model, 0.95)));
  CHECK(pi[0] == doctest::Approx(0.05));
  CHECK(pi[1] == doctest::Approx(0.05));
}

TEST_CASE("a rejected fit says why it was rejected")
{
  vec<double> th, pi;
  level_report_t report;
  CHECK_FALSE(level_thresholds(vec<double>(5, 0.1), {0.1}, th, pi, 10, nullptr, nanx(), &report));
  CHECK(report.skip.find("at least 10") != std::string::npos);

  vec<double> spike(100, 1e-9);
  for (int i = 0; i < 5; ++i)
    spike.push_back(0.01 + 0.01 * static_cast<double>(i));
  CHECK_FALSE(level_thresholds(spike, {0.1}, th, pi, 20, nullptr, 1e-9, &report));
  CHECK(report.skip.find("5 windows above the floor") != std::string::npos);
}

} // TEST_SUITE

TEST_SUITE("hdist threshold")
{
  TEST_CASE("the >20 Mbp cap applies only to an unset --hdist-th")
  {
    const uint64_t small = 19ull * 1000 * 1000;
    const uint64_t large = 20ull * 1000 * 1000;

    // Unset: a large input is capped to 2, a small one keeps the default.
    CHECK(hdist_th_for(large, std::nullopt) == 2);
    CHECK(hdist_th_for(small, std::nullopt) == 3);

    // Explicit: the requested value is honoured whatever the length.
    CHECK(hdist_th_for(large, 3) == 3);
    CHECK(hdist_th_for(large, 5) == 5);

    // An explicit request at or below the cap is unchanged.
    CHECK(hdist_th_for(large, 2) == 2);
    CHECK(hdist_th_for(large, 1) == 1);
  }
}
