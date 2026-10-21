#include "doctest/doctest.h"
#include "dist.hpp"
#include "map.hpp"
#include <algorithm>
#include <boost/math/distributions/gamma.hpp>
#include <cmath>
#include <random>

TEST_SUITE("empirical significance") {

TEST_CASE("apply_empirical_significance reports the two-sided p and the fold of the pool") {
  // Ten samples 0.10 .. 0.19; 0.145 sits at or above exactly five of them, i.e. dead centre.
  vec<double> pool;
  for (uint64_t i = 0; i < 10; ++i)
    pool.push_back(0.10 + 0.01 * static_cast<double>(i));
  const double median = linear_quantile(pool, 0.5);

  record_t mid(0, 32, interval_t{1, 20}, false, 0.145, 1.0, 0);
  apply_empirical_significance(mid, pool);
  CHECK(mid.pvalue == doctest::Approx(1.0));
  CHECK(mid.fold == doctest::Approx(mid.d / median));

  record_t low(0, 32, interval_t{1, 20}, false, 0.05, 1.0, 0);
  apply_empirical_significance(low, pool);
  CHECK(low.pvalue == doctest::Approx(0.0));
  CHECK(low.fold == doctest::Approx(0.05 / median));

  // Both tails land on the same scale, so an unusually divergent record is as significant as an
  // unusually similar one - whatever strand it is on (d_diff is not consulted).
  record_t high(0, 32, interval_t{1, 20}, false, 0.30, 1.0, 0);
  high.d_diff = -0.05;
  apply_empirical_significance(high, pool);
  CHECK(high.pvalue == doctest::Approx(0.0));
}

TEST_CASE("apply_gamma_significance reads the p-value and fold off the fitted gamma") {
  gamma_fit_t fit;
  fit.shape = 4.0;
  fit.scale = 0.02;
  const boost::math::gamma_distribution<double> gamma(fit.shape, fit.scale);

  const double median = boost::math::median(gamma);
  record_t med(0, 32, interval_t{1, 20}, false, median, 1.0, 0);
  apply_gamma_significance(med, fit);
  CHECK(med.pvalue == doctest::Approx(1.0)); // the median sits in both tails
  CHECK(med.fold == doctest::Approx(1.0));

  // Both tails land on the same scale, exactly as the empirical rule does.
  record_t low(0, 32, interval_t{1, 20}, false, boost::math::quantile(gamma, 1e-4), 1.0, 0);
  apply_gamma_significance(low, fit);
  CHECK(low.pvalue == doctest::Approx(2e-4).epsilon(0.02));
  record_t high(0, 32, interval_t{1, 20}, false, boost::math::quantile(gamma, 1.0 - 1e-4), 1.0, 0);
  apply_gamma_significance(high, fit);
  CHECK(high.pvalue == doctest::Approx(2e-4).epsilon(0.02));

  // No usable fit, or no distance: the row stays untested.
  const gamma_fit_t no_fit;
  record_t untested(0, 32, interval_t{1, 20}, false, 0.1, 1.0, 0);
  apply_gamma_significance(untested, no_fit);
  CHECK(std::isnan(untested.pvalue));
  CHECK(std::isnan(untested.fold));
  record_t unmapped(0, 32, interval_t{1, 20}, false, nanx(), 1.0, 0);
  apply_gamma_significance(unmapped, fit);
  CHECK(std::isnan(unmapped.pvalue));
}

TEST_CASE("benjamini_hochberg_correction adjusts within each strand and never below the raw p") {
  vec<record_t> recs;
  for (int i = 1; i <= 4; ++i) {
    record_t r(0, 32, interval_t{1, 20}, false, 0.05, 1.0, 0);
    r.pvalue = 0.01 * static_cast<double>(i);
    recs.push_back(r);
  }
  record_t lone(0, 32, interval_t{1, 20}, true, 0.05, 1.0, 0);
  lone.pvalue = 0.001;
  recs.push_back(lone);

  benjamini_hochberg_correction(recs);

  CHECK(recs[4].qvalue == doctest::Approx(recs[4].pvalue)); // a family of one
  for (const record_t& r : recs) {
    CHECK(std::isfinite(r.qvalue));
    CHECK(r.qvalue >= r.pvalue - 1e-12);
    CHECK(r.qvalue <= 1.0);
  }
  CHECK(recs[0].qvalue <= recs[1].qvalue);
  CHECK(recs[1].qvalue <= recs[2].qvalue);
  CHECK(recs[2].qvalue <= recs[3].qvalue);
}

TEST_CASE("apply_empirical_significance skips short pools and unmapped records") {
  const vec<double> short_pool{0.10, 0.11, 0.12};
  record_t r(0, 32, interval_t{1, 20}, false, 0.115, 1.0, 0);
  apply_empirical_significance(r, short_pool);
  CHECK(std::isnan(r.pvalue));
  CHECK(std::isnan(r.fold));

  // A caller that carries a sample size raises the floor; a pool below it stays untested.
  record_t raised(0, 32, interval_t{1, 20}, false, 0.115, 1.0, 0);
  apply_empirical_significance(raised, vec<double>(50, 0.1), 100);
  CHECK(std::isnan(raised.pvalue));

  const vec<double> pool(min_null_samples, 0.1);
  record_t unmapped(0, 32, interval_t{1, 20}, false, nanx(), 1.0, 0);
  apply_empirical_significance(unmapped, pool);
  CHECK(std::isnan(unmapped.pvalue));
}

TEST_CASE("benjamini_hochberg_correction adjusts records on the reported scale") {
  vec<record_t> records;
  records.emplace_back(0, 100, interval_t{1, 50}, false, 0.1, 10.0, 0);
  records.emplace_back(0, 100, interval_t{51, 100}, false, 0.12, 10.0, 0);
  records[0].pvalue = 0.05;
  records[1].pvalue = 0.10;

  benjamini_hochberg_correction(records);

  // BH with m = 2 over {0.05, 0.10}: the largest rank gives 0.10, and the step-up keeps 0.10.
  CHECK(records[0].qvalue == doctest::Approx(0.10));
  CHECK(records[1].qvalue == doctest::Approx(0.10));
}

TEST_CASE("benjamini_hochberg_correction leaves NaN qvalues untouched") {
  vec<record_t> records;
  records.emplace_back(0, 100, interval_t{1, 50}, false, 0.1, 10.0, 0);
  records.emplace_back(0, 100, interval_t{51, 100}, false, 0.12, 10.0, 0);
  records[0].pvalue = 0.05;
  records[1].pvalue = nanx();

  benjamini_hochberg_correction(records);

  // One tested record: the family size is 1, so q equals the reported p.
  CHECK(records[0].qvalue == doctest::Approx(0.05));
  CHECK(std::isnan(records[1].qvalue));
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

TEST_SUITE("bracket_distance") {

TEST_CASE("NaN input gives the full range") {
  const vec<double> th_v{0.1, 0.2, 0.3};
  const auto [lo, hi] = bracket_distance(nanx(), th_v);
  CHECK(lo == doctest::Approx(d_eps));
  CHECK(hi == doctest::Approx(d_ub));
}

TEST_CASE("distance below all thresholds") {
  const vec<double> th_v{0.1, 0.2, 0.3};
  const auto [lo, hi] = bracket_distance(0.05, th_v);
  CHECK(lo == doctest::Approx(d_eps));
  CHECK(hi == doctest::Approx(0.1));
}

TEST_CASE("distance above all thresholds") {
  const vec<double> th_v{0.1, 0.2, 0.3};
  const auto [lo, hi] = bracket_distance(0.5, th_v);
  CHECK(lo == doctest::Approx(0.3));
  CHECK(hi == doctest::Approx(d_ub));
}

TEST_CASE("distance between thresholds") {
  const vec<double> th_v{0.1, 0.2, 0.3};
  const auto [lo, hi] = bracket_distance(0.15, th_v);
  CHECK(lo == doctest::Approx(0.1));
  CHECK(hi == doctest::Approx(0.2));
}

TEST_CASE("distance exactly at a threshold") {
  const vec<double> th_v{0.1, 0.2, 0.3};
  const auto [lo, hi] = bracket_distance(0.2, th_v);
  CHECK(lo == doctest::Approx(0.1));
  CHECK(hi == doctest::Approx(0.2));
}

TEST_CASE("empty thresholds give the full range") {
  const vec<double> th_v{};
  const auto [lo, hi] = bracket_distance(0.15, th_v);
  CHECK(lo == doctest::Approx(d_eps));
  CHECK(hi == doctest::Approx(d_ub));
}

} // TEST_SUITE

TEST_SUITE("empirical significance boundaries") {

TEST_CASE("ties are counted with a <= rank") {
  const vec<double> pool{0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2};
  record_t r(0, 32, interval_t{1, 20}, false, 0.1, 1.0, 0);
  apply_empirical_significance(r, pool);
  CHECK(r.pvalue == doctest::Approx(1.0)); // the median is in both tails
}

TEST_CASE("exactly min_null_samples is enough") {
  vec<double> pool;
  for (size_t i = 0; i < min_null_samples; ++i)
    pool.push_back(0.1 + 0.01 * static_cast<double>(i));
  record_t r(0, 32, interval_t{1, 20}, false, 0.5, 1.0, 0);
  apply_empirical_significance(r, pool);
  CHECK(r.pvalue == doctest::Approx(0.0)); // above every sample
  CHECK(r.fold > 1.0);
}

TEST_CASE("a zero median leaves the fold undefined") {
  const vec<double> pool(min_null_samples, 0.0);
  record_t r(0, 32, interval_t{1, 20}, false, 0.1, 1.0, 0);
  apply_empirical_significance(r, pool);
  CHECK(r.pvalue == doctest::Approx(0.0));
  CHECK(std::isnan(r.fold));
}

} // TEST_SUITE

TEST_SUITE("gamma background thresholds") {

TEST_CASE("four levels resolve into eight sorted gamma thresholds") {
  vec<double> sample;
  for (int i = 0; i < 1000; ++i)
    sample.push_back(0.05 + 0.0001 * static_cast<double>(i));

  gamma_fit_t fit;
  vec<double> th;
  REQUIRE(thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005}, th, min_null_samples, &fit));
  CHECK(fit.shape > 0.0);
  CHECK(fit.scale > 0.0);
  REQUIRE(th.size() == 8);
  CHECK(std::is_sorted(th.begin(), th.end()));
  for (size_t i = 1; i < th.size(); ++i)
    CHECK(th[i - 1] < th[i]); // a smooth fit on a clean sample separates every level
}

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
  REQUIRE(fit_gamma_mom(sample, fit, min_null_samples));
  CHECK(fit.shape == doctest::Approx(shape).epsilon(0.02));
  CHECK(fit.scale == doctest::Approx(scale).epsilon(0.02));
}

TEST_CASE("a pool without variance yields one band per distinct distance") {
  const vec<double> sample(min_null_samples, 0.1);
  gamma_fit_t fit;
  vec<double> th;
  REQUIRE(thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005}, th, min_null_samples, &fit));
  CHECK_FALSE(std::isfinite(fit.shape)); // NaN marks the empirical fallback
  REQUIRE(th.size() == 8);
  CHECK(th.front() == doctest::Approx(th.back())); // one distance, eight thresholds

  // The eight collapse into that one band, and the two bounds join it.
  CHECK(compose_lanes(th, a_range_t{0.05, 0.2, true}) == 3);
  REQUIRE(th.size() == 3);
  CHECK(th[0] == doctest::Approx(0.05)); // the low bound
  CHECK(th[1] == doctest::Approx(0.1));  // the pool's single distance
  CHECK(th[2] == doctest::Approx(0.2));  // the high bound
}

TEST_CASE("a short sample fails") {
  const vec<double> sample{0.1, 0.2, 0.3};
  vec<double> th;
  CHECK_FALSE(thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005}, th, min_null_samples));

  // The caller's minimum - a tenth of --sample-size - gates the gamma fit as well: twelve values
  // fit when ten are asked for, and do not when more are.
  vec<double> twelve;
  for (int i = 0; i < 12; ++i)
    twelve.push_back(0.02 + 0.02 * static_cast<double>(i));
  gamma_fit_t fit;
  CHECK(thresholds_from_levels(twelve, {0.1, 0.05, 0.01, 0.005}, th, 10, &fit));
  CHECK(std::isfinite(fit.shape));
  CHECK_FALSE(thresholds_from_levels(twelve, {0.1, 0.05, 0.01, 0.005}, th, 20, &fit));
}

TEST_CASE("any number of levels resolves into two thresholds each") {
  vec<double> sample;
  for (int i = 0; i < 1000; ++i)
    sample.push_back(0.05 + 0.0001 * static_cast<double>(i));
  vec<double> th;
  REQUIRE(thresholds_from_levels(sample, {0.05}, th, min_null_samples));
  CHECK(th.size() == 2);
  REQUIRE(thresholds_from_levels(sample, {0.05, 0.01}, th, min_null_samples));
  CHECK(th.size() == 4);
  REQUIRE(thresholds_from_levels(sample, {0.02, 0.01, 0.005}, th, min_null_samples));
  CHECK(th.size() == 6);
  REQUIRE(thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005}, th, min_null_samples));
  CHECK(th.size() == 8);
  CHECK_FALSE(thresholds_from_levels(sample, {}, th, min_null_samples));
}

TEST_CASE("the fit carries the floor spike separately and keeps the gamma on the rest") {
  // A pool with a fifth of its windows pinned at the estimable floor: windows with the same
  // histogram get the same MLE, so that is an exact point mass, and a gamma fitted straight to the
  // pool collapses onto it (shape 1.18 here instead of 2). The hurdle counts the spike and fits the
  // windows above it, which recovers both the share and the shape.
  const double shape = 2.0;
  const double scale = 0.165;
  const double floor_d = 1e-9;
  const boost::math::gamma_distribution<double> gamma(shape, scale);
  const size_t n = 500;
  const size_t n_floor = n / 5;

  vec<double> pool(n_floor, floor_d);
  for (size_t i = 0; i < n - n_floor; ++i)
    pool.push_back(boost::math::quantile(gamma, (static_cast<double>(i) + 0.5) / static_cast<double>(n - n_floor)));
  std::sort(pool.begin(), pool.end());

  gamma_fit_t bare;
  REQUIRE(fit_gamma_mom(pool, bare, min_null_samples));
  CHECK(bare.shape < 1.5); // the spike drags the moment fit down

  gamma_fit_t fit;
  REQUIRE(fit_background_gamma(pool, fit, min_null_samples));
  CHECK(fit.p0 == doctest::Approx(0.2).epsilon(0.02));
  CHECK(fit.d_min == doctest::Approx(floor_d));
  CHECK(fit.shape == doctest::Approx(shape).epsilon(0.05));
  CHECK(fit.via_ml);

  // Every level whose tail sits inside the spike reports the floor itself, as the empirical
  // quantiles would, and the rest come from the fitted gamma.
  vec<double> th;
  REQUIRE(thresholds_from_levels(pool, {0.05, 0.01, 0.005}, th, min_null_samples, &fit));
  CHECK(std::count(th.begin(), th.end(), doctest::Approx(floor_d).epsilon(1e-6)) == 3);
  CHECK(th.front() == doctest::Approx(floor_d));
  CHECK(th.back() > 1.0);

  // The hurdle is a distribution over the whole domain, so it also scores records: one at the floor
  // sits in the spike's mass.
  record_t r(0, 100, interval_t{1, 50}, false, floor_d, 1.0, 0);
  apply_gamma_significance(r, fit);
  CHECK(r.pvalue == doctest::Approx(2.0 * fit.p0).epsilon(0.02));
  CHECK(std::isfinite(r.fold));
}

TEST_CASE("the hurdle anchors at the physical floor when the sample minimum is above it")
{
  const double floor_d = 1e-9;
  vec<double> sample;
  for (int i = 0; i < 500; ++i)
    sample.push_back(floor_d + 0.001 + 0.0001 * static_cast<double>(i));
  gamma_fit_t fit;
  REQUIRE(fit_background_gamma(sample, fit, 10, floor_d));
  CHECK(fit.d_min == doctest::Approx(floor_d));
  CHECK(fit.p0 == doctest::Approx(0.0));
  CHECK(fit.shape > 0.0);

  vec<double> th, alpha;
  REQUIRE(level_thresholds(sample, {0.1}, th, alpha, 10, nullptr, floor_d));
  REQUIRE(!th.empty());
  CHECK(th.front() > floor_d);
}

TEST_CASE("the gamma path fits the same pool the empirical path quantiles") {
  // Whatever the --bg-model, both paths read the background null: the fit is not conditioned on
  // lr_ub, because that test is a hard cut at the high bound and would bias every quantile it
  // places (measured: -10 % on the upper threshold and a 31x inflated p-value at shape 91).
  const vec<double> sample(min_null_samples, 0.1);
  vec<double> th;
  gamma_fit_t fit;
  REQUIRE(thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005}, th, min_null_samples, &fit));
  CHECK_FALSE(std::isfinite(fit.shape)); // zero variance: the empirical fallback
  CHECK_FALSE(thresholds_from_levels({}, {0.1}, th, min_null_samples)); // no pool, no thresholds
}

TEST_CASE("compose_lanes keeps thresholds outside the bounds and always adds the bounds") {
  const a_range_t ar{1e-4, 0.5, true};

  // 1e-12 sits below the range and 0.9 above it: both stay lanes, nothing is clamped.
  vec<double> th{1e-12, 1e-4, 5e-3, 0.9};
  CHECK(compose_lanes(th, ar) == 5);
  REQUIRE(th.size() == 5);
  CHECK(th[0] == doctest::Approx(1e-12));
  CHECK(th[1] == doctest::Approx(1e-4)); // a model threshold that equals the bound is one band
  CHECK(th[2] == doctest::Approx(5e-3));
  CHECK(th[3] == doctest::Approx(0.5));
  CHECK(th[4] == doctest::Approx(0.9));

  // A run whose model placed nothing still has the two bounds to report.
  vec<double> none;
  CHECK(compose_lanes(none, ar) == 2);
  REQUIRE(none.size() == 2);
  CHECK(none[0] == doctest::Approx(1e-4));
  CHECK(none[1] == doctest::Approx(0.5));

  // The widest ladder the CLI allows - three levels, two thresholds each - fits the lanes with room
  // to spare for the two bounds.
  vec<double> full;
  for (int i = 0; i < 2 * max_levels; ++i)
    full.push_back(0.01 + 0.01 * static_cast<double>(i));
  CHECK(compose_lanes(full, ar) == 2 * max_levels + 2);
  CHECK(full.size() <= rwidth);
}

TEST_CASE("empirical thresholds are the pool's own quantiles, floor mass included") {
  vec<double> sample(20, 0.0); // a point mass at the floor
  for (int i = 0; i < 200; ++i)
    sample.push_back(0.02 + 0.0005 * static_cast<double>(i));

  vec<double> th;
  REQUIRE(empirical_thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005}, th, min_null_samples));
  REQUIRE(th.size() == 8);
  CHECK(std::is_sorted(th.begin(), th.end()));
  CHECK(th.front() == doctest::Approx(0.0)); // the floor, reported as the pool has it
}

} // TEST_SUITE

TEST_SUITE("hdist threshold")
{
  TEST_CASE("the >20 Mbp cap applies only to an unset --hdist-th")
  {
    const uint64_t small = 19ull * 1000 * 1000;
    const uint64_t large = 20ull * 1000 * 1000;

    // Unset: a large input is capped to 2, a small one keeps the default.
    CHECK(hdist_th_for(large, 0xFFFFFFFFu) == 2);
    CHECK(hdist_th_for(small, 0xFFFFFFFFu) == 3);

    // Explicit: the requested value is honoured whatever the length.
    CHECK(hdist_th_for(large, 3) == 3);
    CHECK(hdist_th_for(large, 5) == 5);

    // An explicit request at or below the cap is unchanged.
    CHECK(hdist_th_for(large, 2) == 2);
    CHECK(hdist_th_for(large, 1) == 1);
  }
}
