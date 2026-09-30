#include "doctest/doctest.h"
#include "dist.hpp"
#include "map.hpp"
#include <algorithm>
#include <boost/math/distributions/gamma.hpp>
#include <cmath>
#include <random>

TEST_SUITE("empirical significance") {

TEST_CASE("apply_empirical_significance reports the rank and fold of the pool") {
  // Ten samples 0.10 .. 0.19; 0.145 sits at or above exactly five of them.
  vec<double> pool;
  for (uint64_t i = 0; i < 10; ++i)
    pool.push_back(0.10 + 0.01 * static_cast<double>(i));
  const double median = linear_quantile(pool, 0.5);

  record_t r(0, 32, interval_t{1, 20}, false, 0.145, 1.0, 0);
  apply_empirical_significance(r, pool);
  CHECK(r.percentile == doctest::Approx(0.5));
  CHECK(r.fold == doctest::Approx(r.d / median));

  record_t low(0, 32, interval_t{1, 20}, false, 0.05, 1.0, 0);
  apply_empirical_significance(low, pool);
  CHECK(low.percentile == doctest::Approx(0.0));
  CHECK(low.fold == doctest::Approx(0.05 / median));

  // On the reference strand the smaller tail is doubled.
  record_t two_sided(0, 32, interval_t{1, 20}, false, 0.145, 1.0, 0);
  two_sided.d_diff = -0.05;
  apply_empirical_significance(two_sided, pool);
  CHECK(two_sided.percentile == doctest::Approx(1.0));
}

TEST_CASE("apply_empirical_significance skips short pools and unmapped records") {
  const vec<double> short_pool{0.10, 0.11, 0.12};
  record_t r(0, 32, interval_t{1, 20}, false, 0.115, 1.0, 0);
  apply_empirical_significance(r, short_pool);
  CHECK(std::isnan(r.percentile));
  CHECK(std::isnan(r.fold));

  const vec<double> pool(min_null_samples, 0.1);
  record_t unmapped(0, 32, interval_t{1, 20}, false, nanx(), 1.0, 0);
  apply_empirical_significance(unmapped, pool);
  CHECK(std::isnan(unmapped.percentile));
}

TEST_CASE("benjamini_hochberg_correction adjusts canonical records on the two-sided scale") {
  vec<record_t> records;
  records.emplace_back(0, 100, interval_t{1, 50}, false, 0.1, 10.0, 0);
  records.emplace_back(0, 100, interval_t{51, 100}, false, 0.12, 10.0, 0);
  records[0].percentile = 0.05; // canonical records are one-sided (d_diff NaN)
  records[1].percentile = 0.10;

  benjamini_hochberg_correction(records);

  // Two-sided p are 0.10 and 0.20; BH with m = 2 gives 0.20 at both ranks.
  CHECK(records[0].qvalue == doctest::Approx(0.20));
  CHECK(records[1].qvalue == doctest::Approx(0.20));
}

TEST_CASE("benjamini_hochberg_correction leaves NaN qvalues untouched") {
  vec<record_t> records;
  records.emplace_back(0, 100, interval_t{1, 50}, false, 0.1, 10.0, 0);
  records.emplace_back(0, 100, interval_t{51, 100}, false, 0.12, 10.0, 0);
  records[0].percentile = 0.05;
  records[1].percentile = nanx();

  benjamini_hochberg_correction(records);

  // One record: the two-sided p is 0.10 and the family size is 1.
  CHECK(records[0].qvalue == doctest::Approx(0.10));
  CHECK(std::isnan(records[1].qvalue));
}

TEST_CASE("two_sided_p doubles the one-sided percentile and passes two-sided through") {
  record_t q(0, 100, interval_t{1, 50}, false, 0.1, 10.0, 0); // d_diff NaN: one-sided
  q.percentile = 0.05;
  CHECK_FALSE(is_two_sided(q));
  CHECK(two_sided_p(q) == doctest::Approx(0.10));

  q.percentile = 0.90;
  CHECK(two_sided_p(q) == doctest::Approx(0.20)); // 2 * min(0.9, 0.1)

  record_t pref(0, 100, interval_t{1, 50}, false, 0.1, 10.0, 0);
  pref.d_diff = -0.1; // fw is the closer strand, so this fw record is two-sided
  pref.percentile = 0.05;
  CHECK(is_two_sided(pref));
  CHECK(two_sided_p(pref) == doctest::Approx(0.05));
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
  CHECK(r.percentile == doctest::Approx(0.5));
}

TEST_CASE("exactly min_null_samples is enough") {
  vec<double> pool;
  for (size_t i = 0; i < min_null_samples; ++i)
    pool.push_back(0.1 + 0.01 * static_cast<double>(i));
  record_t r(0, 32, interval_t{1, 20}, false, 0.5, 1.0, 0);
  apply_empirical_significance(r, pool);
  CHECK(r.percentile == doctest::Approx(1.0));
  CHECK(r.fold > 1.0);
}

TEST_CASE("a zero median leaves the fold undefined") {
  const vec<double> pool(min_null_samples, 0.0);
  record_t r(0, 32, interval_t{1, 20}, false, 0.1, 1.0, 0);
  apply_empirical_significance(r, pool);
  CHECK(r.percentile == doctest::Approx(1.0));
  CHECK(std::isnan(r.fold));
}

} // TEST_SUITE

TEST_SUITE("gamma background thresholds") {

TEST_CASE("four levels resolve into eight distinct in-range gamma quantiles") {
  vec<double> sample;
  for (int i = 0; i < 1000; ++i)
    sample.push_back(0.05 + 0.0001 * static_cast<double>(i));

  gamma_fit_t fit;
  vec<double> th;
  REQUIRE(thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005}, th, &fit));
  CHECK(fit.shape > 0.0);
  CHECK(fit.scale > 0.0);
  REQUIRE(th.size() == 8);
  for (size_t i = 1; i < th.size(); ++i)
    CHECK(th[i - 1] < th[i]);
  for (const double t : th) {
    CHECK(t > d_eps);
    CHECK(t < d_ub - d_eps);
  }
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
  REQUIRE(fit_gamma_mom(sample, fit));
  CHECK(fit.shape == doctest::Approx(shape).epsilon(0.02));
  CHECK(fit.scale == doctest::Approx(scale).epsilon(0.02));
}

TEST_CASE("a degenerate sample fails instead of duplicating lanes") {
  const vec<double> sample(min_null_samples, 0.1);
  vec<double> th;
  CHECK_FALSE(thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005}, th));
}

TEST_CASE("a short sample fails") {
  const vec<double> sample{0.1, 0.2, 0.3};
  vec<double> th;
  CHECK_FALSE(thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005}, th));
}

TEST_CASE("exactly four levels are required") {
  vec<double> sample;
  for (int i = 0; i < 1000; ++i)
    sample.push_back(0.05 + 0.0001 * static_cast<double>(i));
  vec<double> th;
  CHECK_FALSE(thresholds_from_levels(sample, {0.05, 0.01}, th));
  CHECK_FALSE(thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005, 0.001}, th));
}

TEST_CASE("the lr filter keeps evidence-bearing windows only when enough survive") {
  const auto point = [](double d, double s) {
    sample_point_t p;
    p.d = d;
    p.s = s;
    return p;
  };

  // Three of four windows carry evidence: portion 0.75 > 0.66, so only those are kept.
  vec<sample_point_t> points{point(0.10, 50.0), point(0.20, 40.0), point(0.30, 30.0), point(0.90, 1.0)};
  vec<double> kept = lr_filtered_sample(points, 6.635, 0.66);
  REQUIRE(kept.size() == 3);
  CHECK(kept.back() == doctest::Approx(0.30));

  // Only one of four survives: the fallback keeps every valid window.
  points = {point(0.10, 50.0), point(0.20, 1.0), point(0.30, 1.0), point(0.40, nanx())};
  kept = lr_filtered_sample(points, 6.635, 0.66);
  CHECK(kept.size() == 4);
}

TEST_CASE("a rejected gamma falls back to empirical quantiles") {
  // Twelve distinct values: enough for the empirical quantiles, too few for a two-parameter fit.
  vec<double> sample;
  for (int i = 0; i < 12; ++i)
    sample.push_back(0.02 + 0.02 * static_cast<double>(i));

  gamma_fit_t fit;
  vec<double> th;
  REQUIRE(thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005}, th, &fit));
  CHECK_FALSE(std::isfinite(fit.shape)); // NaN marks the empirical fallback
  REQUIRE(th.size() == 8);
  for (size_t i = 1; i < th.size(); ++i)
    CHECK(th[i - 1] < th[i]);
  for (const double t : th) {
    CHECK(t > d_eps);
    CHECK(t < d_ub - d_eps);
  }
}

TEST_CASE("the empirical fallback lifts floored quantiles to distinct distances") {
  vec<double> sample(20, 0.0); // floored mass, still counted by the quantile
  for (int i = 0; i < 200; ++i)
    sample.push_back(0.02 + 0.0005 * static_cast<double>(i));

  vec<double> th;
  REQUIRE(empirical_thresholds_from_levels(sample, {0.1, 0.05, 0.01, 0.005}, th));
  REQUIRE(th.size() == 8);
  for (size_t i = 1; i < th.size(); ++i)
    CHECK(th[i - 1] < th[i]);
  CHECK(th.front() == doctest::Approx(0.02)); // first distance above the floor
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
