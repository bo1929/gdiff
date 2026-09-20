#include "doctest/doctest.h"
#include "map.hpp"
#include <algorithm>
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

TEST_CASE("benjamini_hochberg_correction with canonical-only records") {
  vec<record_t> records;
  records.emplace_back(0, 100, interval_t{1, 50}, false, 0.1, 10.0, 0);
  records.emplace_back(0, 100, interval_t{51, 100}, false, 0.12, 10.0, 0);
  records[0].percentile = 0.05;
  records[1].percentile = 0.10;

  benjamini_hochberg_correction(records);

  CHECK(records[0].qvalue == doctest::Approx(0.10));
  CHECK(records[1].qvalue == doctest::Approx(0.10));
}

TEST_CASE("benjamini_hochberg_correction leaves NaN qvalues untouched") {
  vec<record_t> records;
  records.emplace_back(0, 100, interval_t{1, 50}, false, 0.1, 10.0, 0);
  records.emplace_back(0, 100, interval_t{51, 100}, false, 0.12, 10.0, 0);
  records[0].percentile = 0.05;
  records[1].percentile = nanx();

  benjamini_hochberg_correction(records);

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

TEST_SUITE("thresholds_from_levels") {

TEST_CASE("four levels resolve into eight distinct in-range thresholds") {
  vec<double> pool;
  for (int i = 0; i < 1000; ++i)
    pool.push_back(0.05 + 0.0001 * static_cast<double>(i));

  vec<double> th;
  REQUIRE(thresholds_from_levels(pool, {0.1, 0.05, 0.01, 0.005}, th));
  REQUIRE(th.size() == 8);
  for (size_t i = 1; i < th.size(); ++i)
    CHECK(th[i - 1] < th[i]);
  for (const double t : th) {
    CHECK(t > d_eps);
    CHECK(t < d_ub - d_eps);
  }
}

TEST_CASE("floored lower quantiles are lifted to distinct distances") {
  vec<double> pool(20, 0.0); // floored mass, still counted by the quantile
  for (int i = 0; i < 200; ++i)
    pool.push_back(0.02 + 0.0005 * static_cast<double>(i));

  vec<double> th;
  REQUIRE(thresholds_from_levels(pool, {0.1, 0.05, 0.01, 0.005}, th));
  REQUIRE(th.size() == 8);
  for (size_t i = 1; i < th.size(); ++i)
    CHECK(th[i - 1] < th[i]);
  CHECK(th.front() == doctest::Approx(0.02)); // first distance above the floor
}

TEST_CASE("a coarse pool fails instead of duplicating lanes") {
  const vec<double> pool(min_null_samples, 0.1);
  vec<double> th;
  CHECK_FALSE(thresholds_from_levels(pool, {0.1, 0.05, 0.01, 0.005}, th));
}

TEST_CASE("a short pool fails") {
  const vec<double> pool{0.1, 0.2, 0.3};
  vec<double> th;
  CHECK_FALSE(thresholds_from_levels(pool, {0.1, 0.05, 0.01, 0.005}, th));
}

TEST_CASE("exactly four levels are required") {
  vec<double> pool;
  for (int i = 0; i < 1000; ++i)
    pool.push_back(0.05 + 0.0001 * static_cast<double>(i));
  vec<double> th;
  CHECK_FALSE(thresholds_from_levels(pool, {0.05, 0.01}, th));
  CHECK_FALSE(thresholds_from_levels(pool, {0.1, 0.05, 0.01, 0.005, 0.001}, th));
}

} // TEST_SUITE
