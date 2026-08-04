#include "doctest/doctest.h"
#include "dist.hpp"
#include <cmath>

TEST_SUITE("dist") {

TEST_CASE("cumulative histograms extract exact half-open intervals")
{
  HDHist hist(5, 2, 0);
  hist.aggregate_mer(0, 0);
  hist.aggregate_mer(1, 2);
  hist.aggregate_mer(2, 4);
  hist.compute_prefhistsum();

  vec<uint64_t> counts;
  uint64_t u = 0, t = 0;
  hist.extract_histogram(1, 5, counts, u, t);

  CHECK(counts[0] == 0);
  CHECK(counts[1] == 1);
  CHECK(counts[2] == 1);
  CHECK(t == 2);
  CHECK(u == 0);
}

TEST_CASE("explicit misses accumulate into u")
{
  HDHist hist(5, 2, 0);
  hist.aggregate_mer(0, 0);
  hist.aggregate_mer(3, 1);
  hist.aggregate_mer(5, 2);
  hist.aggregate_mer(1, 3);
  hist.compute_prefhistsum();

  vec<uint64_t> counts;
  uint64_t u = 0, t = 0;
  hist.extract_histogram(0, 5, counts, u, t);

  CHECK(counts[0] == 1);
  CHECK(counts[1] == 1);
  CHECK(t == 2);
  CHECK(u == 2);

  uint64_t u1 = 0, t1 = 0, u2 = 0, t2 = 0;
  vec<uint64_t> c1, c2;
  hist.extract_histogram(0, 2, c1, u1, t1);
  hist.extract_histogram(2, 5, c2, u2, t2);
  CHECK(u1 + u2 == u);
  CHECK(t1 + t2 == t);
}

TEST_CASE("summary uses sample deviation and linear quantiles")
{
  const dist_summary_t summary = summarize_distances({0.0, 1.0, 2.0, 3.0, 4.0, nanx()});

  CHECK(summary.n == 5);
  CHECK(summary.mean == doctest::Approx(2.0));
  CHECK(summary.sd == doctest::Approx(std::sqrt(2.5)));
  CHECK(summary.quantiles[0] == doctest::Approx(0.04));
  CHECK(summary.quantiles[1] == doctest::Approx(0.20));
  CHECK(summary.quantiles[2] == doctest::Approx(1.00));
  CHECK(summary.quantiles[3] == doctest::Approx(2.00));
  CHECK(summary.quantiles[4] == doctest::Approx(3.00));
  CHECK(summary.quantiles[5] == doctest::Approx(3.80));
  CHECK(summary.quantiles[6] == doctest::Approx(3.96));
}

TEST_CASE("sampling is deterministic and preserves exact coordinates")
{
  std::mt19937 rng_a(17);
  std::mt19937 rng_b(17);
  const vec<uint64_t> starts_a = sample_region_starts(76, 200, rng_a);
  const vec<uint64_t> starts_b = sample_region_starts(76, 200, rng_b);

  CHECK(starts_a == starts_b);
  CHECK(starts_a.size() == 200);
  for (const uint64_t start : starts_a) {
    CHECK(start <= 75);
  }
}

TEST_CASE("sampling supports the only possible start")
{
  std::mt19937 rng(2);
  const vec<uint64_t> starts = sample_region_starts(1, 4, rng);
  CHECK(starts == vec<uint64_t>{0, 0, 0, 0});
}

TEST_CASE("strand selection chooses the lower valid distance")
{
  CHECK(select_strand_distance(0.2, 0.1).first == doctest::Approx(0.1));
  CHECK(select_strand_distance(0.2, 0.1).second == '-');
  CHECK(select_strand_distance(0.1, 0.2).second == '+');
  CHECK(select_strand_distance(0.1, nanx()).second == '+');
  CHECK(select_strand_distance(nanx(), 0.1).second == '-');
  CHECK(std::isnan(select_strand_distance(nanx(), nanx()).first));
}

} // TEST_SUITE
