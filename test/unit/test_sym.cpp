// Rank-pairing and filter rules of the symmetric reconciliation.
#include "doctest/doctest.h"
#include "dist.hpp"
#include <cmath>
#include <utility>

namespace {

  // A bound well above the default cut; a window is kept unless a test says otherwise.
  constexpr double keep_lr = 100.0;
  // Mirrors dist's current CLI defaults.
  constexpr double lr_th_default = 10.828;
  constexpr double min_portion_default = 0.66;

  mle_t mle(double d, double s = keep_lr) { return mle_t{d, s}; }

  summary_t summarize(vec<mle_t> ab, vec<mle_t> ba, double lr_th = lr_th_default, double min_portion = min_portion_default)
  {
    return summarize_symmetric(std::move(ab), std::move(ba), lr_th, min_portion);
  }

  // The reconciled distances, ascending. With every window keeping, this is exactly the merged window
  // set, so the rank-pairing rules below are observed through it.
  vec<double> merged(vec<mle_t> ab, vec<mle_t> ba) { return summarize(std::move(ab), std::move(ba)).d_v; }

} // namespace

TEST_SUITE("summarize_symmetric rank pairing")
{

  TEST_CASE("rank i takes the smaller d of the two sides")
  {
    // Unsorted input: rank i is the i-th closest window, not the i-th sample.
    const vec<mle_t> ab{mle(0.30), mle(0.10)};
    const vec<mle_t> ba{mle(0.05), mle(0.40)};
    const vec<double> got = merged(ab, ba);
    REQUIRE(got.size() == 2);
    CHECK(got[0] == doctest::Approx(0.05)); // min(0.10, 0.05)
    CHECK(got[1] == doctest::Approx(0.30)); // min(0.30, 0.40)
    CHECK(summarize(ab, ba).n_na == 0);
  }

  TEST_CASE("rank pairing is symmetric in its arguments")
  {
    const vec<mle_t> ab{mle(0.2), mle(0.4), mle(0.1)};
    const vec<mle_t> ba{mle(0.3), mle(0.15), mle(0.5)};
    CHECK(merged(ab, ba) == merged(ba, ab));
  }

  TEST_CASE("the reported distances ascend")
  {
    const vec<double> got = merged({mle(0.9), mle(0.1), mle(0.5)}, {mle(0.4), mle(0.8), mle(0.2)});
    REQUIRE(got.size() == 3);
    for (size_t i = 1; i < got.size(); ++i)
      CHECK(got[i - 1] <= got[i]);
  }

  TEST_CASE("a number always beats a NaN at the same rank")
  {
    // NaNs sort last, so rank 1 pairs ab's 0.5 against ba's NaN.
    const vec<mle_t> ab{mle(0.2), mle(0.5)};
    const vec<mle_t> ba{mle(0.7), mle(nanx())};
    const vec<double> got = merged(ab, ba);
    REQUIRE(got.size() == 2);
    CHECK(got[0] == doctest::Approx(0.2)); // min(0.2, 0.7)
    CHECK(got[1] == doctest::Approx(0.5)); // 0.5 vs NaN
    CHECK(summarize(ab, ba).n_na == 0);
  }

  TEST_CASE("NaNs are sorted to the tail before ranks are paired")
  {
    // After sorting, rank 0 is min(0.2, 0.7) and rank 1 is NaN on both sides.
    const vec<mle_t> ab{mle(0.2), mle(nanx())};
    const vec<mle_t> ba{mle(nanx()), mle(0.7)};
    const vec<double> got = merged(ab, ba);
    REQUIRE(got.size() == 1);
    CHECK(got[0] == doctest::Approx(0.2));
    CHECK(summarize(ab, ba).n_na == 1);
  }

  TEST_CASE("ranks where both sides are NaN are dropped and counted")
  {
    const vec<mle_t> ab{mle(0.2), mle(nanx()), mle(nanx())};
    const vec<mle_t> ba{mle(0.3), mle(nanx()), mle(nanx())};
    const vec<double> got = merged(ab, ba);
    REQUIRE(got.size() == 1);
    CHECK(got[0] == doctest::Approx(0.2));
    CHECK(summarize(ab, ba).n_na == 2);
  }

  TEST_CASE("the winning window keeps its own bound, not the other side's")
  {
    // ab wins on distance (0.10 < 0.50), so it is ab's bound of 7.0 that the filter sees. With the
    // cut placed between the two bounds, that window is rejected - were ba's 99.0 travelling instead,
    // nothing would be filtered.
    const summary_t s = summarize({mle(0.10, 7.0)}, {mle(0.50, 99.0)}, 10.0);
    CHECK(s.n_filtered == 1);
    CHECK(s.d_v.size() == 1);            // the window still counts, unfiltered
    CHECK(s.d == doctest::Approx(0.10)); // no survivor: the unfiltered mean
  }

  TEST_CASE("uneven side sizes keep the longer side's tail")
  {
    const vec<double> got = merged({mle(0.1), mle(0.2), mle(0.3)}, {mle(0.25)});
    REQUIRE(got.size() == 3);
    CHECK(got[0] == doctest::Approx(0.1)); // min(0.1, 0.25)
    CHECK(got[1] == doctest::Approx(0.2)); // ba exhausted from here
    CHECK(got[2] == doctest::Approx(0.3));
  }

  TEST_CASE("one empty side degenerates to the other")
  {
    const vec<double> got = merged({mle(0.3), mle(0.1)}, {});
    REQUIRE(got.size() == 2);
    CHECK(got[0] == doctest::Approx(0.1));
    CHECK(got[1] == doctest::Approx(0.3));
  }

  TEST_CASE("two empty sides yield nothing")
  {
    CHECK(merged({}, {}).empty());
    CHECK(summarize({}, {}).n_na == 0);
  }

} // TEST_SUITE

TEST_SUITE("summarize_symmetric summary")
{

  TEST_CASE("with every window significant the estimate is the plain mean")
  {
    const summary_t s = summarize({mle(0.1), mle(0.2)}, {mle(0.1), mle(0.2)});
    CHECK(s.d == doctest::Approx(0.15));
    CHECK(s.d_median == doctest::Approx(0.15));
    CHECK(s.n_filtered == 0);
    CHECK(s.d_highest == doctest::Approx(0.2));
    CHECK(s.d_upper == doctest::Approx(0.2));
    CHECK(s.d_v.size() == 2);
  }

  TEST_CASE("windows at or below the lr cut are filtered out of the reported mean")
  {
    // 3 of 4 survive (0.75 > 0.66): filtered mean, 0.9 outlier excluded.
    const summary_t s = summarize({mle(0.1), mle(0.2), mle(0.3), mle(0.9, 1.0)}, {});
    CHECK(s.d == doctest::Approx(0.2)); // mean(0.1, 0.2, 0.3)
    CHECK(s.n_filtered == 1);
    CHECK(s.d_mean == doctest::Approx(0.375)); // the plain mean, which filtering excluded
    // `d_upper` covers kept windows only; `d_highest` covers all of them.
    CHECK(s.d_upper == doctest::Approx(0.3));
    CHECK(s.d_highest == doctest::Approx(0.9));
  }

  TEST_CASE("the cut is strict: a bound exactly at the threshold is rejected")
  {
    const summary_t s = summarize({mle(0.1, lr_th_default), mle(0.2, lr_th_default)}, {});
    CHECK(s.n_filtered == 2);
    CHECK(s.d == doctest::Approx(0.15)); // falls back to unfiltered
  }

  TEST_CASE("too few survivors falls back to the unfiltered mean")
  {
    // 1 of 4 survives (0.25 < 0.66): fall back to the unfiltered mean.
    const summary_t s = summarize({mle(0.1), mle(0.2, 1.0), mle(0.3, 1.0), mle(0.4, 1.0)}, {});
    CHECK(s.d == doctest::Approx(0.25));      // mean of all four: the fallback
    CHECK(s.d_mean == doctest::Approx(0.25)); // the plain mean, always unfiltered
    CHECK(s.d_v.size() == 4);                 // the distances behind the reported mean
  }

  TEST_CASE("a missing bound carries no evidence and is rejected")
  {
    const summary_t s = summarize({mle(0.1, nanx()), mle(0.2, nanx())}, {});
    CHECK(s.n_filtered == 2);
    CHECK(s.d == doctest::Approx(0.15));
  }

  TEST_CASE("a bound of exactly zero is counted separately")
  {
    const summary_t s = summarize({mle(0.1, 0.0), mle(0.2, 0.0), mle(0.3)}, {});
    CHECK(s.n_ub == 2);
    CHECK(s.n_filtered == 2);
  }

  TEST_CASE("NA ranks come through to the estimate")
  {
    const summary_t s = summarize({mle(0.2), mle(nanx())}, {mle(0.3), mle(nanx())});
    CHECK(s.n_na == 1);
    CHECK(s.d_v.size() == 1);
    CHECK(s.d == doctest::Approx(0.2));
  }

  TEST_CASE("no windows at all yields NaN rather than zero")
  {
    const summary_t s = summarize({}, {});
    CHECK(std::isnan(s.d));
    CHECK(std::isnan(s.d_median));
    CHECK(s.d_v.empty());
  }

  TEST_CASE("the median is taken over the windows behind the reported estimate")
  {
    // Filtered path: median of the kept windows (0.1, 0.2, 0.3), not of all four.
    const summary_t s = summarize({mle(0.1), mle(0.2), mle(0.3), mle(9.0, 1.0)}, {}, lr_th_default, 0.5);
    CHECK(s.d_median == doctest::Approx(0.2));
    CHECK(s.d_mean == doctest::Approx(2.4)); // the plain mean, over all four windows
  }

  TEST_CASE("min_portion of zero always takes the filtered mean when any window survives")
  {
    const summary_t s = summarize({mle(0.1), mle(0.9, 1.0), mle(0.8, 1.0)}, {}, lr_th_default, 0.0);
    CHECK(s.d == doctest::Approx(0.1));      // the filtered mean
    CHECK(s.d_mean == doctest::Approx(0.6)); // the plain mean, over all three windows
  }

  TEST_CASE("min_portion of one always takes the unfiltered mean")
  {
    const summary_t s = summarize({mle(0.1), mle(0.2, 1.0), mle(0.3, 1.0)}, {}, lr_th_default, 1.0);
    CHECK(s.d == doctest::Approx(0.2));      // no portion can exceed 1: fallback
    CHECK(s.d_mean == doctest::Approx(0.2)); // the plain mean, which is what it reports
  }

} // TEST_SUITE
