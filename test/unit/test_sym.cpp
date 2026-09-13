// Tie and filter rules of the symmetric merge (port of estimate_sample_ani.c).
#include "doctest/doctest.h"
#include "sym.hpp"
#include <cmath>

namespace {

  // lr_ub well above the default cut; a row is kept unless a test says otherwise.
  constexpr double keep_lr = 100.0;

  dpoint_t row(double d, double lr = keep_lr) { return dpoint_t{d, lr}; }

  vec<double> distances(const sym_merge_t& mg)
  {
    vec<double> out;
    for (const dpoint_t& r : mg.rows) out.push_back(r.d);
    return out;
  }

} // namespace

TEST_SUITE("sym_merge")
{

  TEST_CASE("rank i takes the smaller d of the two directions")
  {
    // Unsorted input: rank i is the i-th closest window, not the i-th sample.
    const sym_merge_t mg = sym_merge({row(0.30), row(0.10)}, {row(0.05), row(0.40)});
    REQUIRE(mg.rows.size() == 2);
    CHECK(mg.rows[0].d == doctest::Approx(0.05)); // min(0.10, 0.05)
    CHECK(mg.rows[1].d == doctest::Approx(0.30)); // min(0.30, 0.40)
    CHECK(mg.n_na == 0);
  }

  TEST_CASE("the merge is symmetric in its arguments")
  {
    const vec<dpoint_t> ab{row(0.2), row(0.4), row(0.1)};
    const vec<dpoint_t> ba{row(0.3), row(0.15), row(0.5)};
    CHECK(distances(sym_merge(ab, ba)) == distances(sym_merge(ba, ab)));
  }

  TEST_CASE("output is ascending by d")
  {
    const sym_merge_t mg =
      sym_merge({row(0.9), row(0.1), row(0.5)}, {row(0.4), row(0.8), row(0.2)});
    REQUIRE(mg.rows.size() == 3);
    for (size_t i = 1; i < mg.rows.size(); ++i)
      CHECK(mg.rows[i - 1].d <= mg.rows[i].d);
  }

  TEST_CASE("a number always beats a NaN at the same rank")
  {
    // NaNs sort last, so rank 1 pairs ab's 0.5 against ba's NaN.
    const sym_merge_t mg = sym_merge({row(0.2), row(0.5)}, {row(0.7), row(nanx())});
    REQUIRE(mg.rows.size() == 2);
    CHECK(mg.rows[0].d == doctest::Approx(0.2)); // min(0.2, 0.7)
    CHECK(mg.rows[1].d == doctest::Approx(0.5)); // 0.5 vs NaN
    CHECK(mg.n_na == 0);
  }

  TEST_CASE("NaNs are sorted to the tail before ranks are paired")
  {
    // After sort, rank 0 is min(0.2, 0.7) and rank 1 is NaN on both sides.
    const sym_merge_t mg = sym_merge({row(0.2), row(nanx())}, {row(nanx()), row(0.7)});
    REQUIRE(mg.rows.size() == 1);
    CHECK(mg.rows[0].d == doctest::Approx(0.2));
    CHECK(mg.n_na == 1);
  }

  TEST_CASE("ranks where both sides are NaN are dropped and counted")
  {
    const sym_merge_t mg =
      sym_merge({row(0.2), row(nanx()), row(nanx())}, {row(0.3), row(nanx()), row(nanx())});
    REQUIRE(mg.rows.size() == 1);
    CHECK(mg.rows[0].d == doctest::Approx(0.2));
    CHECK(mg.n_na == 2);
  }

  TEST_CASE("the row carries its own lr_ub, not the other direction's")
  {
    // The winning side's lr_ub travels with its d.
    const sym_merge_t mg = sym_merge({row(0.10, 7.0)}, {row(0.50, 99.0)});
    REQUIRE(mg.rows.size() == 1);
    CHECK(mg.rows[0].d == doctest::Approx(0.10));
    CHECK(mg.rows[0].lr_ub == doctest::Approx(7.0));
  }

  TEST_CASE("uneven direction sizes keep the longer side's tail")
  {
    const sym_merge_t mg = sym_merge({row(0.1), row(0.2), row(0.3)}, {row(0.25)});
    REQUIRE(mg.rows.size() == 3);
    CHECK(mg.rows[0].d == doctest::Approx(0.1)); // min(0.1, 0.25)
    CHECK(mg.rows[1].d == doctest::Approx(0.2)); // ba exhausted from here
    CHECK(mg.rows[2].d == doctest::Approx(0.3));
    CHECK(mg.n_na == 0);
  }

  TEST_CASE("one empty direction degenerates to the other")
  {
    const sym_merge_t mg = sym_merge({row(0.3), row(0.1)}, {});
    REQUIRE(mg.rows.size() == 2);
    CHECK(mg.rows[0].d == doctest::Approx(0.1));
    CHECK(mg.rows[1].d == doctest::Approx(0.3));
  }

  TEST_CASE("two empty directions yield nothing")
  {
    const sym_merge_t mg = sym_merge({}, {});
    CHECK(mg.rows.empty());
    CHECK(mg.n_na == 0);
  }

} // TEST_SUITE

TEST_SUITE("sym_estimate")
{

  TEST_CASE("with every row significant the estimate is the plain mean")
  {
    const sym_merge_t mg = sym_merge({row(0.1), row(0.2)}, {row(0.1), row(0.2)});
    const sym_est_t est = sym_estimate(mg, lr_th_default, min_portion_default);

    CHECK(est.distance == doctest::Approx(0.15));
    CHECK(est.median == doctest::Approx(0.15));
    CHECK(est.num_filtered == 0);
    CHECK(est.n_total == 2);
    CHECK(est.n_kept == 2);
    CHECK(est.used_filtered);
    CHECK(est.max_unfiltered == doctest::Approx(0.2));
    CHECK(est.max_distance == doctest::Approx(0.2));
    CHECK(est.null_d_v.size() == 2);
  }

  TEST_CASE("rows at or below the lr cut are filtered out of the reported mean")
  {
    // 3 of 4 survive (0.75 > 0.66): filtered mean, 0.9 outlier excluded.
    const vec<dpoint_t> ab{row(0.1), row(0.2), row(0.3), row(0.9, 1.0)};
    const sym_est_t est = sym_estimate(sym_merge(ab, {}), lr_th_default, min_portion_default);

    CHECK(est.used_filtered);
    CHECK(est.distance == doctest::Approx(0.2)); // mean(0.1, 0.2, 0.3)
    CHECK(est.num_filtered == 1);
    CHECK(est.n_total == 4);
    CHECK(est.n_kept == 3);
    CHECK(est.alternative_mean == doctest::Approx(0.375)); // the unfiltered mean
    // max_distance covers kept rows only; max_unfiltered covers all of them.
    CHECK(est.max_distance == doctest::Approx(0.3));
    CHECK(est.max_unfiltered == doctest::Approx(0.9));
  }

  TEST_CASE("the cut is strict: lr_ub exactly at the threshold is rejected")
  {
    const vec<dpoint_t> ab{row(0.1, lr_th_default), row(0.2, lr_th_default)};
    const sym_est_t est = sym_estimate(sym_merge(ab, {}), lr_th_default, min_portion_default);
    CHECK(est.n_kept == 0);
    CHECK(est.num_filtered == 2);
    CHECK_FALSE(est.used_filtered);
    CHECK(est.distance == doctest::Approx(0.15)); // falls back to unfiltered
  }

  TEST_CASE("too few survivors falls back to the unfiltered mean")
  {
    // 1 of 4 survives (0.25 < 0.66): fall back to the unfiltered mean.
    const vec<dpoint_t> ab{row(0.1), row(0.2, 1.0), row(0.3, 1.0), row(0.4, 1.0)};
    const sym_est_t est = sym_estimate(sym_merge(ab, {}), lr_th_default, min_portion_default);

    CHECK_FALSE(est.used_filtered);
    CHECK(est.distance == doctest::Approx(0.25)); // mean of all four
    CHECK(est.alternative_mean == doctest::Approx(0.1));
    CHECK(est.n_kept == 1);
    CHECK(est.null_d_v.size() == 4); // the null follows the reported mean
  }

  TEST_CASE("a missing lr_ub carries no evidence and is rejected")
  {
    const vec<dpoint_t> ab{row(0.1, nanx()), row(0.2, nanx())};
    const sym_est_t est = sym_estimate(sym_merge(ab, {}), lr_th_default, min_portion_default);
    CHECK(est.num_filtered == 2);
    CHECK(est.n_kept == 0);
    CHECK_FALSE(est.used_filtered);
    CHECK(est.distance == doctest::Approx(0.15));
  }

  TEST_CASE("lr_ub of exactly zero is counted separately")
  {
    const vec<dpoint_t> ab{row(0.1, 0.0), row(0.2, 0.0), row(0.3)};
    const sym_est_t est = sym_estimate(sym_merge(ab, {}), lr_th_default, min_portion_default);
    CHECK(est.n_lr_zero == 2);
    CHECK(est.num_filtered == 2);
  }

  TEST_CASE("NA ranks are carried through from the merge")
  {
    const sym_merge_t mg = sym_merge({row(0.2), row(nanx())}, {row(0.3), row(nanx())});
    const sym_est_t est = sym_estimate(mg, lr_th_default, min_portion_default);
    CHECK(est.num_na == 1);
    CHECK(est.n_total == 1);
    CHECK(est.distance == doctest::Approx(0.2));
  }

  TEST_CASE("no rows at all yields NaN rather than zero")
  {
    const sym_est_t est = sym_estimate(sym_merge({}, {}), lr_th_default, min_portion_default);
    CHECK(std::isnan(est.distance));
    CHECK(std::isnan(est.median));
    CHECK(est.n_total == 0);
    CHECK(est.n_kept == 0);
  }

  TEST_CASE("the median is taken over the rows behind the reported estimate")
  {
    // Filtered path: median of the kept rows (0.1, 0.2, 0.3), not of all five.
    const vec<dpoint_t> ab{row(0.1), row(0.2), row(0.3), row(9.0, 1.0)};
    const sym_est_t est = sym_estimate(sym_merge(ab, {}), lr_th_default, 0.5);
    CHECK(est.used_filtered);
    CHECK(est.median == doctest::Approx(0.2));
  }

  TEST_CASE("min_portion of zero always takes the filtered mean when any row survives")
  {
    const vec<dpoint_t> ab{row(0.1), row(0.9, 1.0), row(0.8, 1.0)};
    const sym_est_t est = sym_estimate(sym_merge(ab, {}), lr_th_default, 0.0);
    CHECK(est.used_filtered);
    CHECK(est.distance == doctest::Approx(0.1));
  }

  TEST_CASE("min_portion of one always takes the unfiltered mean")
  {
    const vec<dpoint_t> ab{row(0.1), row(0.2), row(0.3)};
    const sym_est_t est = sym_estimate(sym_merge(ab, {}), lr_th_default, 1.0);
    CHECK_FALSE(est.used_filtered);
    CHECK(est.distance == doctest::Approx(0.2));
  }

} // TEST_SUITE
