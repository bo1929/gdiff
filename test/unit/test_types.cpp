#include "doctest/doctest.h"
#include "records.hpp"
#include "intext.hpp"
#include "types.hpp"
#include <sstream>

TEST_SUITE("map_params") {

TEST_CASE("bin_size computed from bin_shift") {
  map_params<double> p(0.1, 4, 1000, 0, 33.0);
  CHECK(p.bin_size == 1);      // 2^0 = 1
  CHECK(p.bin_shift == 0);

  map_params<double> p2(0.1, 4, 1000, 3, 33.0);
  CHECK(p2.bin_size == 8);     // 2^3 = 8

  map_params<double> p3(0.1, 4, 1000, 10, 33.0);
  CHECK(p3.bin_size == 1024);  // 2^10 = 1024
}

TEST_CASE("tau_bin is ceil(tau / bin_size)") {
  // tau=1000, bin_size=1 -> tau_bin=1000
  map_params<double> p1(0.1, 4, 1000, 0, 33.0);
  CHECK(p1.tau_bin == 1000);

  // tau=1000, bin_size=8 -> tau_bin=ceil(1000/8)=125
  map_params<double> p2(0.1, 4, 1000, 3, 33.0);
  CHECK(p2.tau_bin == 125);

  // tau=1001, bin_size=8 -> tau_bin=ceil(1001/8)=126
  map_params<double> p3(0.1, 4, 1001, 3, 33.0);
  CHECK(p3.tau_bin == 126);

  // tau=1024, bin_size=1024 -> tau_bin=1
  map_params<double> p4(0.1, 4, 1024, 10, 33.0);
  CHECK(p4.tau_bin == 1);
}

TEST_CASE("map_params with cmlane_t") {
  cmlane_t dths{};
  for (int i = 0; i < 8; ++i) dths[i] = 0.05 * (i + 1);
  map_params<cmlane_t> p(dths, 4, 5000, 2, 33.0);
  CHECK(p.bin_size == 4);
  // tau_bin = ceil(5000/4) = 1250
  CHECK(p.tau_bin == 1250);
}

} // TEST_SUITE

TEST_SUITE("record_t") {

TEST_CASE("reference strand gets two-sided significance percentile") {
  const auto pct_ref = [](double prob, bool on_ref) {
    return on_ref ? (2.0 * std::min(prob, 1.0 - prob)) : prob;
  };
  CHECK(pct_ref(0.01, true) == doctest::Approx(0.02));
  CHECK(pct_ref(0.99, true) == doctest::Approx(0.02));
  CHECK(pct_ref(0.25, false) == doctest::Approx(0.25));
}

TEST_CASE("query-wide MLE selects the reference strand for background sampling") {
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  // Query-level winner used for background samples / add_to_acc: smaller finite strand MLE; fw on tie/NaN.
  const auto winner_rc = [](double fw, double rc) { return strand_diff(fw, rc) > 0.0; };
  CHECK_FALSE(winner_rc(0.2, 0.3)); // fw lower
  CHECK_FALSE(winner_rc(0.3, 0.3)); // tie -> fw
  CHECK(winner_rc(0.3, 0.2));       // rc lower
  CHECK_FALSE(winner_rc(0.2, nan)); // only fw finite -> fw
  CHECK(winner_rc(nan, 0.2));       // only rc finite -> rc
  CHECK_FALSE(winner_rc(nan, nan)); // both NaN -> fw
}

TEST_CASE("map_sequences is_rc matches strand_diff encoding") {
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  const auto pinf = std::numeric_limits<double>::infinity();
  const auto ninf = -std::numeric_limits<double>::infinity();
  const auto is_rc = [](double d_diff) { return (!std::isnan(d_diff)) && d_diff > 0.0; };
  CHECK_FALSE(is_rc(nan));
  CHECK_FALSE(is_rc(ninf));
  CHECK(is_rc(pinf));
  CHECK_FALSE(is_rc(-0.1));
  CHECK_FALSE(is_rc(0.0));
  CHECK(is_rc(0.1));
}

TEST_CASE("strand_diff encodes four strand MLE cases") {
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  const auto pinf = std::numeric_limits<double>::infinity();
  const auto ninf = -std::numeric_limits<double>::infinity();

  CHECK(std::isnan(strand_diff(nan, nan)));
  CHECK(strand_diff(0.2, nan) == ninf);  // only fw finite
  CHECK(strand_diff(nan, 0.2) == pinf);   // only rc finite
  CHECK(strand_diff(0.2, 0.3) == doctest::Approx(-0.1));
}

TEST_CASE("two-sided test and STRAND from d_diff encoding") {
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  const auto pinf = std::numeric_limits<double>::infinity();
  const auto ninf = -std::numeric_limits<double>::infinity();
  const auto two_sided = [](bool rec_is_rc, double d_diff) {
    return !std::isnan(d_diff) && (rec_is_rc == (d_diff > 0.0));
  };

  CHECK(two_sided(false, -0.1));
  CHECK(report_strand(false, -0.1) == '+');
  CHECK_FALSE(two_sided(true, -0.1));
  CHECK(report_strand(true, -0.1) == '-');
  CHECK(two_sided(true, 0.1));
  CHECK(report_strand(true, 0.1) == '+');
  CHECK_FALSE(two_sided(false, 0.1));
  CHECK(report_strand(false, 0.1) == '-');

  CHECK(two_sided(false, 0.0));
  CHECK(report_strand(false, 0.0) == '+');
  CHECK_FALSE(two_sided(true, 0.0));
  CHECK(report_strand(true, 0.0) == '-');

  CHECK_FALSE(two_sided(false, nan));
  CHECK_FALSE(two_sided(true, nan));
  CHECK(report_strand(false, nan) == '.');
  CHECK(report_strand(true, nan) == '.');

  CHECK(two_sided(false, ninf));
  CHECK_FALSE(two_sided(true, ninf));
  CHECK(report_strand(false, ninf) == '+');
  CHECK(report_strand(true, ninf) == '.');

  CHECK_FALSE(two_sided(false, pinf));
  CHECK(two_sided(true, pinf));
  CHECK(report_strand(false, pinf) == '.');
  CHECK(report_strand(true, pinf) == '+');
}

TEST_CASE("canonical NaN d_diff uses one-sided test") {
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  const auto two_sided = [](bool rec_is_rc, double d_diff) {
    return !std::isnan(d_diff) && (rec_is_rc == (d_diff > 0.0));
  };
  CHECK_FALSE(two_sided(false, nan));
  CHECK_FALSE(two_sided(true, nan));
}

TEST_CASE("validate_distance rejects out-of-range and non-finite values") {
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  const auto inf = std::numeric_limits<double>::infinity();
  CHECK(is_valid_distance(0.0));
  CHECK(is_valid_distance(0.5));
  CHECK(is_valid_distance(d_ub - 2.0 * eps));
  CHECK_FALSE(is_valid_distance(d_ub - eps));
  CHECK_FALSE(is_valid_distance(d_ub));
  CHECK_FALSE(is_valid_distance(-0.1));
  CHECK_FALSE(is_valid_distance(nan));
  CHECK_FALSE(is_valid_distance(inf));
  CHECK_FALSE(is_valid_distance(-inf));
  CHECK(std::isnan(validate_distance(d_ub)));
  CHECK(std::isnan(validate_distance(d_ub - eps)));
  CHECK(std::isnan(validate_distance(-0.1)));
  CHECK(std::isnan(validate_distance(nan)));
  CHECK(validate_distance(0.5) == doctest::Approx(0.5));
  CHECK(validate_distance(0.0) == doctest::Approx(0.0));
}

TEST_CASE("coordinate helpers preserve output conventions") {
  const uint64_t enmers = 10;
  const uint32_t k = 5;

  auto iv = get_coordinates({1, 2}, 2, enmers, k);
  CHECK(iv.a == 1);
  CHECK(iv.b == 8); // full bp coverage: (2-1)<<2 + k - 1

  iv = get_coordinates({2, 4}, 2, enmers, k);
  CHECK(iv.a == 5);
  CHECK(iv.b == enmers + k - 1);

  iv = get_coordinates({1, 3}, 0, 2, k);
  CHECK(iv.a == 1);
  CHECK(iv.b == k + 1);
}

} // TEST_SUITE

TEST_SUITE("write_tsv") {

TEST_CASE("basic tab-separated output") {
  std::ostringstream oss;
  write_tsv(oss, "hello", 42, 3.14, "world");
  CHECK(oss.str() == "hello\t42\t3.14\tworld");
}

TEST_CASE("single value") {
  std::ostringstream oss;
  write_tsv(oss, 99);
  CHECK(oss.str() == "99");
}

TEST_CASE("mixed types") {
  std::ostringstream oss;
  write_tsv(oss, "seq1", 1000UL, 200UL, 800UL, "+", "ref.gs", 0.1);
  std::string result = oss.str();
  CHECK(result.find("seq1") == 0);
  // Count tabs
  int tabs = 0;
  for (char c : result) if (c == '\t') tabs++;
  CHECK(tabs == 6);
}

} // TEST_SUITE
