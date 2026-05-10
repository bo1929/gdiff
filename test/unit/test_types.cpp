#include "doctest/doctest.h"
#include "types.hpp"
#include "map.hpp"
#include <sstream>

TEST_SUITE("params_t") {

TEST_CASE("bin_size computed from bin_shift") {
  params_t<double> p(1, 0.1, 4, 1000, 33.0, 0, 1000, true);
  CHECK(p.bin_size == 1);      // 2^0 = 1
  CHECK(p.bin_shift == 0);

  params_t<double> p2(1, 0.1, 4, 1000, 33.0, 3, 1000, true);
  CHECK(p2.bin_size == 8);     // 2^3 = 8

  params_t<double> p3(1, 0.1, 4, 1000, 33.0, 10, 1000, true);
  CHECK(p3.bin_size == 1024);  // 2^10 = 1024
}

TEST_CASE("tau_bin is ceil(tau / bin_size)") {
  // tau=1000, bin_size=1 -> tau_bin=1000
  params_t<double> p1(1, 0.1, 4, 1000, 33.0, 0, 1000, true);
  CHECK(p1.tau_bin == 1000);

  // tau=1000, bin_size=8 -> tau_bin=ceil(1000/8)=125
  params_t<double> p2(1, 0.1, 4, 1000, 33.0, 3, 1000, true);
  CHECK(p2.tau_bin == 125);

  // tau=1001, bin_size=8 -> tau_bin=ceil(1001/8)=126
  params_t<double> p3(1, 0.1, 4, 1001, 33.0, 3, 1000, true);
  CHECK(p3.tau_bin == 126);

  // tau=1024, bin_size=1024 -> tau_bin=1
  params_t<double> p4(1, 0.1, 4, 1024, 33.0, 10, 1000, true);
  CHECK(p4.tau_bin == 1);
}

TEST_CASE("params_t with cm512_t") {
  cm512_t dths{};
  for (int i = 0; i < 8; ++i) dths[i] = 0.05 * (i + 1);
  params_t<cm512_t> p(8, dths, 4, 5000, 33.0, 2, 1000, false);
  CHECK(p.n == 8);
  CHECK(p.bin_size == 4);
  CHECK(p.enum_only == false);
  CHECK(p.tau == 5000);
  // tau_bin = ceil(5000/4) = 1250
  CHECK(p.tau_bin == 1250);
}

} // TEST_SUITE

TEST_SUITE("record_t") {

TEST_CASE("intact detection is record-local") {
  record_t full(0, 100, 1, 100, 1, 11, false, true, 0.1, 0.1, 10.0, 1, {0.0, 0.2});
  CHECK(full.is_intact());
  CHECK(full.nbins == 10);
  CHECK(full.get_interval().first == 0);
  CHECK(full.get_interval().second == 10);

  record_t partial(0, 100, 1, 90, 1, 10, true, false, 0.1, 0.1, 10.0, 1, {0.0, 0.2});
  CHECK_FALSE(partial.is_intact());
  CHECK(partial.get_interval().first == 0);
  CHECK(partial.get_interval().second == 9);
}

TEST_CASE("reference strand gets two-sided significance percentile") {
  const auto pct_ref = [](double prob, bool is_ref) {
    return is_ref ? (2.0 * std::min(prob, 1.0 - prob)) : prob;
  };
  CHECK(pct_ref(0.01, true) == doctest::Approx(0.02));
  CHECK(pct_ref(0.99, true) == doctest::Approx(0.02));
  CHECK(pct_ref(0.25, false) == doctest::Approx(0.25));
}

TEST_CASE("null overlap uses half-open bin boundaries") {
  const auto overlaps_half_open = [](const interval_t& lhs, const interval_t& rhs) {
    return lhs.first < rhs.second && rhs.first < lhs.second;
  };
  CHECK(overlaps_half_open({0, 5}, {4, 6}));
  CHECK_FALSE(overlaps_half_open({0, 5}, {5, 8}));
  CHECK_FALSE(overlaps_half_open({5, 8}, {0, 5}));
}

TEST_CASE("coordinate helpers preserve output conventions") {
  const uint64_t enmers = 10;
  const uint32_t k = 5;

  auto row = get_rinterval({1, 2}, 2, enmers, k, false);
  CHECK(row.first == 1);
  CHECK(row.second == 5);

  row = get_rinterval({2, 4}, 2, enmers, k, true);
  CHECK(row.first == 5);
  CHECK(row.second == enmers + k - 1);

  auto iv = get_einterval({1, 3}, 2, enmers, k);
  CHECK(iv.first == 1);
  CHECK(iv.second == 14);

  iv = get_einterval({3, 3}, 2, enmers, k);
  CHECK(iv.first == 9);
  CHECK(iv.second == 14);

  iv = get_einterval({1, 1}, 0, 1, k);
  CHECK(iv.first == 1);
  CHECK(iv.second == k);

  row = get_rinterval({1, 3}, 0, 2, k, true);
  CHECK(row.first == 1);
  CHECK(row.second == k + 1);
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
  write_tsv(oss, "seq1", 1000UL, 200UL, 800UL, "+", "ref.skc", 0.1);
  std::string result = oss.str();
  CHECK(result.find("seq1") == 0);
  // Count tabs
  int tabs = 0;
  for (char c : result) if (c == '\t') tabs++;
  CHECK(tabs == 6);
}

} // TEST_SUITE
