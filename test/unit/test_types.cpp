#include "doctest/doctest.h"
#include "types.hpp"
#include "map.hpp"
#include <sstream>

TEST_SUITE("params_t") {

TEST_CASE("bin_size computed from bin_shift") {
  params_t<double> p(1, 0.1, 4, 1000, 33.0, 0, 1000, false, true);
  CHECK(p.bin_size == 1);      // 2^0 = 1
  CHECK(p.bin_shift == 0);

  params_t<double> p2(1, 0.1, 4, 1000, 33.0, 3, 1000, false, true);
  CHECK(p2.bin_size == 8);     // 2^3 = 8

  params_t<double> p3(1, 0.1, 4, 1000, 33.0, 10, 1000, false, true);
  CHECK(p3.bin_size == 1024);  // 2^10 = 1024
}

TEST_CASE("tau_bin is ceil(tau / bin_size)") {
  // tau=1000, bin_size=1 -> tau_bin=1000
  params_t<double> p1(1, 0.1, 4, 1000, 33.0, 0, 1000, false, true);
  CHECK(p1.tau_bin == 1000);

  // tau=1000, bin_size=8 -> tau_bin=ceil(1000/8)=125
  params_t<double> p2(1, 0.1, 4, 1000, 33.0, 3, 1000, false, true);
  CHECK(p2.tau_bin == 125);

  // tau=1001, bin_size=8 -> tau_bin=ceil(1001/8)=126
  params_t<double> p3(1, 0.1, 4, 1001, 33.0, 3, 1000, false, true);
  CHECK(p3.tau_bin == 126);

  // tau=1024, bin_size=1024 -> tau_bin=1
  params_t<double> p4(1, 0.1, 4, 1024, 33.0, 10, 1000, false, true);
  CHECK(p4.tau_bin == 1);
}

TEST_CASE("params_t with cm512_t") {
  cm512_t dths{};
  for (int i = 0; i < 8; ++i) dths[i] = 0.05 * (i + 1);
  params_t<cm512_t> p(8, dths, 4, 5000, 33.0, 2, 1000, false, false);
  CHECK(p.n == 8);
  CHECK(p.bin_size == 4);
  CHECK(p.enum_only == false);
  CHECK(p.tau == 5000);
  // tau_bin = ceil(5000/4) = 1250
  CHECK(p.tau_bin == 1250);
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
