// Integration tests for QIE (query interval estimation) using real test genomes.
#include "doctest/doctest.h"
#include "map.hpp"
#include "sketch.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

static const std::string TEST_DIR = "test/";
static const std::string GENOMES_DIR = TEST_DIR + "genomes/";
static const std::string SKETCHES_DIR = TEST_DIR + "sketches/";

// Unified TSV: 16 columns (15 tabs)
static constexpr int k_out_cols = 16;
static constexpr int k_out_tabs = 15;
static constexpr int k_col_strand = 4;
static constexpr int k_col_is_rc = 5;
static constexpr int k_col_dist = 7;
static constexpr int k_col_mask = 8;
static constexpr int k_col_strand_diff = 12;
static constexpr int k_col_percentile = 13;
static constexpr int k_col_fold = 14;
static constexpr int k_col_qvalue = 15;

static bool test_data_available()
{
  return std::filesystem::exists("gdiff") &&
         std::filesystem::exists(GENOMES_DIR + "G000016665.fna.gz");
}

static bool sketch_available(const std::string& name)
{
  return std::filesystem::exists(SKETCHES_DIR + name + ".skc");
}

static int count_tabs(const std::string& line)
{
  int n = 0;
  for (char c : line) {
    if (c == '\t') ++n;
  }
  return n;
}

static std::vector<std::string> split_tsv(const std::string& line)
{
  std::vector<std::string> fields;
  std::string field;
  for (char c : line) {
    if (c == '\t') {
      fields.push_back(field);
      field.clear();
    } else {
      field += c;
    }
  }
  fields.push_back(field);
  return fields;
}

static bool token_is_finite(const std::string& tok)
{
  if (tok.empty()) return false;
  char* end = nullptr;
  const double v = std::strtod(tok.c_str(), &end);
  return end != tok.c_str() && std::isfinite(v);
}

static bool token_is_nan(const std::string& tok)
{
  if (tok.empty()) return true;
  char* end = nullptr;
  const double v = std::strtod(tok.c_str(), &end);
  return end == tok.c_str() || std::isnan(v);
}

struct qie_fixture_t
{
  sketch_sptr_t sketch;
  qseq_sptr_t qs;

  static qie_fixture_t load(const std::string& ref_name, const std::string& query_name)
  {
    if (!sketch_available(ref_name)) {
      std::filesystem::create_directories(SKETCHES_DIR);
      std::string cmd = "./gdiff sketch -k 27 -w 31 -h 11 -m 2 -r 1 --frac -i " +
                        GENOMES_DIR + ref_name + ".fna.gz -o " +
                        SKETCHES_DIR + ref_name + ".skc 2>/dev/null";
      if (std::system(cmd.c_str()) != 0) {
        throw std::runtime_error("sketch creation failed");
      }
    }

    const std::string sketch_path = SKETCHES_DIR + ref_name + ".skc";
    std::ifstream sketch_stream(sketch_path, std::ifstream::binary);
    uint32_t nsketches = 0;
    sketch_stream.read(reinterpret_cast<char*>(&nsketches), sizeof(uint32_t));

    qie_fixture_t fx;
    fx.sketch = std::make_shared<Sketch>(sketch_path);
    fx.sketch->load_from_offset(sketch_stream, 0);
    sketch_stream.close();
    fx.sketch->make_rho_partial();

    fx.qs = std::make_shared<QSeq>(GENOMES_DIR + query_name + ".fna.gz");
    while (fx.qs->read_next_batch()) {}
    return fx;
  }
};

static std::string run_qie(const qie_fixture_t& fx, const params_t<double>& params)
{
  QIE<double> qie(params, fx.sketch, fx.sketch->get_lshf(), fx.qs->get_seq_batch(), fx.qs->get_qid_batch());
  std::ostringstream sout;
  qie.map_sequences(sout, fx.sketch->get_rid());
  return sout.str();
}

static void check_output_shape(const std::string& output)
{
  std::istringstream iss(output);
  std::string line;
  while (std::getline(iss, line)) {
    if (line.empty()) continue;
    CHECK(count_tabs(line) == k_out_tabs);
    CHECK(static_cast<int>(split_tsv(line).size()) == k_out_cols);
  }
}

TEST_SUITE("QIE integration") {

TEST_CASE("end-to-end with real sketch and query" * doctest::skip(!test_data_available())) {
  const auto fx = qie_fixture_t::load("G000018865", "G000016665");
  REQUIRE(!fx.qs->is_empty());

  SUBCASE("enum lite") {
    params_t<double> params(1, 0.1, 4, 9900, 10000.0, 0, 0, true);
    const std::string output = run_qie(fx, params);
    if (!output.empty()) check_output_shape(output);
  }

  SUBCASE("continuous") {
    params_t<double> params(1, 0.1, 4, 9900, 33.0, 0, 200, false);
    const std::string output = run_qie(fx, params);
    if (!output.empty()) check_output_shape(output);
  }
}

TEST_CASE("known pair: three operating modes" * doctest::skip(!test_data_available())) {
  const auto fx = qie_fixture_t::load("G000341695", "G000025025");
  REQUIRE(!fx.qs->is_empty());

  const params_t<double> p_lite(1, 0.1, 4, 9900, 10000.0, 0, 0, true);
  const params_t<double> p_enum_test(1, 0.1, 4, 9900, 10000.0, 0, 200, true);
  const params_t<double> p_cont(1, 0.1, 4, 9900, 33.0, 0, 200, false);

  const std::string out_lite = run_qie(fx, p_lite);
  const std::string out_enum = run_qie(fx, p_enum_test);
  const std::string out_cont = run_qie(fx, p_cont);

  CHECK(!out_lite.empty());
  CHECK(!out_enum.empty());
  CHECK(!out_cont.empty());

  check_output_shape(out_lite);
  check_output_shape(out_enum);
  check_output_shape(out_cont);

  auto count_lines = [](const std::string& s) {
    int n = 0;
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
      if (!line.empty()) ++n;
    }
    return n;
  };

  auto any_finite_dist = [](const std::string& s) {
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
      if (line.empty()) continue;
      const auto fields = split_tsv(line);
      if (static_cast<int>(fields.size()) != k_out_cols) continue;
      if (!token_is_nan(fields[k_col_dist])) return true;
    }
    return false;
  };

  auto all_dist_nan = [](const std::string& s) {
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
      if (line.empty()) continue;
      const auto fields = split_tsv(line);
      if (static_cast<int>(fields.size()) != k_out_cols) continue;
      if (!token_is_nan(fields[k_col_dist])) return false;
    }
    return true;
  };

  auto has_background_gap = [](const std::string& s) {
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
      if (line.empty()) continue;
      const auto fields = split_tsv(line);
      if (static_cast<int>(fields.size()) != k_out_cols) continue;
      if (fields[k_col_mask] == "0") return true;
    }
    return false;
  };

  CHECK(all_dist_nan(out_lite));
  CHECK_FALSE(any_finite_dist(out_lite));

  CHECK(any_finite_dist(out_enum));
  CHECK(any_finite_dist(out_cont));

  auto any_finite_col = [](const std::string& s, int col) {
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
      if (line.empty()) continue;
      const auto fields = split_tsv(line);
      if (static_cast<int>(fields.size()) != k_out_cols) continue;
      if (token_is_finite(fields[col])) return true;
    }
    return false;
  };

  auto strand_chars_valid = [](const std::string& s) {
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
      if (line.empty()) continue;
      const auto fields = split_tsv(line);
      if (static_cast<int>(fields.size()) != k_out_cols) continue;
      const char c = fields[k_col_strand].empty() ? '.' : fields[k_col_strand][0];
      CHECK((c == '+' || c == '-' || c == '.'));
    }
  };

  strand_chars_valid(out_enum);
  strand_chars_valid(out_cont);

  auto is_rc_valid = [](const std::string& s) {
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
      if (line.empty()) continue;
      const auto fields = split_tsv(line);
      if (static_cast<int>(fields.size()) != k_out_cols) continue;
      CHECK((fields[k_col_is_rc] == "0" || fields[k_col_is_rc] == "1"));
    }
  };

  is_rc_valid(out_enum);
  is_rc_valid(out_cont);
  is_rc_valid(out_lite);

  CHECK(any_finite_col(out_enum, k_col_percentile));
  CHECK(any_finite_col(out_cont, k_col_percentile));
  CHECK(any_finite_col(out_enum, k_col_qvalue));
  CHECK(any_finite_col(out_cont, k_col_qvalue));
  CHECK(any_finite_col(out_enum, k_col_strand_diff));
  CHECK(any_finite_col(out_cont, k_col_strand_diff));

  // Continuous mode only emits interval records (no background gaps).
  CHECK_FALSE(has_background_gap(out_cont));
  CHECK_FALSE(has_background_gap(out_lite));

  const int n_lite = count_lines(out_lite);
  const int n_cont = count_lines(out_cont);
  CHECK(n_cont >= n_lite);
  MESSAGE("lines: enum_lite=", n_lite, " enum+test=", count_lines(out_enum), " continuous=", n_cont);
}

TEST_CASE("QIE with multiple thresholds (cm512_t)" * doctest::skip(!test_data_available())) {
  const auto fx = qie_fixture_t::load("G000341695", "G000025025");
  REQUIRE(!fx.qs->is_empty());

  cm512_t dths{};
  dths[0] = 0.05; dths[1] = 0.10; dths[2] = 0.15; dths[3] = 0.20;
  dths[4] = 0.25; dths[5] = 0.30; dths[6] = 0.35; dths[7] = 0.40;

  params_t<cm512_t> params(8, dths, 4, 9900, 10000.0, 0, 200, true);
  QIE<cm512_t> qie(params, fx.sketch, fx.sketch->get_lshf(), fx.qs->get_seq_batch(), fx.qs->get_qid_batch());

  std::ostringstream sout;
  qie.map_sequences(sout, fx.sketch->get_rid());
  const std::string output = sout.str();

  CHECK(!output.empty());
  check_output_shape(output);
}

} // TEST_SUITE
