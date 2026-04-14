// Integration tests for QIE (query interval estimation) using real test genomes.
#include "doctest/doctest.h"
#include "map.hpp"
#include "sketch.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>

// Path to test data relative to project root
static const std::string TEST_DIR = "test/";
static const std::string GENOMES_DIR = TEST_DIR + "genomes/";
static const std::string SKETCHES_DIR = TEST_DIR + "sketches/";

// Check if the gdiff binary and test data exist
static bool test_data_available()
{
  return std::filesystem::exists("gdiff") &&
         std::filesystem::exists(GENOMES_DIR + "G000016665.fna.gz");
}

// Check if a sketch exists (created by the regression test)
static bool sketch_available(const std::string& name)
{
  return std::filesystem::exists(SKETCHES_DIR + name + ".skc");
}

TEST_SUITE("QIE integration") {

TEST_CASE("end-to-end with real sketch and query" * doctest::skip(!test_data_available())) {
  // This test requires sketches to have been created first (e.g. by running make test-regression)
  // If they don't exist, we create one on the fly using the binary
  const std::string ref_name = "G000018865";
  const std::string query_name = "G000016665";

  if (!sketch_available(ref_name)) {
    // Create sketch directory
    std::filesystem::create_directories(SKETCHES_DIR);
    std::string cmd = "./gdiff sketch -k 27 -w 31 -h 11 -m 2 -r 1 --frac -i " +
                      GENOMES_DIR + ref_name + ".fna.gz -o " +
                      SKETCHES_DIR + ref_name + ".skc 2>/dev/null";
    int ret = std::system(cmd.c_str());
    REQUIRE(ret == 0);
  }

  // Load sketch
  std::string sketch_path = SKETCHES_DIR + ref_name + ".skc";
  REQUIRE(std::filesystem::exists(sketch_path));

  std::ifstream sketch_stream(sketch_path, std::ifstream::binary);
  uint32_t nsketches;
  sketch_stream.read(reinterpret_cast<char*>(&nsketches), sizeof(uint32_t));
  REQUIRE(nsketches >= 1);

  auto sketch = std::make_shared<Sketch>(sketch_path);
  sketch->load_from_offset(sketch_stream, 0);
  sketch_stream.close();
  sketch->make_rho_partial();

  // Load query sequences
  auto qs = std::make_shared<QSeq>(GENOMES_DIR + query_name + ".fna.gz");
  while (qs->read_next_batch()) {}

  REQUIRE(!qs->is_empty());

  SUBCASE("enum_only mode") {
    params_t<double> params(1, 0.1, 4, 9900, 10000.0, 0, 1000, false, true);
    QIE<double> qie(params, sketch, sketch->get_lshf(), qs->get_seq_batch(), qs->get_qid_batch());

    std::ostringstream sout;
    qie.map_sequences(sout, sketch->get_rid());
    std::string output = sout.str();

    // Output should contain tab-separated lines
    // Could be empty (no intervals) for this pair — that's also valid
    // Just verify no crash and well-formed output
    if (!output.empty()) {
      // Each line should have exactly 6 tabs (7 columns)
      std::istringstream iss(output);
      std::string line;
      while (std::getline(iss, line)) {
        if (line.empty()) continue;
        int tabs = 0;
        for (char c : line) if (c == '\t') tabs++;
        CHECK(tabs == 6);
      }
    }
  }

  SUBCASE("continuous mode (non enum_only)") {
    params_t<double> params(1, 0.1, 4, 9900, 33.0, 0, 1000, false, false);
    QIE<double> qie(params, sketch, sketch->get_lshf(), qs->get_seq_batch(), qs->get_qid_batch());

    std::ostringstream sout;
    qie.map_sequences(sout, sketch->get_rid());
    std::string output = sout.str();

    if (!output.empty()) {
      // Continuous mode: 12 tabs (13 columns)
      std::istringstream iss(output);
      std::string line;
      while (std::getline(iss, line)) {
        if (line.empty()) continue;
        int tabs = 0;
        for (char c : line) if (c == '\t') tabs++;
        CHECK(tabs == 12);
      }
    }
  }
}

TEST_CASE("QIE with real data: known pair produces intervals" * doctest::skip(!test_data_available())) {
  // G000025025 vs G000341695 is known to produce intervals (from ground truth)
  const std::string ref_name = "G000341695";
  const std::string query_name = "G000025025";

  if (!sketch_available(ref_name)) {
    std::filesystem::create_directories(SKETCHES_DIR);
    std::string cmd = "./gdiff sketch -k 27 -w 31 -h 11 -m 2 -r 1 --frac -i " +
                      GENOMES_DIR + ref_name + ".fna.gz -o " +
                      SKETCHES_DIR + ref_name + ".skc 2>/dev/null";
    int ret = std::system(cmd.c_str());
    REQUIRE(ret == 0);
  }

  std::string sketch_path = SKETCHES_DIR + ref_name + ".skc";
  std::ifstream sketch_stream(sketch_path, std::ifstream::binary);
  uint32_t nsketches;
  sketch_stream.read(reinterpret_cast<char*>(&nsketches), sizeof(uint32_t));
  auto sketch = std::make_shared<Sketch>(sketch_path);
  sketch->load_from_offset(sketch_stream, 0);
  sketch_stream.close();
  sketch->make_rho_partial();

  auto qs = std::make_shared<QSeq>(GENOMES_DIR + query_name + ".fna.gz");
  while (qs->read_next_batch()) {}
  REQUIRE(!qs->is_empty());

  params_t<double> params(1, 0.1, 4, 9900, 10000.0, 0, 1000, false, true);
  QIE<double> qie(params, sketch, sketch->get_lshf(), qs->get_seq_batch(), qs->get_qid_batch());

  std::ostringstream sout;
  qie.map_sequences(sout, sketch->get_rid());
  std::string output = sout.str();

  // This pair should produce intervals (ground truth has 9 lines = 8 intervals)
  CHECK(!output.empty());

  // Count lines
  int line_count = 0;
  std::istringstream iss(output);
  std::string line;
  while (std::getline(iss, line)) {
    if (!line.empty()) line_count++;
  }
  CHECK(line_count > 0);
  MESSAGE("Found " << line_count << " intervals for " << query_name << " vs " << ref_name);
}

TEST_CASE("QIE with multiple thresholds (cm512_t)" * doctest::skip(!test_data_available())) {
  const std::string ref_name = "G000341695";
  const std::string query_name = "G000025025";

  if (!sketch_available(ref_name)) {
    std::filesystem::create_directories(SKETCHES_DIR);
    std::string cmd = "./gdiff sketch -k 27 -w 31 -h 11 -m 2 -r 1 --frac -i " +
                      GENOMES_DIR + ref_name + ".fna.gz -o " +
                      SKETCHES_DIR + ref_name + ".skc 2>/dev/null";
    int ret = std::system(cmd.c_str());
    REQUIRE(ret == 0);
  }

  std::string sketch_path = SKETCHES_DIR + ref_name + ".skc";
  std::ifstream sketch_stream(sketch_path, std::ifstream::binary);
  uint32_t nsketches;
  sketch_stream.read(reinterpret_cast<char*>(&nsketches), sizeof(uint32_t));
  auto sketch = std::make_shared<Sketch>(sketch_path);
  sketch->load_from_offset(sketch_stream, 0);
  sketch_stream.close();
  sketch->make_rho_partial();

  auto qs = std::make_shared<QSeq>(GENOMES_DIR + query_name + ".fna.gz");
  while (qs->read_next_batch()) {}
  REQUIRE(!qs->is_empty());

  // 8 distance thresholds
  cm512_t dths{};
  dths[0] = 0.05; dths[1] = 0.10; dths[2] = 0.15; dths[3] = 0.20;
  dths[4] = 0.25; dths[5] = 0.30; dths[6] = 0.35; dths[7] = 0.40;

  params_t<cm512_t> params(8, dths, 4, 9900, 10000.0, 0, 1000, false, true);
  QIE<cm512_t> qie(params, sketch, sketch->get_lshf(), qs->get_seq_batch(), qs->get_qid_batch());

  std::ostringstream sout;
  qie.map_sequences(sout, sketch->get_rid());
  std::string output = sout.str();

  // Should produce intervals at least for the smallest threshold
  CHECK(!output.empty());
}

} // TEST_SUITE
