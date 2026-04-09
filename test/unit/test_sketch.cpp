#include "doctest/doctest.h"
#include "sketch.hpp"
#include <fstream>
#include <filesystem>
#include <cstring>

// Helper: write a minimal but valid sketch file byte-by-byte
static std::filesystem::path write_tiny_sketch(const std::string& name = "test_sketch")
{
  auto tmp = std::filesystem::temp_directory_path() / (name + ".skc");

  const uint8_t k = 27, w = 33, h = 11;
  const uint32_t m = 2, r = 1;
  const bool frac = true;

  // Build an LSHF just to get valid ppos/npos
  auto lshf_obj = std::make_shared<LSHF>(k, h, m);
  auto ppos = lshf_obj->get_ppos();
  auto npos = lshf_obj->get_npos();

  std::ofstream sout(tmp, std::ofstream::binary);

  uint32_t nsketches = 1;
  sout.write(reinterpret_cast<const char*>(&nsketches), sizeof(uint32_t));

  // Header: rid
  std::string rid = "tiny.skc";
  uint64_t rid_len = rid.size();
  sout.write(reinterpret_cast<const char*>(&rid_len), sizeof(uint64_t));
  sout.write(rid.data(), rid_len);

  // timestamp
  uint64_t timestamp = 1234567890;
  sout.write(reinterpret_cast<const char*>(&timestamp), sizeof(uint64_t));

  // Config
  sout.write(reinterpret_cast<const char*>(&k), sizeof(uint8_t));
  sout.write(reinterpret_cast<const char*>(&w), sizeof(uint8_t));
  sout.write(reinterpret_cast<const char*>(&h), sizeof(uint8_t));
  sout.write(reinterpret_cast<const char*>(&m), sizeof(uint32_t));
  sout.write(reinterpret_cast<const char*>(&r), sizeof(uint32_t));
  sout.write(reinterpret_cast<const char*>(&frac), sizeof(bool));

  uint32_t nrows = 4;
  sout.write(reinterpret_cast<const char*>(&nrows), sizeof(uint32_t));

  // ppos and npos
  sout.write(reinterpret_cast<const char*>(ppos.data()), h * sizeof(uint8_t));
  sout.write(reinterpret_cast<const char*>(npos.data()), (k - h) * sizeof(uint8_t));

  // rho
  double rho = 0.8;
  sout.write(reinterpret_cast<const char*>(&rho), sizeof(double));

  // SFHM data: nkmers, enc_v, nrows, inc_v
  uint64_t nkmers = 3;
  sout.write(reinterpret_cast<const char*>(&nkmers), sizeof(uint64_t));
  enc_t enc_data[3] = {100, 200, 300};
  sout.write(reinterpret_cast<const char*>(enc_data), 3 * sizeof(enc_t));
  sout.write(reinterpret_cast<const char*>(&nrows), sizeof(uint32_t));
  inc_t inc_data[4] = {1, 2, 3, 3};
  sout.write(reinterpret_cast<const char*>(inc_data), 4 * sizeof(inc_t));

  sout.close();
  return tmp;
}

// Helper: load a sketch from a tiny file, skipping the nsketches header
static void load_sketch_from_file(Sketch& sketch, const std::filesystem::path& path)
{
  std::ifstream stream(path, std::ifstream::binary);
  uint32_t ns;
  stream.read(reinterpret_cast<char*>(&ns), sizeof(uint32_t));
  // Use the current stream position (after nsketches) as the offset
  uint64_t offset = static_cast<uint64_t>(stream.tellg());
  sketch.load_from_offset(stream, offset);
  stream.close();
}

TEST_SUITE("Sketch I/O") {

TEST_CASE("load sketch and verify fields") {
  auto path = write_tiny_sketch("test_load");

  Sketch sketch(path);
  load_sketch_from_file(sketch, path);

  CHECK(sketch.get_rid() == "tiny.skc");
  CHECK(sketch.get_timestamp() == 1234567890);
  CHECK(sketch.get_rho() == doctest::Approx(0.8));

  std::filesystem::remove(path);
}

TEST_CASE("partial_offset returns valid values for frac=true, m=2, r=1") {
  auto path = write_tiny_sketch("test_partial");
  Sketch sketch(path);
  load_sketch_from_file(sketch, path);

  // With m=2, r=1, frac=true: rix%m <= r is always true (0<=1, 1<=1)
  uint32_t off0 = sketch.partial_offset(0);
  CHECK(off0 != std::numeric_limits<uint32_t>::max());

  uint32_t off1 = sketch.partial_offset(1);
  CHECK(off1 != std::numeric_limits<uint32_t>::max());

  std::filesystem::remove(path);
}

TEST_CASE("scan_bucket returns false for invalid offset") {
  auto path = write_tiny_sketch("test_scan");
  Sketch sketch(path);
  load_sketch_from_file(sketch, path);

  uint32_t hdist_min;
  bool found = sketch.scan_bucket(std::numeric_limits<uint32_t>::max(), 0, hdist_min);
  CHECK(found == false);

  std::filesystem::remove(path);
}

TEST_CASE("make_rho_partial adjusts rho") {
  auto path = write_tiny_sketch("test_rho");
  Sketch sketch(path);
  load_sketch_from_file(sketch, path);

  double rho_before = sketch.get_rho();
  sketch.make_rho_partial();
  double rho_after = sketch.get_rho();

  // For frac=true, m=2, r=1: rho *= (r+1)/m = 2/2 = 1.0 -> unchanged
  CHECK(rho_after == doctest::Approx(rho_before));

  std::filesystem::remove(path);
}

TEST_CASE("seek_past correctly advances stream position") {
  auto path = write_tiny_sketch("test_seek");

  std::ifstream stream(path, std::ifstream::binary);
  uint32_t ns;
  stream.read(reinterpret_cast<char*>(&ns), sizeof(uint32_t));

  auto pos_before = stream.tellg();
  Sketch::seek_past(stream);
  auto pos_after = stream.tellg();

  // After seeking past one sketch, we should be at the end of file
  CHECK(pos_after > pos_before);
  // Try reading one more byte — should fail (EOF)
  char c;
  stream.read(&c, 1);
  CHECK(stream.eof());

  stream.close();
  std::filesystem::remove(path);
}

} // TEST_SUITE
