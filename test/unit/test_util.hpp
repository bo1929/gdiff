#ifndef _TEST_UTIL_HPP
#define _TEST_UTIL_HPP

// Shared scratch space and fixture generation for the unit tests.
//
// One private scratch directory per test process, removed when the process exits. Tests must
// never write into the repository or into a fixed /tmp name, so concurrent runs cannot collide
// and a finished run leaves nothing behind.

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <unistd.h>

namespace test_util {

struct cleanup_t
{
  ~cleanup_t()
  {
    std::error_code ec;
    std::filesystem::remove_all(dir(), ec);
  }

  static const std::filesystem::path& dir();
};

inline const std::filesystem::path& cleanup_t::dir()
{
  static const std::filesystem::path d = [] {
    auto base = std::filesystem::temp_directory_path() / ("gdiff-test-" + std::to_string(::getpid()));
    std::filesystem::remove_all(base);
    std::filesystem::create_directories(base);
    return base;
  }();
  static const cleanup_t cleanup; // destroyed before `d`, so the path is still valid
  return d;
}

// A path inside the scratch directory; the directory is created on first use.
inline std::filesystem::path path(const std::string& name) { return cleanup_t::dir() / name; }

// The bytes of a file, for comparing two outputs exactly.
inline std::string read_file(const std::filesystem::path& p)
{
  std::ifstream in(p, std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// A deterministic pseudo-random FASTA, wrapped at 80 columns. Tests assert structure and
// geometry, not which k-mers survive, so the composition only has to be reproducible.
inline std::filesystem::path write_fasta(const std::string& name, uint64_t len, uint64_t seed)
{
  const auto p = path(name + ".fa");
  std::mt19937_64 rng(seed);
  const char bases[] = "ACGT";
  std::ofstream out(p);
  out << ">" << name << "\n";
  for (uint64_t i = 0; i < len; ++i) {
    out << bases[rng() & 3];
    if ((i + 1) % 80 == 0) out << "\n";
  }
  out << "\n";
  return p;
}

} // namespace test_util

#endif
