// roll: window geometry, the -s fallback, and the canonical / strand-aware header shapes.
#include "doctest/doctest.h"
#include "test_util.hpp"
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

  const std::string GDIFF_BIN = "./gdiff";

  bool gdiff_available() { return std::filesystem::exists(GDIFF_BIN); }

  // Reverse complement of the single sequence in `src`.
  std::filesystem::path write_revcomp(const std::string& name, const std::filesystem::path& src)
  {
    std::ifstream in(src);
    std::string line, seq;
    while (std::getline(in, line)) {
      if (!line.empty() && line[0] == '>') continue;
      seq += line;
    }
    std::string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
      switch (*it) {
        case 'A': rc.push_back('T'); break;
        case 'C': rc.push_back('G'); break;
        case 'G': rc.push_back('C'); break;
        case 'T': rc.push_back('A'); break;
        default: rc.push_back('N'); break;
      }
    }
    const auto path = test_util::path(name + ".fa");
    std::ofstream out(path);
    out << ">" << name << "\n";
    for (size_t i = 0; i < rc.size(); i += 80)
      out << rc.substr(i, 80) << "\n";
    return path;
  }

  std::filesystem::path make_sketch(const std::filesystem::path& fa, const std::string& tag, const std::string& extra)
  {
    const auto skc = test_util::path(tag + ".gs");
    const std::string cmd =
      GDIFF_BIN + " sketch -i " + fa.string() + " -o " + skc.string() + " " + extra + " >/dev/null 2>&1";
    if (std::system(cmd.c_str()) != 0) return {};
    return skc;
  }

  // Rolls `query` against `ref` into one output file and returns its non-empty lines.
  std::vector<std::string> roll_lines(const std::filesystem::path& query,
                                     const std::filesystem::path& ref,
                                     const std::string& args)
  {
    const auto out = test_util::path("roll_out.tsv");
    const std::string cmd =
      GDIFF_BIN + " roll " + query.string() + " " + ref.string() + " " + args + " > " + out.string() + " 2>/dev/null";
    if (std::system(cmd.c_str()) != 0) return {};
    std::ifstream in(out);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line))
      if (!line.empty() && line[0] != '#') lines.push_back(line);
    return lines;
  }

  std::vector<std::string> fields(const std::string& line)
  {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string f;
    while (std::getline(ss, f, '\t')) out.push_back(f);
    return out;
  }

  double as_double(const std::string& s) { return std::stod(s); }

} // namespace

TEST_SUITE("roll")
{

  TEST_CASE("canonical: header, window geometry and the -s fallback" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("roll_geom", 4000, 3);
    const auto skc = make_sketch(fa, "roll_geom", "");
    REQUIRE_FALSE(skc.empty());

    // k defaults to 27, so a window of 100 k-mers spans 100 + 27 - 1 = 126 bases.
    const auto stepped = roll_lines(fa, skc, "-l 100 -s 50");
    REQUIRE_FALSE(stepped.empty());
    CHECK(stepped[0] == "seq\tstart\tend\tstrand\treference\td");
    CHECK(stepped.size() == 1 + 78); // (enmers - tau) / step + 1 = (3974 - 100) / 50 + 1

    const auto first = fields(stepped[1]);
    REQUIRE(first.size() == 6);
    CHECK(first[0] == "roll_geom");
    CHECK(first[1] == "1");
    CHECK(first[2] == "126");
    CHECK(first[3] == ".");
    CHECK(first[4] == "roll_geom.fa");
    CHECK(as_double(first[5]) < 1e-6); // self-match

    const auto second = fields(stepped[2]);
    CHECK(second[1] == "51");
    CHECK(second[2] == "176");

    // No -s: the step falls back to -l, so the windows are disjoint.
    const auto disjoint = roll_lines(fa, skc, "-l 100");
    REQUIRE_FALSE(disjoint.empty());
    CHECK(disjoint.size() == 1 + 39); // (3974 - 100) / 100 + 1
    CHECK(fields(disjoint[2])[1] == "101");
  }

  TEST_CASE("strand-aware: both directions, one column each" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("roll_aware", 4000, 4);
    const auto skc = make_sketch(fa, "roll_aware", "--strand-aware");
    REQUIRE_FALSE(skc.empty());

    const auto lines = roll_lines(fa, skc, "-l 100 -s 500");
    REQUIRE_FALSE(lines.empty());
    CHECK(lines[0] == "seq\tstart\tend\treference\td_fw\td_rc");
    CHECK(lines.size() == 1 + 8); // (3974 - 100) / 500 + 1

    // A forward query matches the forward strand only; its reverse complement is the mirror.
    const auto rc = write_revcomp("roll_aware_rc", fa);
    const auto rc_lines = roll_lines(rc, skc, "-l 100 -s 500");
    REQUIRE_FALSE(rc_lines.empty());
    REQUIRE(rc_lines.size() == lines.size());

    const auto fwd = fields(lines[1]);
    const auto rev = fields(rc_lines[1]);
    REQUIRE(fwd.size() == 6);
    REQUIRE(rev.size() == 6);
    CHECK(fwd[0] == "roll_aware");
    CHECK(rev[0] == "roll_aware_rc");
    CHECK(fwd[3] == "roll_aware.fa"); // no strand column here, so reference is field 3

    CHECK(as_double(fwd[4]) < 1e-6);      // forward query -> d_fw
    CHECK_FALSE(as_double(fwd[5]) < 1e-6); //     and not d_rc
    CHECK(as_double(rev[5]) < 1e-6);      // reverse-complement query -> d_rc
    CHECK_FALSE(as_double(rev[4]) < 1e-6); //     and not d_fw
  }

  TEST_CASE("roll -o writes exactly the stdout report" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("roll_o", 4000, 5);
    const auto skc = make_sketch(fa, "roll_o", "");
    REQUIRE_FALSE(skc.empty());

    const auto via_stdout = test_util::path("roll_stdout.tsv");
    const std::string stdout_cmd =
      GDIFF_BIN + " roll -l 100 -s 500 " + fa.string() + " " + skc.string() + " > " + via_stdout.string() + " 2>/dev/null";
    REQUIRE(std::system(stdout_cmd.c_str()) == 0);

    const auto out = test_util::path("roll_o.tsv");
    const std::string file_cmd =
      GDIFF_BIN + " roll -l 100 -s 500 -o " + out.string() + " " + fa.string() + " " + skc.string() + " >/dev/null 2>&1";
    REQUIRE(std::system(file_cmd.c_str()) == 0);
    // The two runs carry different provenance (argv and timestamp), so compare the data.
    CHECK(test_util::strip_comments(test_util::read_file(out)) ==
          test_util::strip_comments(test_util::read_file(via_stdout)));
  }

  TEST_CASE("roll --hdist-th changes matching, not window geometry" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("roll_hd", 4000, 6);
    const auto skc = make_sketch(fa, "roll_hd", "");
    REQUIRE_FALSE(skc.empty());

    const auto strict = roll_lines(fa, skc, "-l 100 -s 500 --hdist-th 0");
    const auto loose = roll_lines(fa, skc, "-l 100 -s 500 --hdist-th 3");
    REQUIRE_FALSE(strict.empty());
    REQUIRE(strict.size() == loose.size());

    for (size_t i = 0; i < strict.size(); ++i) {
      const auto a = fields(strict[i]);
      const auto b = fields(loose[i]);
      REQUIRE(a.size() == 6);
      CHECK(a[1] == b[1]); // start
      CHECK(a[2] == b[2]); // end
    }
  }

} // TEST_SUITE
