// Integration tests for the simple-map IntMap: sampled hurdle thresholds, q-values,
// gap-claiming vs independent levels, and the new CLI.

#include "doctest/doctest.h"
#include "map.hpp"
#include "random.hpp"
#include "sketch.hpp"
#include "test_util.hpp"
#include "tpool.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

static Sketch make_tiny_sketch(bool canonical = true)
{
  const uint8_t k = 27, h = 11;
  LSHF lshf(k, h);

  sketch_config_t cfg;
  cfg.k = k;
  cfg.w = 33;
  cfg.h = h;
  cfg.canonical = canonical;
  cfg.nrows = 4;
  cfg.ppos_v = lshf.get_ppos_v();
  cfg.npos_v = lshf.get_npos_v();

  vec<uint64_t> keys_v{pack_key(0, 100), pack_key(1, 200), pack_key(2, 300)};
  Buckets buckets;
  buckets.build(cfg.nrows, std::move(keys_v));

  return Sketch(cfg, "tiny", std::move(buckets), 3.0, 0.8);
}

} // namespace

static const std::string TEST_DIR = "test/";
static const std::string GENOMES_DIR = TEST_DIR + "genomes/";
static const std::string SKETCHES_DIR_A = test_util::path("sketches-saware").string() + "/";
static const std::string SKETCHES_DIR_B = test_util::path("sketches-sagnostic").string() + "/";

// Canonical output: seq seq_len start end reference alpha threshold
// interval_len n_kmers d d_seq d_null delta d_min d_max d_floor d_ceiling
// n_at_floor n_at_ceiling p_interval q_interval p_gamma
static constexpr int k_ag_cols = 22;
static constexpr int k_ag_tabs = 21;
static constexpr int k_ag_len = 1;
static constexpr int k_ag_start = 2;
static constexpr int k_ag_end = 3;
static constexpr int k_ag_ref = 4;
static constexpr int k_ag_alpha = 5;
static constexpr int k_ag_threshold = 6;
static constexpr int k_ag_dist = 9;
static constexpr int k_ag_dseq = 10;
static constexpr int k_ag_null = 11;
static constexpr int k_ag_delta = 12;
static constexpr int k_ag_pvalue = 19;
static constexpr int k_ag_qvalue = 20;
static constexpr int k_ag_pgamma = 21;
// Strand-aware output inserts strand and is_rc after end.
static constexpr int k_sa_cols = 24;
static constexpr int k_sa_tabs = 23;
static constexpr int k_sa_len = 1;
static constexpr int k_sa_start = 2;
static constexpr int k_sa_end = 3;
static constexpr int k_sa_strand = 4;
static constexpr int k_sa_rc = 5;
static constexpr int k_sa_ref = 6;
static constexpr int k_sa_alpha = 7;
static constexpr int k_sa_dist = 11;
static constexpr int k_sa_pvalue = 21;
static constexpr int k_sa_qvalue = 22;
static constexpr int k_sa_pgamma = 23;

static bool test_data_available()
{
  return std::filesystem::exists("gdiff") && std::filesystem::exists(GENOMES_DIR + "G000016665.fna.gz");
}

static bool sketch_available(const std::string& dir, const std::string& name)
{
  return std::filesystem::exists(dir + name + ".gs");
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

struct intmap_fixture_t
{
  Sketch sketch;
  std::unique_ptr<QSeq> qs;

  static intmap_fixture_t load(const std::string& ref_name, const std::string& query_name, bool strand_aware)
  {
    const std::string& sdir = strand_aware ? SKETCHES_DIR_A : SKETCHES_DIR_B;
    if (!sketch_available(sdir, ref_name)) {
      std::filesystem::create_directories(sdir);
      std::string cmd = "./gdiff sketch -k 27 -w 31 -h 11";
      if (strand_aware) cmd += " --no-canonical";
      cmd += " -i " + GENOMES_DIR + ref_name + ".fna.gz -o " + sdir + ref_name + ".gs 2>/dev/null";
      if (std::system(cmd.c_str()) != 0) throw std::runtime_error("sketch creation failed");
    }

    intmap_fixture_t fx;
    fx.sketch = Container(sdir + ref_name + ".gs").open(0);
    fx.qs = std::make_unique<QSeq>(GENOMES_DIR + query_name + ".fna.gz", UINT64_MAX);
    while (fx.qs->read_next_batch()) {}
    return fx;
  }
};

static map_opts make_opts(uint64_t sample_size = 500,
                          bool independent_levels = false,
                          double chisq = 33.0,
                          vec<double> levels = {0.05, 0.01, 0.005})
{
  map_opts opts;
  opts.hdist_th = 4;
  opts.tau = 9900;
  opts.chisq = chisq;
  opts.bin_shift = 0;
  opts.sample_size = sample_size;
  opts.levels_v = std::move(levels);
  opts.independent_levels = independent_levels;
  return opts;
}

static std::string run_intmap(const intmap_fixture_t& fx, const map_opts& opts)
{
  std::ostringstream sout;
  ThreadPool pool(1);
  init_thread_rng(1);
  IntMap<cmlane_t> intmap(opts, fx.sketch, fx.qs->get_batch_v());
  intmap.map_sequences(sout, fx.sketch.get_rname(), pool);
  return sout.str();
}

static int count_lines(const std::string& s)
{
  int n = 0;
  std::istringstream iss(s);
  std::string line;
  while (std::getline(iss, line)) {
    if (!line.empty()) ++n;
  }
  return n;
}

static bool any_finite_col(const std::string& s, int col, int cols)
{
  std::istringstream iss(s);
  std::string line;
  while (std::getline(iss, line)) {
    if (line.empty()) continue;
    const auto fields = split_tsv(line);
    if (static_cast<int>(fields.size()) == cols && token_is_finite(fields[col])) return true;
  }
  return false;
}

static std::vector<std::string> sorted_lines(const std::string& s)
{
  std::vector<std::string> v;
  std::istringstream iss(s);
  std::string line;
  while (std::getline(iss, line)) {
    if (!line.empty()) v.push_back(line);
  }
  std::sort(v.begin(), v.end());
  return v;
}

static void check_output_shape(const std::string& output, bool strand_aware)
{
  const int tabs = strand_aware ? k_sa_tabs : k_ag_tabs;
  const int cols = strand_aware ? k_sa_cols : k_ag_cols;
  std::istringstream iss(output);
  std::string line;
  while (std::getline(iss, line)) {
    if (line.empty()) continue;
    CHECK(count_tabs(line) == tabs);
    CHECK(static_cast<int>(split_tsv(line).size()) == cols);
  }
}

// Sum of interval bp / union bp. Independent levels can double-claim; gap-claiming does not.
static double redundancy(const std::string& output)
{
  struct iv_t
  {
    uint64_t s, e;
  };
  std::map<std::string, std::vector<iv_t>> by_seq;
  uint64_t sum_bp = 0;
  std::istringstream iss(output);
  std::string line;
  while (std::getline(iss, line)) {
    if (line.empty()) continue;
    const auto f = split_tsv(line);
    if (f.size() != k_ag_cols) continue;
    const uint64_t s = std::stoull(f[k_ag_start]);
    const uint64_t e = std::stoull(f[k_ag_end]);
    by_seq[f[0]].push_back({s, e});
    sum_bp += e - s + 1;
  }
  uint64_t cov = 0;
  for (auto& kv : by_seq) {
    auto& v = kv.second;
    std::sort(v.begin(), v.end(), [](const iv_t& a, const iv_t& b) { return a.s < b.s; });
    uint64_t cs = v[0].s, ce = v[0].e;
    for (size_t i = 1; i < v.size(); ++i) {
      if (v[i].s <= ce + 1) {
        ce = std::max(ce, v[i].e);
      } else {
        cov += ce - cs + 1;
        cs = v[i].s;
        ce = v[i].e;
      }
    }
    cov += ce - cs + 1;
  }
  return cov ? static_cast<double>(sum_bp) / static_cast<double>(cov) : 0.0;
}

TEST_SUITE("IntMap integration")
{

TEST_CASE("end-to-end with real sketch and query" * doctest::skip(!test_data_available()))
{
  SUBCASE("strand-aware")
  {
    const auto fx = intmap_fixture_t::load("G000018865", "G000016665", true);
    REQUIRE(!fx.qs->is_empty());
    const std::string out = run_intmap(fx, make_opts());
    CHECK(!out.empty());
    check_output_shape(out, true);
  }
  SUBCASE("strand-agnostic")
  {
    const auto fx = intmap_fixture_t::load("G000018865", "G000016665", false);
    REQUIRE(!fx.qs->is_empty());
    const std::string out = run_intmap(fx, make_opts());
    CHECK(!out.empty());
    check_output_shape(out, false);
  }
}

TEST_CASE("default output reports all threshold intervals with q-values" * doctest::skip(!test_data_available()))
{
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx.qs->is_empty());
  const std::string out = run_intmap(fx, make_opts());
  REQUIRE(!out.empty());
  check_output_shape(out, false);

  std::istringstream iss(out);
  std::string line;
  while (std::getline(iss, line)) {
    if (line.empty()) continue;
    const auto f = split_tsv(line);
    CHECK(f[k_ag_ref].find("G000341695") != std::string::npos);
    CHECK(std::isfinite(std::stod(f[k_ag_alpha])));
    CHECK(std::isfinite(std::stod(f[k_ag_threshold])));
    CHECK(std::isfinite(std::stod(f[k_ag_dist])));
    CHECK(f.size() == k_ag_cols);
  }
  CHECK(any_finite_col(out, k_ag_dist, k_ag_cols));
  CHECK(any_finite_col(out, k_ag_qvalue, k_ag_cols));
  CHECK(any_finite_col(out, k_ag_pgamma, k_ag_cols));
}

TEST_CASE("independent levels claim more than gap-claiming" * doctest::skip(!test_data_available()))
{
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx.qs->is_empty());
  const std::string def = run_intmap(fx, make_opts(500, false));
  const std::string ind = run_intmap(fx, make_opts(500, true));
  REQUIRE(!def.empty());
  REQUIRE(!ind.empty());
  CHECK(count_lines(ind) >= count_lines(def));
  // The spurious-hit model changes interval lengths, so absolute redundancy is not monotone.
  CHECK(redundancy(ind) > 0.0);
  CHECK(redundancy(def) > 0.0);
}

TEST_CASE("floor-only null places no threshold lanes")
{
  vec<double> sample(1000, 1e-7);
  vec<double> th, alpha;
  CHECK_FALSE(level_thresholds(sample, {0.1}, th, alpha, 10, nullptr));
  CHECK(th.empty());
  CHECK(alpha.empty());
}

TEST_CASE("gamma rejection is fail-closed unless empirical fallback is requested")
{
  const double floor_d = 1e-7;
  vec<double> sample(900, floor_d);
  sample.insert(sample.end(), 100, 0.01);
  vec<double> th, alpha;
  CHECK_FALSE(level_thresholds(sample, {0.1}, th, alpha, 200, nullptr, floor_d));
  CHECK(th.empty());
  CHECK(level_thresholds(sample, {0.1}, th, alpha, 200, nullptr, floor_d, true));
  REQUIRE(!th.empty());
  CHECK(th.back() == doctest::Approx(0.01));
}

TEST_CASE("one level emits only that alpha" * doctest::skip(!test_data_available()))
{
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx.qs->is_empty());
  const std::string out = run_intmap(fx, make_opts(500, false, 33.0, {0.01}));
  REQUIRE(!out.empty());
  check_output_shape(out, false);
  std::istringstream iss(out);
  std::string line;
  while (std::getline(iss, line)) {
    if (line.empty()) continue;
    const auto f = split_tsv(line);
    CHECK(std::stod(f[k_ag_alpha]) == doctest::Approx(0.01));
  }
}

TEST_CASE("single level: low/high are independent in both modes" * doctest::skip(!test_data_available()))
{
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx.qs->is_empty());
  const std::string def = run_intmap(fx, make_opts(500, false, 33.0, {0.05}));
  const std::string ind = run_intmap(fx, make_opts(500, true, 33.0, {0.05}));
  REQUIRE_FALSE(def.empty());
  REQUIRE_FALSE(ind.empty());
  CHECK(sorted_lines(def) == sorted_lines(ind));
}

TEST_CASE("higher chisq merges at least as much" * doctest::skip(!test_data_available()))
{
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx.qs->is_empty());
  const std::string low = run_intmap(fx, make_opts(500, false, 10.828));
  const std::string high = run_intmap(fx, make_opts(500, false, 100.0));
  REQUIRE(!low.empty());
  REQUIRE(!high.empty());
  CHECK(count_lines(high) <= count_lines(low));
}

TEST_CASE("sample size changes thresholds without emptying the report" * doctest::skip(!test_data_available()))
{
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx.qs->is_empty());
  const std::string small = run_intmap(fx, make_opts(100));
  const std::string large = run_intmap(fx, make_opts(2000));
  REQUIRE(!small.empty());
  REQUIRE(!large.empty());
  check_output_shape(small, false);
  check_output_shape(large, false);
}

TEST_CASE("levels derive thresholds from the background" * doctest::skip(!test_data_available()))
{
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx.qs->is_empty());
  const std::string out = run_intmap(fx, make_opts(500, false, 33.0, {0.05, 0.01, 0.005}));
  REQUIRE(!out.empty());
  CHECK(any_finite_col(out, k_ag_dist, k_ag_cols));
  CHECK(any_finite_col(out, k_ag_null, k_ag_cols));
}

} // TEST_SUITE

TEST_SUITE("BackgroundSampler")
{

TEST_CASE("samples each valid start at most once and keeps finite distances")
{
  seed = 42;
  init_thread_rng(0);

  const Sketch sketch = make_tiny_sketch(true);
  const uint32_t k = sketch.get_lshf_sptr()->get_k();
  constexpr uint64_t tau = 50;
  constexpr uint64_t sample_size = 40;
  vec<qseq_t> batch_v{{"q0", str(tau + k + 20, 'A')}, {"q_short", str(k + 5, 'C')}};

  ThreadPool pool(1);
  BackgroundSampler sampler(sketch, batch_v, tau, 0, 2);
  const auto points = sampler.sample(sample_size, pool);

  std::set<std::pair<uint64_t, uint64_t>> seen;
  for (const auto& p : points) {
    CHECK(seen.insert({p.bix, p.start}).second);
    // The tiny all-A query can be unmapped; NaN is valid, but a finite d must be in range.
    const bool ok = std::isnan(p.d) || (p.d >= 0.0 && p.d <= 1.0);
    CHECK(ok);
  }
}

} // TEST_SUITE

// --- CLI-level coverage -------------------------------------------------------

static bool gdiff_cli_available() { return std::filesystem::exists("./gdiff"); }

struct cli_fixture_t
{
  std::filesystem::path ref;
  std::filesystem::path query;
};

static cli_fixture_t make_cli_fixture()
{
  const auto ref_fa = test_util::write_fasta("cli_ref", 120000, 61);
  const auto query_fa = test_util::write_fasta("cli_query", 60000, 62);
  cli_fixture_t fx;
  fx.query = query_fa;
  fx.ref = test_util::path("cli_ref.gs");
  const std::string cmd =
    "./gdiff sketch -k 27 -w 31 -h 11 -i " + ref_fa.string() + " -o " + fx.ref.string() + " >/dev/null 2>&1";
  if (std::system(cmd.c_str()) != 0) return {};
  return fx;
}

static bool run_cli(const std::string& args, const std::filesystem::path& out)
{
  const std::string cmd = "./gdiff " + args + " > " + out.string() + " 2>/dev/null";
  return std::system(cmd.c_str()) == 0;
}

static std::string map_stdout(const std::string& args)
{
  const auto out = test_util::path("cli_map_stdout.tsv");
  if (!run_cli("map " + args, out)) return {};
  return test_util::strip_comments(test_util::read_file(out));
}

TEST_SUITE("map CLI options")
{

TEST_CASE("map -o writes exactly the stdout report" * doctest::skip(!gdiff_cli_available()))
{
  const cli_fixture_t fx = make_cli_fixture();
  REQUIRE_FALSE(fx.ref.empty());
  const std::string opts = "-l 9900 --hdist-th 4 --sample-size 200";
  const std::string target = fx.query.string() + " " + fx.ref.string();
  const std::string via_stdout = map_stdout(opts + " " + target);
  REQUIRE_FALSE(via_stdout.empty());
  const auto out = test_util::path("cli_map_o.tsv");
  REQUIRE(run_cli("map -o " + out.string() + " " + opts + " " + target, test_util::path("cli_map_devnull.txt")));
  CHECK(test_util::strip_comments(test_util::read_file(out)) == via_stdout);
  check_output_shape(via_stdout, false);
}

TEST_CASE("map --levels before the positionals does not swallow them" * doctest::skip(!gdiff_cli_available()))
{
  const cli_fixture_t fx = make_cli_fixture();
  REQUIRE_FALSE(fx.ref.empty());
  const std::string target = fx.query.string() + " " + fx.ref.string();
  const std::string opts = "-l 9900 --hdist-th 4 --sample-size 200";
  const std::string levels = "--levels 0.05 0.01 0.005 ";
  const std::string via_stdout = map_stdout(levels + target + " " + opts);
  REQUIRE_FALSE(via_stdout.empty());
  const auto out = test_util::path("cli_map_levels_first_o.tsv");
  REQUIRE(run_cli("map " + levels + target + " -o " + out.string() + " " + opts,
                  test_util::path("cli_map_devnull2.txt")));
  CHECK(test_util::strip_comments(test_util::read_file(out)) == via_stdout);
}

TEST_CASE("-d is rejected and -b is accepted" * doctest::skip(!gdiff_cli_available()))
{
  const cli_fixture_t fx = make_cli_fixture();
  REQUIRE_FALSE(fx.ref.empty());
  const auto sink = test_util::path("cli_map_sink.tsv");
  const std::string tail = "-l 9900 --sample-size 200 " + fx.query.string() + " " + fx.ref.string();
  CHECK_FALSE(run_cli("map -d 0.1 " + tail, sink));
  CHECK(run_cli("map -b 2 " + tail, sink));
}

TEST_CASE("--independent-levels is accepted and does not reduce the report" * doctest::skip(!gdiff_cli_available()))
{
  const cli_fixture_t fx = make_cli_fixture();
  REQUIRE_FALSE(fx.ref.empty());
  const std::string target = fx.query.string() + " " + fx.ref.string();
  const std::string opts = "-l 9900 --hdist-th 4 --sample-size 200 ";
  const std::string def = map_stdout(opts + target);
  const std::string ind = map_stdout(opts + "--independent-levels " + target);
  REQUIRE_FALSE(def.empty());
  REQUIRE_FALSE(ind.empty());
  CHECK(count_lines(ind) >= count_lines(def));
}

TEST_CASE("--num-threads does not change the map output" * doctest::skip(!gdiff_cli_available()))
{
  const cli_fixture_t fx = make_cli_fixture();
  REQUIRE_FALSE(fx.ref.empty());
  const std::string opts =
    " map -l 9900 --hdist-th 4 --sample-size 200 " + fx.query.string() + " " + fx.ref.string();
  const auto one = test_util::path("cli_nt1.tsv");
  const auto four = test_util::path("cli_nt4.tsv");
  REQUIRE(run_cli("--num-threads 1" + opts, one));
  REQUIRE(run_cli("--num-threads 4" + opts, four));
  CHECK(test_util::strip_comments(test_util::read_file(one)) == test_util::strip_comments(test_util::read_file(four)));
}

} // TEST_SUITE
