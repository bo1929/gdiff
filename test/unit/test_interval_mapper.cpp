// Integration tests for IntMap (the interval mapper) using real test genomes.
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
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

// Tiny sketch: one encoding in each of the first three buckets of a 4-row table.
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
static const std::string SKETCHES_DIR_A = test_util::path("sketches-saware").string() + "/";  // strand-aware
static const std::string SKETCHES_DIR_B = test_util::path("sketches-sagnostic").string() + "/"; // strand-agnostic

// Strand-aware output: 18 columns (17 tabs); trailing two are info/lr_ub
static constexpr int k_sa_cols = 18;
static constexpr int k_sa_tabs = 17;
static constexpr int k_sa_len = 1;
static constexpr int k_sa_start = 2;
static constexpr int k_sa_end = 3;
static constexpr int k_sa_strand = 4;
static constexpr int k_sa_rc = 5;
static constexpr int k_sa_dist = 7;
static constexpr int k_sa_mask = 8;
static constexpr int k_sa_diff = 11;
static constexpr int k_sa_percentile = 13;
static constexpr int k_sa_qvalue = 15;
static constexpr int k_sa_info = 16;
static constexpr int k_sa_lru = 17;

// Strand-agnostic output: 15 columns (14 tabs); trailing two are info/lr_ub
static constexpr int k_ag_cols = 15;
static constexpr int k_ag_tabs = 14;
static constexpr int k_ag_dist = 5;
static constexpr int k_ag_mask = 6;
static constexpr int k_ag_percentile = 10;
static constexpr int k_ag_qvalue = 12;

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

static bool token_is_nan(const std::string& tok)
{
  if (tok.empty()) return true;
  char* end = nullptr;
  const double v = std::strtod(tok.c_str(), &end);
  return end == tok.c_str() || std::isnan(v);
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
      if (std::system(cmd.c_str()) != 0) {
        throw std::runtime_error("sketch creation failed");
      }
    }

    intmap_fixture_t fx;
    // Sketch keeps the mapping alive; Container is only needed for open().
    fx.sketch = Container(sdir + ref_name + ".gs").open(0);

    // One batch holding the whole query.
    fx.qs = std::make_unique<QSeq>(GENOMES_DIR + query_name + ".fna.gz", UINT64_MAX);
    while (fx.qs->read_next_batch()) {}
    return fx;
  }
};

static map_opts make_opts(double d, uint64_t sample_size, bool enum_only, double chisq = 33.0)
{
  map_opts opts;
  opts.thresholds_v = {d};
  opts.hdist_th = 4;
  opts.tau = 9900;
  opts.chisq = chisq;
  opts.bin_shift = 0;
  opts.sample_size = sample_size;
  opts.enum_only = enum_only;
  return opts;
}

static std::string run_intmap(const intmap_fixture_t& fx, const map_opts& opts)
{
  std::ostringstream sout;
  ThreadPool pool(1);
  const bool scalar = opts.levels_v.empty() && opts.thresholds_v.size() == 1;
  if (scalar) {
    IntMap<double> intmap(opts, fx.sketch, fx.qs->get_batch_v());
    intmap.map_sequences(sout, fx.sketch.get_rname(), pool);
  } else {
    IntMap<cmlane_t> intmap(opts, fx.sketch, fx.qs->get_batch_v());
    intmap.map_sequences(sout, fx.sketch.get_rname(), pool);
  }
  return sout.str();
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
    if (static_cast<int>(fields.size()) != cols) continue;
    if (token_is_finite(fields[col])) return true;
  }
  return false;
}

static bool has_gap_ag(const std::string& s)
{
  std::istringstream iss(s);
  std::string line;
  while (std::getline(iss, line)) {
    if (line.empty()) continue;
    const auto fields = split_tsv(line);
    if (static_cast<int>(fields.size()) != k_ag_cols) continue;
    if (fields[k_ag_mask] == "0") {
      const bool is_full_query = (fields[2] == "1") && (fields[3] == fields[1]);
      if (!is_full_query) return true;
    }
  }
  return false;
}

static void check_ag_cols(const std::string& output, const std::string& ref_name)
{
  std::istringstream iss(output);
  std::string line;
  while (std::getline(iss, line)) {
    if (line.empty()) continue;
    const auto fields = split_tsv(line);
    REQUIRE(static_cast<int>(fields.size()) == k_ag_cols);
    // AG format: no strand / is_rc / d_diff columns; ref id is at index 4.
    CHECK(fields[4].find(ref_name) != std::string::npos);
  }
}

TEST_SUITE("IntMap integration") {

TEST_CASE("end-to-end with real sketch and query" * doctest::skip(!test_data_available())) {
  SUBCASE("strand-aware") {
    const auto fx = intmap_fixture_t::load("G000018865", "G000016665", true);
    REQUIRE(!fx.qs->is_empty());
    const map_opts params = make_opts(0.1, 200, false, 33.0);
    const std::string output = run_intmap(fx, params);
    if (!output.empty()) check_output_shape(output, true);
  }
  SUBCASE("strand-agnostic") {
    const auto fx = intmap_fixture_t::load("G000018865", "G000016665", false);
    REQUIRE(!fx.qs->is_empty());
    const map_opts params = make_opts(0.1, 200, false, 33.0);
    const std::string output = run_intmap(fx, params);
    if (!output.empty()) check_output_shape(output, false);
  }
}

TEST_CASE("known pair: three operating modes, strand-aware" * doctest::skip(!test_data_available())) {
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", true);
  REQUIRE(!fx.qs->is_empty());

  const map_opts p_lite = make_opts(0.1, 0, true, 10000.0);
  const map_opts p_enum_test = make_opts(0.1, 200, true, 10000.0);
  const map_opts p_cont = make_opts(0.1, 200, false, 33.0);

  const std::string out_lite = run_intmap(fx, p_lite);
  const std::string out_enum = run_intmap(fx, p_enum_test);
  const std::string out_cont = run_intmap(fx, p_cont);

  CHECK(!out_lite.empty());
  CHECK(!out_enum.empty());
  CHECK(!out_cont.empty());

  check_output_shape(out_lite, true);
  check_output_shape(out_enum, true);
  check_output_shape(out_cont, true);

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
      if (static_cast<int>(fields.size()) != k_sa_cols) continue;
      if (!token_is_nan(fields[k_sa_dist])) return true;
    }
    return false;
  };

  auto has_background_gap = [](const std::string& s) {
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
      if (line.empty()) continue;
      const auto fields = split_tsv(line);
      if (static_cast<int>(fields.size()) != k_sa_cols) continue;
      if (fields[k_sa_mask] == "0") {
        const bool is_full_query =
          (fields[k_sa_start] == "1") && (fields[k_sa_end] == fields[k_sa_len]);
        if (!is_full_query) return true;
      }
    }
    return false;
  };

  CHECK(any_finite_dist(out_lite));
  CHECK_FALSE(any_finite_col(out_lite, k_sa_percentile, k_sa_cols));

  CHECK(any_finite_dist(out_enum));
  CHECK(any_finite_dist(out_cont));

  auto any_finite_col = [](const std::string& s, int col) {
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
      if (line.empty()) continue;
      const auto fields = split_tsv(line);
      if (static_cast<int>(fields.size()) != k_sa_cols) continue;
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
      if (static_cast<int>(fields.size()) != k_sa_cols) continue;
      const char c = fields[k_sa_strand].empty() ? '.' : fields[k_sa_strand][0];
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
      if (static_cast<int>(fields.size()) != k_sa_cols) continue;
      CHECK((fields[k_sa_rc] == "0" || fields[k_sa_rc] == "1"));
    }
  };

  is_rc_valid(out_enum);
  is_rc_valid(out_cont);
  is_rc_valid(out_lite);

  CHECK(any_finite_col(out_enum, k_sa_percentile));
  CHECK(any_finite_col(out_cont, k_sa_percentile));
  CHECK(any_finite_col(out_enum, k_sa_qvalue));
  CHECK(any_finite_col(out_cont, k_sa_qvalue));
  CHECK(any_finite_col(out_enum, k_sa_diff));
  CHECK(any_finite_col(out_cont, k_sa_diff));
  CHECK(any_finite_col(out_cont, k_sa_info));
  CHECK(any_finite_col(out_cont, k_sa_lru));

  // Continuous mode only emits interval records (no background gaps).
  CHECK_FALSE(has_background_gap(out_cont));
  CHECK_FALSE(has_background_gap(out_lite));

  const int n_lite = count_lines(out_lite);
  const int n_cont = count_lines(out_cont);
  CHECK(n_cont >= n_lite);
  MESSAGE("lines: enum_lite=", n_lite, " enum+test=", count_lines(out_enum), " continuous=", n_cont);
}

TEST_CASE("known pair: three operating modes, strand-agnostic" * doctest::skip(!test_data_available())) {
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx.qs->is_empty());
  REQUIRE(fx.sketch.is_canonical());

  const map_opts p_lite = make_opts(0.1, 0, true, 10000.0);
  const map_opts p_enum_test = make_opts(0.1, 200, true, 10000.0);
  const map_opts p_cont = make_opts(0.1, 200, false, 33.0);

  const std::string out_lite = run_intmap(fx, p_lite);
  const std::string out_enum = run_intmap(fx, p_enum_test);
  const std::string out_cont = run_intmap(fx, p_cont);

  CHECK(!out_lite.empty());
  CHECK(!out_enum.empty());
  CHECK(!out_cont.empty());

  check_output_shape(out_lite, false);
  check_output_shape(out_enum, false);
  check_output_shape(out_cont, false);
  check_ag_cols(out_cont, "G000341695");
  check_ag_cols(out_enum, "G000341695");

  CHECK(any_finite_col(out_lite, k_ag_dist, k_ag_cols));
  CHECK(any_finite_col(out_enum, k_ag_dist, k_ag_cols));
  CHECK(any_finite_col(out_cont, k_ag_dist, k_ag_cols));
  CHECK(any_finite_col(out_enum, k_ag_percentile, k_ag_cols));
  CHECK(any_finite_col(out_cont, k_ag_percentile, k_ag_cols));
  CHECK(any_finite_col(out_enum, k_ag_qvalue, k_ag_cols));
  CHECK(any_finite_col(out_cont, k_ag_qvalue, k_ag_cols));

  CHECK_FALSE(has_gap_ag(out_cont));
  CHECK_FALSE(has_gap_ag(out_lite));

  const int n_lite = count_lines(out_lite);
  const int n_cont = count_lines(out_cont);
  CHECK(n_cont >= n_lite);
}

TEST_CASE("strand-agnostic mode shape and significance" * doctest::skip(!test_data_available())) {
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx.qs->is_empty());
  REQUIRE(fx.sketch.is_canonical());

  const map_opts p_cont = make_opts(0.1, 200, false, 33.0);
  const map_opts p_enum = make_opts(0.1, 200, true, 10000.0);
  const std::string out_cont = run_intmap(fx, p_cont);
  const std::string out_enum = run_intmap(fx, p_enum);

  CHECK(!out_cont.empty());
  CHECK(!out_enum.empty());
  check_output_shape(out_cont, false);
  check_output_shape(out_enum, false);

  auto any_finite_col = [](const std::string& s, int col) {
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
      if (line.empty()) continue;
      const auto fields = split_tsv(line);
      if (static_cast<int>(fields.size()) != k_ag_cols) continue;
      if (token_is_finite(fields[col])) return true;
    }
    return false;
  };

  CHECK(any_finite_col(out_cont, k_ag_percentile));
  CHECK(any_finite_col(out_cont, k_ag_qvalue));
  CHECK(any_finite_col(out_enum, k_ag_percentile));
}

TEST_CASE("SA vs AG same pair both produce output" * doctest::skip(!test_data_available())) {
  auto fx_sa = intmap_fixture_t::load("G000341695", "G000025025", true);
  const auto fx_ag = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx_sa.qs->is_empty());
  REQUIRE(!fx_ag.qs->is_empty());

  const map_opts p = make_opts(0.1, 0, true, 10000.0);
  const std::string out_sa = run_intmap(fx_sa, p);
  const std::string out_ag = run_intmap(fx_ag, p);

  CHECK(!out_sa.empty());
  CHECK(!out_ag.empty());
  CHECK(count_lines(out_sa) > 0);
  CHECK(count_lines(out_ag) > 0);
}

TEST_CASE("IntMap with multiple thresholds (cmlane_t), strand-aware" * doctest::skip(!test_data_available())) {
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", true);
  REQUIRE(!fx.qs->is_empty());

  map_opts opts = make_opts(0.05, 200, true, 10000.0);
  opts.thresholds_v = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40};
  IntMap<cmlane_t> intmap(opts, fx.sketch, fx.qs->get_batch_v());

  std::ostringstream sout;
  ThreadPool pool(1);
  intmap.map_sequences(sout, fx.sketch.get_rname(), pool);
  const std::string output = sout.str();

  CHECK(!output.empty());
  check_output_shape(output, true);
}

TEST_CASE("IntMap with multiple thresholds (cmlane_t), strand-agnostic" * doctest::skip(!test_data_available())) {
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx.qs->is_empty());
  REQUIRE(fx.sketch.is_canonical());

  map_opts opts = make_opts(0.05, 200, true, 10000.0);
  opts.thresholds_v = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40};
  IntMap<cmlane_t> intmap(opts, fx.sketch, fx.qs->get_batch_v());

  std::ostringstream sout;
  ThreadPool pool(1);
  intmap.map_sequences(sout, fx.sketch.get_rname(), pool);
  const std::string output = sout.str();

  CHECK(!output.empty());
  check_output_shape(output, false);
}

TEST_CASE("canonicalize yields a well-formed canonical sketch" * doctest::skip(!test_data_available())) {
  auto fx_sa = intmap_fixture_t::load("G000341695", "G000025025", true);
  const auto fx_ag = intmap_fixture_t::load("G000341695", "G000025025", false);

  REQUIRE_FALSE(fx_sa.sketch.is_canonical());
  CHECK(fx_ag.sketch.is_canonical());

  fx_sa.sketch.canonicalize();
  CHECK(fx_sa.sketch.is_canonical());

  const Buckets& sa = fx_sa.sketch.get_buckets();
  const Buckets& ag = fx_ag.sketch.get_buckets();
  // A forward-strand sketch canonicalized after the fact is not bit-identical to one
  // sketched canonical from the start: minimizer/LSH selection is orientation sensitive,
  // so the retained k-mer sets differ. Both must still be well formed.
  CHECK(sa.get_nrows() == ag.get_nrows());
  CHECK(sa.get_nkmers() > 0);
  CHECK(sa.get_nnonempty() > 0);
  for (uint32_t rix = 0; rix < sa.get_nrows(); ++rix) {
    const enc_t* beg = nullptr;
    const enc_t* end = nullptr;
    if (!sa.range(rix, beg, end)) continue;
    CHECK(end > beg);
  }
}

} // TEST_SUITE

TEST_CASE("levels derive thresholds from the background" * doctest::skip(!test_data_available())) {
  const auto fx = intmap_fixture_t::load("G000341695", "G000025025", false);
  REQUIRE(!fx.qs->is_empty());

  map_opts opts = make_opts(0.10, 500, false, 33.0);
  opts.thresholds_v.clear();
  opts.levels_v = {0.1, 0.05, 0.01, 0.005};

  const std::string out = run_intmap(fx, opts);
  CHECK(!out.empty());
  check_output_shape(out, false);
  CHECK(any_finite_col(out, k_ag_dist, k_ag_cols));
}

TEST_SUITE("BackgroundSampler") {

TEST_CASE("samples each valid start at most once and keeps finite distances") {
  seed = 42;
  init_thread_rng(0);

  const Sketch sketch = make_tiny_sketch(true);
  const uint32_t k = sketch.get_lshf_sptr()->get_k();
  constexpr uint64_t tau = 50;
  constexpr uint64_t sample_size = 40;
  // Long enough for full windows: L >= nwinmers + k - 1 with bin_shift=0.
  vec<qseq_t> batch_v{{"q0", str(tau + k + 20, 'A')}, {"q_short", str(k + 5, 'C')}};

  ThreadPool pool(1);
  BackgroundSampler sampler(sketch, batch_v, tau, 0, 4);
  const vec<sample_point_t> pts = sampler.sample(sample_size, false, pool);

  const uint64_t npos = batch_v[0].seq.size() - k + 1 - tau + 1;
  CHECK(pts.size() == npos);

  uint64_t nfinite = 0;
  std::set<uint64_t> starts;
  for (const sample_point_t& p : pts) {
    CHECK(p.bix == 0); // the short sequence contributes nothing
    if (std::isfinite(p.d)) ++nfinite;
    starts.insert(p.start);
  }
  CHECK(starts.size() == npos);

  vec<double> d_v;
  for (const sample_point_t& p : pts)
    if (is_valid_distance(p.d)) d_v.push_back(p.d);
  CHECK(d_v.size() == nfinite);
  for (const double d : d_v)
    CHECK(std::isfinite(d));
}

TEST_CASE("skips all sequences shorter than the window") {
  seed = 1;
  init_thread_rng(0);

  const Sketch sketch = make_tiny_sketch(true);
  const uint32_t k = sketch.get_lshf_sptr()->get_k();
  vec<qseq_t> batch_v{{"tiny", str(k + 5, 'A')}};

  ThreadPool pool(1);
  BackgroundSampler sampler(sketch, batch_v, 50, 0, 4);
  CHECK(sampler.sample(30, false, pool).empty());
}

TEST_CASE("bin_shift rounds window length up to whole bins") {
  seed = 2;
  init_thread_rng(0);

  const Sketch sketch = make_tiny_sketch(true);
  const uint32_t k = sketch.get_lshf_sptr()->get_k();
  constexpr uint64_t tau = 50;
  constexpr uint64_t bin_shift = 2; // bin_size = 4
  // tau_bin = ceil(50/4) = 13, nwinmers = 52
  constexpr uint64_t expect_nwinmers = 52;
  vec<qseq_t> batch_v{{"q0", str(expect_nwinmers + k + 20, 'A')}};

  ThreadPool pool(1);
  BackgroundSampler sampler(sketch, batch_v, tau, bin_shift, 4);
  const vec<sample_point_t> pts = sampler.sample(16, false, pool);

  const uint64_t enmers = batch_v[0].seq.size() - k + 1;
  const uint64_t npos = (enmers - expect_nwinmers) / (uint64_t(1) << bin_shift) + 1;
  CHECK(pts.size() == npos);
  for (const sample_point_t& p : pts)
    CHECK((p.start << bin_shift) + expect_nwinmers <= enmers);
}

TEST_CASE("per-sequence sampling applies sample_size to each eligible query") {
  seed = 8;
  init_thread_rng(0);

  const Sketch sketch = make_tiny_sketch(true);
  const uint32_t k = sketch.get_lshf_sptr()->get_k();
  constexpr uint64_t tau = 30;
  constexpr uint64_t sample_size = 8;
  vec<qseq_t> batch_v{{"q0", str(tau + k + 20, 'A')}, {"q1", str(tau + k + 30, 'C')}};

  ThreadPool pool(1);
  BackgroundSampler sampler(sketch, batch_v, tau, 0, 4);
  const vec<sample_point_t> pts = sampler.sample(sample_size, true, pool);

  vec<uint64_t> counts(batch_v.size(), 0);
  for (const sample_point_t& p : pts) ++counts[p.bix];
  CHECK(counts[0] == sample_size);
  CHECK(counts[1] == sample_size);
  CHECK(pts.size() == 2 * sample_size);
}

TEST_CASE("non-canonical sketches reconcile both strands into valid distances") {
  seed = 4;
  init_thread_rng(0);

  const Sketch sketch = make_tiny_sketch(false);
  const uint32_t k = sketch.get_lshf_sptr()->get_k();
  constexpr uint64_t tau = 40;
  vec<qseq_t> batch_v{{"q0", str(tau + k + 10, 'A')}};

  ThreadPool pool(1);
  BackgroundSampler sampler(sketch, batch_v, tau, 0, 4);
  const vec<sample_point_t> pts = sampler.sample(12, false, pool);
  REQUIRE(pts.size() == 12);
  // Every window is scored on both strands and reduced to one distance by
  // select_strand_distance; the losing strand must never leak through.
  for (const sample_point_t& p : pts)
    CHECK((std::isnan(p.d) || is_valid_distance(p.d)));
}

} // TEST_SUITE

// --- CLI-level option coverage -----------------------------------------------
// Synthetic inputs, so these run even when the committed genomes are absent.

static bool gdiff_cli_available() { return std::filesystem::exists("./gdiff"); }

struct cli_fixture_t
{
  std::filesystem::path ref;   // reference container
  std::filesystem::path query; // query FASTA
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

// Runs `./gdiff <args>` with stdout captured at `out`; true on a zero exit.
static bool run_cli(const std::string& args, const std::filesystem::path& out)
{
  const std::string cmd = "./gdiff " + args + " > " + out.string() + " 2>/dev/null";
  return std::system(cmd.c_str()) == 0;
}

// The map report `args` produce on stdout, or an empty string when map fails.
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

  const std::string opts = "-d 0.1 -l 9900 --hdist-th 4 --sample-size 200";
  const std::string target = fx.query.string() + " " + fx.ref.string();
  const std::string via_stdout = map_stdout(opts + " " + target);
  REQUIRE_FALSE(via_stdout.empty());

  const auto out = test_util::path("cli_map_o.tsv");
  REQUIRE(run_cli("map -o " + out.string() + " " + opts + " " + target, test_util::path("cli_map_devnull.txt")));
  CHECK(test_util::strip_comments(test_util::read_file(out)) == via_stdout);
  check_output_shape(via_stdout, false);
}

TEST_CASE("map accepts exactly 1 or 8 distance lanes" * doctest::skip(!gdiff_cli_available()))
{
  const cli_fixture_t fx = make_cli_fixture();
  REQUIRE_FALSE(fx.ref.empty());

  const auto sink = test_util::path("cli_map_sink.tsv");
  const std::string tail = "-l 9900 --sample-size 0 " + fx.query.string() + " " + fx.ref.string();

  CHECK(run_cli("map -d 0.1 " + tail, sink));
  CHECK_FALSE(run_cli("map -d 0.1 -d 0.2 -d 0.3 " + tail, sink)); // three lanes is not a valid width
  CHECK(run_cli("map -d 0.05 -d 0.1 -d 0.15 -d 0.2 -d 0.25 -d 0.3 -d 0.35 -d 0.4 " + tail, sink));
}

TEST_CASE("map -b bins the query and --per-sequence takes a per-query null" *
          doctest::skip(!gdiff_cli_available()))
{
  const cli_fixture_t fx = make_cli_fixture();
  REQUIRE_FALSE(fx.ref.empty());

  const std::string tail = "-l 9900 --hdist-th 4 --sample-size 200 " + fx.query.string() + " " + fx.ref.string();

  const std::string binned = map_stdout("-d 0.1 -b 2 " + tail);
  REQUIRE_FALSE(binned.empty());
  check_output_shape(binned, false);

  const std::string per_sequence = map_stdout("-d 0.1 --per-sequence " + tail);
  REQUIRE_FALSE(per_sequence.empty());
  check_output_shape(per_sequence, false);
}

TEST_CASE("--num-threads does not change the map output" * doctest::skip(!gdiff_cli_available()))
{
  const cli_fixture_t fx = make_cli_fixture();
  REQUIRE_FALSE(fx.ref.empty());

  const std::string opts =
    " map -d 0.1 -l 9900 --hdist-th 4 --sample-size 200 " + fx.query.string() + " " + fx.ref.string();
  const auto one = test_util::path("cli_nt1.tsv");
  const auto four = test_util::path("cli_nt4.tsv");
  REQUIRE(run_cli("--num-threads 1" + opts, one));
  REQUIRE(run_cli("--num-threads 4" + opts, four));
  CHECK(test_util::strip_comments(test_util::read_file(one)) == test_util::strip_comments(test_util::read_file(four)));
}

} // TEST_SUITE

