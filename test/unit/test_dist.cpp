#include "doctest/doctest.h"
#include "dim.hpp"
#include "dist.hpp"
#include "gamma.hpp"
#include "random.hpp"
#include "sketch.hpp"
#include "tpool.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>

namespace {

static params_t<double> make_params(uint64_t sample_size = 32)
{
  return params_t<double>(0.1, 4, 2, 33.0, 0, sample_size, true, false);
}

static llh_sptr_t<double> make_map_llhf()
{
  return std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
}

static LLH<double> make_null_llhf()
{
  return LLH<double>(27, 11, 0.5, 4, 0.0, false);
}

static void inject_hits_all_bins(DIM<double>& dim, uint64_t nbins, uint32_t hdist = 0, int reps = 4)
{
  for (uint64_t i = 0; i < nbins; ++i) {
    for (int r = 0; r < reps; ++r) {
      dim.aggregate_mer(hdist, i);
    }
  }
}

static void sample_null(DIM<double>& dim, uint64_t sample_size, vec<sample_t>& out,
                        uint64_t curr_bins = 3, uint64_t tau_bins = 3, uint64_t bix = 0)
{
  dim.sample_intervals(make_null_llhf(), curr_bins, tau_bins, sample_size, bix, out);
}

// Minimal sketch file (same layout as test_sketch.cpp) for DistanceSampler tests.
static sketch_sptr_t load_tiny_sketch(bool canonical = true)
{
  const auto path = std::filesystem::temp_directory_path() / "test_distance_sampler.skc";
  const uint8_t k = 27, w = 33, h = 11;
  auto lshf_obj = std::make_shared<LSHF>(k, h);
  auto ppos = lshf_obj->get_ppos();
  auto npos = lshf_obj->get_npos();

  {
    std::ofstream sout(path, std::ofstream::binary);
    uint32_t nsketches = 1;
    sout.write(reinterpret_cast<const char*>(&nsketches), sizeof(uint32_t));
    const str rid = "tiny.skc";
    uint64_t rid_len = rid.size();
    sout.write(reinterpret_cast<const char*>(&rid_len), sizeof(uint64_t));
    sout.write(rid.data(), rid_len);
    uint64_t timestamp = 1;
    sout.write(reinterpret_cast<const char*>(&timestamp), sizeof(uint64_t));
    sout.write(reinterpret_cast<const char*>(&k), sizeof(uint8_t));
    sout.write(reinterpret_cast<const char*>(&w), sizeof(uint8_t));
    sout.write(reinterpret_cast<const char*>(&h), sizeof(uint8_t));
    sout.write(reinterpret_cast<const char*>(&canonical), sizeof(bool));
    uint32_t nrows = 4;
    sout.write(reinterpret_cast<const char*>(&nrows), sizeof(uint32_t));
    sout.write(reinterpret_cast<const char*>(ppos.data()), h * sizeof(uint8_t));
    sout.write(reinterpret_cast<const char*>(npos.data()), (k - h) * sizeof(uint8_t));
    double rho = 0.8;
    sout.write(reinterpret_cast<const char*>(&rho), sizeof(double));
    uint64_t nkmers = 3;
    sout.write(reinterpret_cast<const char*>(&nkmers), sizeof(uint64_t));
    enc_t enc_data[3] = {100, 200, 300};
    sout.write(reinterpret_cast<const char*>(enc_data), 3 * sizeof(enc_t));
    sout.write(reinterpret_cast<const char*>(&nrows), sizeof(uint32_t));
    inc_t inc_data[4] = {1, 2, 3, 3};
    sout.write(reinterpret_cast<const char*>(inc_data), 4 * sizeof(inc_t));
  }

  auto sketch = std::make_shared<Sketch>(path);
  std::ifstream stream(path, std::ifstream::binary);
  uint32_t ns = 0;
  stream.read(reinterpret_cast<char*>(&ns), sizeof(uint32_t));
  sketch->load_from_offset(stream, static_cast<uint64_t>(stream.tellg()));
  stream.close();
  std::filesystem::remove(path);
  return sketch;
}

} // namespace

TEST_SUITE("null sampling and significance") {

TEST_CASE("sample_intervals collects windows from hit-rich query") {
  seed = 7;
  init_thread_rng(0);

  auto params = make_params(20);
  auto llhf = make_map_llhf();
  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();

  vec<sample_t> null_samples;
  sample_null(dim, params.sample_size, null_samples);
  REQUIRE(null_samples.size() >= GammaModel::min_nsamples);
  for (const auto& s : null_samples) {
    CHECK(std::isfinite(s.d));
    CHECK(std::isfinite(s.I));
    CHECK(s.bix == 0);
    CHECK(s.bin_iv.a < s.bin_iv.b);
  }
}

TEST_CASE("sample_intervals yields empty pool when query has no hits") {
  seed = 3;
  init_thread_rng(0);

  auto params = make_params(20);
  DIM<double> dim(params, make_map_llhf(), 64, 64);
  dim.compute_prefhistsum();

  vec<sample_t> null_samples;
  sample_null(dim, params.sample_size, null_samples);
  CHECK(null_samples.empty());
}

TEST_CASE("sample_intervals falls back to tau_bins when curr_bins cannot fill") {
  seed = 5;
  init_thread_rng(0);

  auto params = make_params(20);
  auto llhf = make_map_llhf();
  DIM<double> dim(params, llhf, 10, 10);
  inject_hits_all_bins(dim, 10);
  dim.compute_prefhistsum();

  vec<sample_t> null_samples;
  sample_null(dim, params.sample_size, null_samples, 16, 3, 0);
  REQUIRE(null_samples.size() >= GammaModel::min_nsamples);
  for (const auto& s : null_samples)
    CHECK(s.bin_iv.b - s.bin_iv.a == 3);
}

TEST_CASE("sample_intervals skips when window exceeds query mers") {
  auto params = make_params(20);
  auto llhf = make_map_llhf();
  DIM<double> dim(params, llhf, 2, 2);
  inject_hits_all_bins(dim, 2, 0, 2);
  dim.compute_prefhistsum();

  vec<sample_t> null_samples;
  sample_null(dim, params.sample_size, null_samples, 3, 3, 0);
  CHECK(null_samples.empty());
}

TEST_CASE("test_significance scores a single record") {
  seed = 11;
  init_thread_rng(0);

  auto params = make_params(20);
  auto llhf = make_map_llhf();
  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();

  vec<sample_t> null_samples;
  sample_null(dim, params.sample_size, null_samples);
  REQUIRE(null_samples.size() >= GammaModel::min_nsamples);

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(4, 8, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(v.data(), u, d_obs);

  record_t r(0, 32, interval_t{1, 20}, interval_t{5, 9}, false, d_obs, I_obs, 0);
  r.d_diff = -0.05;
  r.d_q = 0.2;

  CHECK(test_significance(r, null_samples, params.sample_size, "q0"));
  CHECK(std::isfinite(r.percentile));
  CHECK(r.percentile >= 0.0);
  CHECK(r.percentile <= 1.0);

  const bool two_sided = !std::isnan(r.d_diff) && (r.is_rc == (r.d_diff > 0.0));
  CHECK(two_sided);
}

TEST_CASE("benjamini_hochberg_correction assigns qvalues per strand") {
  seed = 11;
  init_thread_rng(0);

  auto params = make_params(20);
  auto llhf = make_map_llhf();
  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();

  vec<sample_t> null_samples;
  sample_null(dim, params.sample_size, null_samples);
  REQUIRE(null_samples.size() >= GammaModel::min_nsamples);

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(4, 8, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(v.data(), u, d_obs);

  vec<record_t> records;
  records.emplace_back(0, 32, interval_t{1, 20}, interval_t{5, 9}, false, d_obs, I_obs, 0);
  records.emplace_back(0, 32, interval_t{10, 25}, interval_t{12, 16}, false, d_obs, I_obs, 0);
  records.emplace_back(0, 32, interval_t{1, 20}, interval_t{6, 10}, true, d_obs, I_obs, 0);
  for (auto& r : records) {
    r.d_diff = -0.05;
    r.d_q = 0.2;
    CHECK(test_significance(r, null_samples, params.sample_size, "q0"));
  }

  benjamini_hochberg_correction(records);

  for (const auto& r : records) {
    CHECK(std::isfinite(r.qvalue));
    CHECK(r.qvalue >= r.percentile);
    CHECK(r.qvalue <= 1.0);
  }
}

TEST_CASE("overlapping null windows on same query are excluded from scoring") {
  seed = 19;
  init_thread_rng(0);

  auto params = make_params(20);
  auto llhf = make_map_llhf();
  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();

  vec<sample_t> null_samples;
  sample_null(dim, params.sample_size, null_samples);
  REQUIRE(null_samples.size() >= GammaModel::min_nsamples);

  const sample_t& anchor = null_samples.front();
  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(anchor.bin_iv.a - 1, anchor.bin_iv.b - 1, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(v.data(), u, d_obs);

  record_t r(anchor.bix, 20, interval_t{1, 20}, anchor.bin_iv, false, d_obs, I_obs, 0);
  r.d_diff = 0.1;
  r.d_q = 0.15;

  CHECK(test_significance(r, null_samples, params.sample_size, "q0"));
  CHECK(std::isfinite(r.percentile));
}

TEST_CASE("filter_null_samples drops overlaps and invalid entries") {
  const vec<sample_t> pool{
    {0.10, 1.0, 0, {1, 3}},  // no overlap with {4, 7}
    {0.20, 1.0, 0, {4, 7}},  // overlaps record bin_iv
    {nanx(), 1.0, 0, {10, 14}},
    {0.30, nanx(), 0, {10, 14}},
    {0.40, 1.0, 1, {4, 7}}, // other query: keep despite coord overlap
    {0.50, 2.0, 0, {10, 14}},
  };
  record_t r(0, 32, interval_t{1, 20}, interval_t{4, 7}, false, 0.1, 1.0, 0);
  vec<sample_t> out;
  CHECK(filter_null_samples(pool, r, 100, out));
  REQUIRE(out.size() == 3);
  CHECK(out[0].d == doctest::Approx(0.10));
  CHECK(out[1].d == doctest::Approx(0.40));
  CHECK(out[2].d == doctest::Approx(0.50));
}

TEST_CASE("filter_null_samples reservoirs down to sample_size") {
  seed = 99;
  init_thread_rng(0);

  vec<sample_t> pool;
  for (uint64_t i = 0; i < 20; ++i)
    pool.push_back({0.1 + 0.01 * static_cast<double>(i), 1.0, 0, {100 + i, 104 + i}});
  record_t r(0, 32, interval_t{1, 20}, interval_t{1, 3}, false, 0.1, 1.0, 0);
  vec<sample_t> out;
  CHECK_FALSE(filter_null_samples(pool, r, 5, out));
  CHECK(out.size() == 5);
  for (const auto& s : out) {
    CHECK(std::isfinite(s.d));
    CHECK(std::isfinite(s.I));
  }
}

TEST_CASE("test_significance reuses gamma_fit_cache when no overlaps excluded") {
  vec<sample_t> null_samples;
  for (uint64_t i = 0; i < 16; ++i)
    null_samples.push_back({0.05 + 0.01 * static_cast<double>(i % 5), 5.0, 0, {20 + i, 23 + i}});

  record_t r0(0, 32, interval_t{1, 20}, interval_t{1, 4}, false, 0.08, 5.0, 0);
  record_t r1 = r0;
  r0.d_diff = nanx();
  r1.d_diff = nanx();

  gamma_fit_cache_t cache;
  CHECK(test_significance(r0, null_samples, null_samples.size(), "q0", &cache));
  REQUIRE(cache.ok);
  CHECK(cache.bix == 0);
  const auto cached = cache.params;
  const double p0 = r0.percentile;

  CHECK(test_significance(r1, null_samples, null_samples.size(), "q0", &cache));
  CHECK(cache.params.shape == doctest::Approx(cached.shape));
  CHECK(cache.params.scale == doctest::Approx(cached.scale));
  CHECK(r1.percentile == doctest::Approx(p0));
}

TEST_CASE("test_significance fails when null pool is too small") {
  auto params = make_params(20);
  vec<sample_t> null_samples;
  record_t r(0, 32, interval_t{1, 20}, interval_t{5, 9}, false, 0.1, 10.0, 0);
  r.d_diff = -0.05;
  r.d_q = 0.2;

  CHECK_FALSE(test_significance(r, null_samples, params.sample_size, "q0"));
  CHECK(std::isnan(r.percentile));
  CHECK(std::isnan(r.fold));
}

TEST_CASE("canonical records use one-sided test when d_diff is NaN") {
  seed = 11;
  init_thread_rng(0);

  auto params = make_params(20);
  auto llhf = make_map_llhf();
  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();

  vec<sample_t> null_samples;
  sample_null(dim, params.sample_size, null_samples);
  REQUIRE(null_samples.size() >= GammaModel::min_nsamples);

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(4, 8, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(v.data(), u, d_obs);

  record_t r(0, 32, interval_t{1, 20}, interval_t{5, 9}, false, d_obs, I_obs, 0);
  r.d_diff = nanx();

  CHECK(test_significance(r, null_samples, params.sample_size, "q0"));
  CHECK(std::isfinite(r.percentile));
  CHECK(r.percentile >= 0.0);
  CHECK(r.percentile <= 1.0);
  const bool two_sided = !std::isnan(r.d_diff) && (r.is_rc == (r.d_diff > 0.0));
  CHECK_FALSE(two_sided);
}

TEST_CASE("reference strand d_diff zero uses two-sided test") {
  seed = 11;
  init_thread_rng(0);

  auto params = make_params(20);
  auto llhf = make_map_llhf();
  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();

  vec<sample_t> null_samples;
  sample_null(dim, params.sample_size, null_samples);
  REQUIRE(null_samples.size() >= GammaModel::min_nsamples);

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(4, 8, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(v.data(), u, d_obs);

  record_t r(0, 32, interval_t{1, 20}, interval_t{5, 9}, false, d_obs, I_obs, 0);
  r.d_diff = 0.0;

  CHECK(test_significance(r, null_samples, params.sample_size, "q0"));
  CHECK(std::isfinite(r.percentile));
  CHECK((!std::isnan(r.d_diff) && (r.is_rc == (r.d_diff > 0.0))));
}

TEST_CASE("query strand uses one-sided test") {
  seed = 11;
  init_thread_rng(0);

  auto params = make_params(20);
  auto llhf = make_map_llhf();
  DIM<double> dim(params, llhf, 64, 64);
  inject_hits_all_bins(dim, 64);
  dim.compute_prefhistsum();

  vec<sample_t> null_samples;
  sample_null(dim, params.sample_size, null_samples);
  REQUIRE(null_samples.size() >= GammaModel::min_nsamples);

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(4, 8, v, u, t);
  const double d_obs = llhf->mle(v.data(), u);
  const double I_obs = llhf->compute_fisher_info(v.data(), u, d_obs);

  record_t r(0, 32, interval_t{1, 20}, interval_t{5, 9}, true, d_obs, I_obs, 0);
  r.d_diff = -0.1;
  r.d_q = 0.15;

  CHECK(test_significance(r, null_samples, params.sample_size, "q0"));
  CHECK(std::isfinite(r.percentile));
  const bool two_sided = !std::isnan(r.d_diff) && (r.is_rc == (r.d_diff > 0.0));
  CHECK_FALSE(two_sided);
}

TEST_CASE("benjamini_hochberg_correction with canonical-only records") {
  vec<record_t> records;
  records.emplace_back(0, 100, interval_t{1, 50}, interval_t{1, 5}, false, 0.1, 10.0, 0);
  records.emplace_back(0, 100, interval_t{51, 100}, interval_t{6, 10}, false, 0.12, 10.0, 0);
  records[0].percentile = 0.05;
  records[1].percentile = 0.10;

  benjamini_hochberg_correction(records);

  CHECK(records[0].qvalue == doctest::Approx(0.10));
  CHECK(records[1].qvalue == doctest::Approx(0.10));
}

TEST_CASE("benjamini_hochberg_correction leaves NaN qvalues untouched") {
  vec<record_t> records;
  records.emplace_back(0, 100, interval_t{1, 50}, interval_t{1, 5}, false, 0.1, 10.0, 0);
  records.emplace_back(0, 100, interval_t{51, 100}, interval_t{6, 10}, false, 0.12, 10.0, 0);
  records[0].percentile = 0.05;
  records[1].percentile = nanx();

  benjamini_hochberg_correction(records);

  CHECK(records[0].qvalue == doctest::Approx(0.05));
  CHECK(std::isnan(records[1].qvalue));
}

} // TEST_SUITE

TEST_SUITE("DistanceSampler") {

TEST_CASE("places sample_size windows and collect_distances keeps finite only") {
  seed = 42;
  init_thread_rng(0);

  const sketch_sptr_t sketch = load_tiny_sketch(true);
  const uint32_t k = sketch->get_lshf_sptr()->get_k();
  constexpr uint64_t tau = 50;
  constexpr uint64_t sample_size = 40;
  // Long enough for full windows: L >= nwinmers + k - 1 with bin_shift=0.
  vec<qseq_t> batch_v{{"q0", str(tau + k + 20, 'A')}, {"q_short", str(k + 5, 'C')}};

  ThreadPool pool(1);
  DistanceSampler sampler(sketch, batch_v, tau, 0, 4);
  sampler.run(sample_size, false, pool);

  CHECK(sampler.get_nsamples() == sample_size);
  CHECK(sampler.get_nwinmers() == tau);

  uint64_t n_all = 0, n_finite = 0, n_on_short = 0;
  sampler.for_each_sample([&](uint64_t bix, uint64_t, uint64_t, double d, char) {
    ++n_all;
    if (bix == 1) ++n_on_short;
    if (std::isfinite(d)) ++n_finite;
  });
  CHECK(n_all == sample_size);
  CHECK(n_on_short == 0);

  vec<double> d_v;
  sampler.collect_distances(d_v);
  CHECK(d_v.size() == n_finite);
  for (const double d : d_v)
    CHECK(std::isfinite(d));

  vec<vec<double>> d_vvec(batch_v.size());
  sampler.collect_distances(d_vvec);
  CHECK(d_vvec[0].size() == n_finite);
  CHECK(d_vvec[1].empty());
  CHECK(d_vvec[0].size() + d_vvec[1].size() == d_v.size());
}

TEST_CASE("skips all sequences shorter than the window") {
  seed = 1;
  init_thread_rng(0);

  const sketch_sptr_t sketch = load_tiny_sketch(true);
  const uint32_t k = sketch->get_lshf_sptr()->get_k();
  vec<qseq_t> batch_v{{"tiny", str(k + 5, 'A')}};

  ThreadPool pool(1);
  DistanceSampler sampler(sketch, batch_v, 50, 0, 4);
  sampler.run(30, false, pool);

  CHECK(sampler.get_nsamples() == 0);
  vec<double> d_v;
  sampler.collect_distances(d_v);
  CHECK(d_v.empty());
}

TEST_CASE("bin_shift rounds window length up to whole bins") {
  seed = 2;
  init_thread_rng(0);

  const sketch_sptr_t sketch = load_tiny_sketch(true);
  const uint32_t k = sketch->get_lshf_sptr()->get_k();
  constexpr uint64_t tau = 50;
  constexpr uint64_t bin_shift = 2; // bin_size = 4
  // tau_bin = ceil(50/4) = 13, nwinmers = 52
  constexpr uint64_t expect_nwinmers = 52;
  vec<qseq_t> batch_v{{"q0", str(expect_nwinmers + k + 20, 'A')}};

  ThreadPool pool(1);
  DistanceSampler sampler(sketch, batch_v, tau, bin_shift, 4);
  sampler.run(16, false, pool);

  CHECK(sampler.get_nwinmers() == expect_nwinmers);
  CHECK(sampler.get_nsamples() == 16);
  sampler.for_each_sample([&](uint64_t, uint64_t enmers, uint64_t start_bin, double, char) {
    const uint64_t jx = start_bin << bin_shift;
    CHECK(jx + expect_nwinmers <= enmers);
  });
}

TEST_CASE("canonical samples always report strand '.'") {
  seed = 3;
  init_thread_rng(0);

  const sketch_sptr_t sketch = load_tiny_sketch(true);
  const uint32_t k = sketch->get_lshf_sptr()->get_k();
  constexpr uint64_t tau = 40;
  vec<qseq_t> batch_v{{"q0", str(tau + k + 10, 'A')}};

  ThreadPool pool(1);
  DistanceSampler sampler(sketch, batch_v, tau, 0, 4);
  sampler.run(12, false, pool);
  REQUIRE(sampler.get_nsamples() == 12);

  sampler.for_each_sample([&](uint64_t, uint64_t, uint64_t, double, char strand) { CHECK(strand == '.'); });
}

TEST_CASE("non-canonical samples report +, -, or .") {
  seed = 4;
  init_thread_rng(0);

  const sketch_sptr_t sketch = load_tiny_sketch(false);
  const uint32_t k = sketch->get_lshf_sptr()->get_k();
  constexpr uint64_t tau = 40;
  vec<qseq_t> batch_v{{"q0", str(tau + k + 10, 'A')}};

  ThreadPool pool(1);
  DistanceSampler sampler(sketch, batch_v, tau, 0, 4);
  sampler.run(12, true, pool);
  REQUIRE(sampler.get_nsamples() == 12);
  CHECK(sampler.get_llhf().k == k);

  sampler.for_each_sample([&](uint64_t, uint64_t, uint64_t, double, char strand) {
    CHECK((strand == '+' || strand == '-' || strand == '.'));
  });
}

} // TEST_SUITE

TEST_SUITE("dist helpers") {

TEST_CASE("select_strand_distance prefers finite lower distance") {
  {
    const auto [d, s] = select_strand_distance(0.1, 0.2);
    CHECK(d == doctest::Approx(0.1));
    CHECK(s == '+');
  }
  {
    const auto [d, s] = select_strand_distance(0.3, 0.1);
    CHECK(d == doctest::Approx(0.1));
    CHECK(s == '-');
  }
  {
    const auto [d, s] = select_strand_distance(0.2, 0.2);
    CHECK(d == doctest::Approx(0.2));
    CHECK(s == '+');
  }
  {
    const auto [d, s] = select_strand_distance(nanx(), 0.2);
    CHECK(d == doctest::Approx(0.2));
    CHECK(s == '-');
  }
  {
    const auto [d, s] = select_strand_distance(0.2, nanx());
    CHECK(d == doctest::Approx(0.2));
    CHECK(s == '+');
  }
  {
    const auto [d, s] = select_strand_distance(nanx(), nanx());
    CHECK(std::isnan(d));
    CHECK(s == '.');
  }
}

TEST_CASE("linear_quantile on sorted values") {
  CHECK(std::isnan(linear_quantile({}, 0.5)));
  CHECK(linear_quantile({0.4}, 0.5) == doctest::Approx(0.4));
  const vec<double> v{0.0, 0.5, 1.0};
  CHECK(linear_quantile(v, 0.0) == doctest::Approx(0.0));
  CHECK(linear_quantile(v, 1.0) == doctest::Approx(1.0));
  CHECK(linear_quantile(v, 0.5) == doctest::Approx(0.5));
  CHECK(linear_quantile({0.0, 10.0}, 0.25) == doctest::Approx(2.5));
}

TEST_CASE("validate_binning rejects oversized bins") {
  CHECK(validate_binning(0, 10));
  CHECK(validate_binning(3, 8));
  CHECK_FALSE(validate_binning(4, 8)); // bin_size=16 > tau
  CHECK_FALSE(validate_binning(17, 100));
}

} // TEST_SUITE

TEST_SUITE("bracket_distance") {

TEST_CASE("NaN input gives the full range") {
  const vec<double> th_v{0.1, 0.2, 0.3};
  const auto [lo, hi] = bracket_distance(nanx(), th_v);
  CHECK(lo == doctest::Approx(d_eps));
  CHECK(hi == doctest::Approx(d_ub));
}

TEST_CASE("distance below all thresholds") {
  const vec<double> th_v{0.1, 0.2, 0.3};
  const auto [lo, hi] = bracket_distance(0.05, th_v);
  CHECK(lo == doctest::Approx(d_eps));
  CHECK(hi == doctest::Approx(0.1));
}

TEST_CASE("distance above all thresholds") {
  const vec<double> th_v{0.1, 0.2, 0.3};
  const auto [lo, hi] = bracket_distance(0.5, th_v);
  CHECK(lo == doctest::Approx(0.3));
  CHECK(hi == doctest::Approx(d_ub));
}

TEST_CASE("distance between thresholds") {
  const vec<double> th_v{0.1, 0.2, 0.3};
  const auto [lo, hi] = bracket_distance(0.15, th_v);
  CHECK(lo == doctest::Approx(0.1));
  CHECK(hi == doctest::Approx(0.2));
}

TEST_CASE("distance exactly at a threshold") {
  const vec<double> th_v{0.1, 0.2, 0.3};
  const auto [lo, hi] = bracket_distance(0.2, th_v);
  CHECK(lo == doctest::Approx(0.1));
  CHECK(hi == doctest::Approx(0.2));
}

TEST_CASE("empty thresholds give the full range") {
  const vec<double> th_v{};
  const auto [lo, hi] = bracket_distance(0.15, th_v);
  CHECK(lo == doctest::Approx(d_eps));
  CHECK(hi == doctest::Approx(d_ub));
}

} // TEST_SUITE
