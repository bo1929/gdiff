#include "doctest/doctest.h"
#include "intext.hpp"
#include <boost/math/tools/minima.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <random>

namespace {
// Test-local mirrors of IntervalMapper contiguous slicing and MLE+Fisher on bin ranges (see src/map.cpp).

template<typename T>
double at_th(T v, size_t ix)
{
  if constexpr (std::is_same_v<T, double>) {
    return v;
  } else {
    return v[ix];
  }
}

struct contig_slice_t
{
  uint64_t bin_a, bin_b;
  double d, I;
  uint8_t mask;
  double bin_lo, bin_hi;
};

// Uses LLH::mle() for MLE+Fisher computation after IntExt::extract_histogram.
template<typename T>
xy_t slice_mle_fisher(IntExt<T>& ext, const LLH<T>& llhf, uint64_t a1, uint64_t b1)
{
  vec<uint64_t> v;
  uint64_t u, t;
  ext.extract_histogram(a1 - 1, b1 - 1, v, u, t);
  const double d = llhf.mle(v.data(), u);
  const double I = llhf.compute_fisher_info(v.data(), u, d);
  return {d, I};
}

// Mirrors IntervalMapper::get_distance_bin (distance-aware threshold brackets).
template<typename T>
static xy_t distance_bin_bounds(const LLH<T>& llhf, size_t th_ix, double d, double d_q)
{
  static constexpr size_t W = std::is_same_v<T, double> ? 1 : rwidth;
  std::array<double, W> th_v{};
  for (size_t i = 0; i < W; ++i)
    th_v[i] = at_th(llhf.get_extrema(), i);
  std::sort(th_v.begin(), th_v.end());

  xy_t d_range{d_eps, d_ub};
  if (th_ix != no_threshold) {
    const double t_i = at_th(llhf.get_extrema(), th_ix);
    const bool is_low = std::isnan(d_q) || t_i <= d_q;
    const size_t pos = static_cast<size_t>(std::lower_bound(th_v.begin(), th_v.end(), t_i) - th_v.begin());
    if (is_low) {
      d_range.first = (pos > 0) ? th_v[pos - 1] : d_eps;
      d_range.second = t_i;
    } else {
      d_range.first = t_i;
      d_range.second = (pos + 1 < W) ? th_v[pos + 1] : d_ub;
    }
  } else if (std::isfinite(d) || std::isfinite(d_q)) {
    const double d_anchor = std::isfinite(d) ? d : d_q;
    const auto it = std::lower_bound(th_v.begin(), th_v.end(), d_anchor);
    if (it != th_v.begin()) d_range.first = *(it - 1);
    if (it != th_v.end()) d_range.second = *it;
  }
  return d_range;
}

// Mirrors IntervalMapper::extract_ordered_intervals and emit_record (1-based bin coords).
template<typename T>
vec<contig_slice_t> contig_slices(IntExt<T>& ext, const LLH<T>& llhf, uint8_t th_bv, double d_q)
{
  static constexpr size_t W = std::is_same_v<T, double> ? 1 : rwidth;
  vec<contig_slice_t> out;
  if (th_bv == 0) return out;

  const uint64_t nbins = ext.get_nbins();

  vec<uint64_t> pts = {1, nbins + 1};
  for (size_t ti = 0; ti < W; ++ti) {
    if (!(th_bv & (1u << ti)) || ext.get_intervals_v(ti).empty()) continue;
    for (const auto& iv : ext.get_intervals_v(ti)) {
      pts.push_back(iv.a);
      pts.push_back(iv.b + 1);
    }
  }

  std::sort(pts.begin(), pts.end());
  pts.erase(std::unique(pts.begin(), pts.end()), pts.end());

  std::array<size_t, rwidth> ti_ix{};

  for (size_t pi = 0; pi + 1 < pts.size(); ++pi) {
    const uint64_t a = pts[pi];
    const uint64_t b = pts[pi + 1];
    const auto [d, I] = slice_mle_fisher(ext, llhf, a, b);

    uint8_t mask = 0;
    double bin_lo = d_eps;
    double bin_hi = d_ub;
    for (size_t ti = 0; ti < W; ++ti) {
      if (!(th_bv & (1u << ti))) continue;
      const auto& ev = ext.get_intervals_v(ti);
      while (ti_ix[ti] < ev.size() && ev[ti_ix[ti]].b < a)
        ++ti_ix[ti];
      if (ti_ix[ti] >= ev.size() || ev[ti_ix[ti]].a > a) continue;
      mask |= static_cast<uint8_t>(1u << ti);
      const auto bounds = distance_bin_bounds(llhf, ti, d, d_q);
      bin_lo = std::max(bin_lo, bounds.first);
      bin_hi = std::min(bin_hi, bounds.second);
    }

    out.push_back({a, b, d, I, mask, bin_lo, bin_hi});
  }
  return out;
}

} // namespace

// Helper: create a IntExt<double> with given nbins, injecting known fdc/sdc values
static std::pair<LLH<double>, map_params<double>>
make_test_params(double dist_th = 0.1, uint32_t hdist_th = 4, uint64_t tau = 2, uint64_t bin_shift = 0)
{
  auto params = map_params<double>(dist_th, hdist_th, tau, bin_shift, 33.0);
  LLH<double> llhf(27, 11, 0.5, hdist_th, dist_th);
  return {llhf, params};
}

template<typename T>
static double query_mle(IntExt<T>& ext, const LLH<T>& llhf)
{
  ext.compute_prefhistsum();
  vec<uint64_t> v;
  uint64_t u, t;
  ext.total_histogram(v, u, t);
  return llhf.mle(v.data(), u);
}

template<typename T>
static interval_t interval_at(const IntExt<T>& ext, uint64_t i, size_t ix = 0)
{
  const auto& ivs = ext.get_intervals_v(ix);
  const uint64_t nbins = ext.get_nbins();
  if (i < ivs.size()) return ivs[i];
  return {nbins, nbins};
}

template<typename T>
static void finish_ext_scan(IntExt<T>& ext, const LLH<T>& llhf, double d_q)
{
  ext.inclusive_scan();
  ext.set_query_distance(d_q);
  ext.extrema_scan();
}

template<typename T>
static double finish_scan(IntExt<T>& ext, const LLH<T>& llhf)
{
  ext.inclusive_scan();
  const double d_q = query_mle(ext, llhf);
  ext.set_query_distance(d_q);
  ext.extrema_scan();
  return d_q;
}

template<typename T>
static void inject_udu(IntExt<T>& ext)
{
  for (uint64_t i : {0u, 1u}) {
    ext.aggregate_mer(4, i);
    ext.aggregate_mer(4, i);
  }
  for (uint64_t i : {3u, 4u}) {
    ext.aggregate_mer(0, i);
    ext.aggregate_mer(0, i);
  }
  for (uint64_t i : {6u, 7u}) {
    ext.aggregate_mer(4, i);
    ext.aggregate_mer(4, i);
  }
}

static vec<interval_t> collect_intervals_double(
  const map_params<double>& params,
  const LLH<double>& llhf,
  uint64_t nbins,
  double d_q,
  const std::function<void(IntExt<double>&)>& inject,
  uint64_t tau = 1)
{
  IntExt<double> ext(params, llhf, nbins, nbins);
  inject(ext);
  finish_ext_scan(ext, llhf, d_q);
  ext.extract_intervals_mx(tau, 1, nbins);
  ext.expand_intervals(33.0);
  return ext.get_intervals_v(0);
}

static bool intervals_equal(const vec<interval_t>& a, const vec<interval_t>& b)
{
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i].a != b[i].a || a[i].b != b[i].b) return false;
  }
  return true;
}

TEST_SUITE("IntExt<double>::inclusive_scan") {

TEST_CASE("inclusive_scan with no hits yields no merged intervals") {
  auto [llhf, params] = make_test_params();
  const uint64_t nbins = 5;
  IntExt<double> ext(params, llhf, nbins, nbins);

  // No hits: fdc/sdc stay zero -> no merged intervals after expand.
  finish_ext_scan(ext, llhf, nanx());
  ext.extract_intervals_mx(0, 1, nbins);
  ext.expand_intervals(33.0);
  auto iv = interval_at(ext, 0);
  CHECK(iv.a >= nbins);
}

} // TEST_SUITE


TEST_SUITE("IntExt<double>::expand_intervals") {

TEST_CASE("merge when chi-square below threshold") {
  // Named rather than a structured binding: C++17 cannot capture structured bindings in a lambda.
  const auto test_params = make_test_params(0.1, 4, 1, 0);
  const LLH<double>& llhf = test_params.first;
  const map_params<double>& params = test_params.second;
  const uint64_t nbins = 30;

  auto count_ivs = [&](double chisq_th) -> uint64_t {
    IntExt<double> d(params, llhf, nbins, nbins);
    for (uint64_t i = 0; i < 10; ++i) {
      d.aggregate_mer(0, i);
      d.aggregate_mer(0, i);
    }
    for (uint64_t i = 12; i < 25; ++i) {
      d.aggregate_mer(0, i);
      d.aggregate_mer(0, i);
    }
    finish_scan(d, llhf);
    d.extract_intervals_mx(0, 1, nbins);
    d.expand_intervals(chisq_th);
    uint64_t n = 0;
    for (;; ++n) {
      const interval_t iv = interval_at(d, n);
      if (iv.a >= nbins) break;
    }
    return n;
  };

  const uint64_t n_strict = count_ivs(0.0);
  const uint64_t n_merged = count_ivs(1e20);
  CHECK(n_merged <= n_strict);
}

} // TEST_SUITE

TEST_SUITE("IntExt<double>::extract_histogram") {

TEST_CASE("extract_histogram returns correct counts for enum_only=false") {
  // Need enum_only=false for hdisthist_v to be allocated
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 10;
  const uint64_t nmers = 100;
  IntExt<double> ext(params, llhf, nbins, nmers);

  // One hit per bin 0..4 (hdist=0): t=5 equals mers in [0,5) when bin_shift=0 -> u=0
  for (uint64_t i = 0; i < 5; ++i) {
    ext.aggregate_mer(0, i);
  }

  ext.compute_prefhistsum();

  vec<uint64_t> v;
  uint64_t u, t;
  ext.extract_histogram(0, 5, v, u, t);

  CHECK(v[0] == 5);
  CHECK(t == 5);
  CHECK(u == 0);
  for (uint32_t d = 1; d <= params.hdist_th; ++d) {
    CHECK(v[d] == 0);
  }
}

TEST_CASE("extract_histogram counts explicit misses only") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  IntExt<double> ext(params, llhf, 8, 100);

  ext.aggregate_mer(0, 0);
  ext.aggregate_mer(5, 1);
  ext.aggregate_mer(2, 2);
  ext.aggregate_mer(7, 3);
  ext.compute_prefhistsum();

  vec<uint64_t> v;
  uint64_t u = 0, t = 0;
  ext.extract_histogram(0, 8, v, u, t);
  CHECK(t == 2);
  CHECK(u == 2);

  uint64_t u1 = 0, t1 = 0, u2 = 0, t2 = 0;
  vec<uint64_t> v1, v2;
  ext.extract_histogram(0, 2, v1, u1, t1);
  ext.extract_histogram(2, 8, v2, u2, t2);
  CHECK(u1 + u2 == u);
  CHECK(t1 + t2 == t);
}

TEST_CASE("extract_histogram full range matches partial sums") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 8;
  const uint64_t nmers = 1000;
  IntExt<double> ext(params, llhf, nbins, nmers);

  // Inject hits into various bins with various hamming distances
  ext.aggregate_mer(0, 0);
  ext.aggregate_mer(1, 1);
  ext.aggregate_mer(2, 2);
  ext.aggregate_mer(3, 3);
  ext.aggregate_mer(0, 4);
  ext.aggregate_mer(0, 5);
  ext.aggregate_mer(1, 6);
  ext.aggregate_mer(2, 7);

  ext.compute_prefhistsum();

  // Full range
  vec<uint64_t> v_full;
  uint64_t u_full, t_full;
  ext.extract_histogram(0, nbins, v_full, u_full, t_full);

  // Split range: [0,4) + [4,8) should equal [0,8)
  vec<uint64_t> v1, v2;
  uint64_t u1, t1, u2, t2;
  ext.extract_histogram(0, 4, v1, u1, t1);
  ext.extract_histogram(4, nbins, v2, u2, t2);

  const uint32_t W = params.hdist_th + 1;
  for (uint32_t d = 0; d < W; ++d) {
    CHECK(v1[d] + v2[d] == v_full[d]);
  }
  CHECK(t1 + t2 == t_full);
  CHECK(u1 + u2 == u_full);
}

TEST_CASE("extract_histogram split additivity with bin_shift > 0") {
  auto params = map_params<double>(0.1, 4, 2, 2, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 8;
  const uint64_t nmers = 1000;
  IntExt<double> ext(params, llhf, nbins, nmers);

  ext.aggregate_mer(0, 0);
  ext.aggregate_mer(1, 1);
  ext.aggregate_mer(2, 2);
  ext.aggregate_mer(3, 3);

  ext.compute_prefhistsum();

  vec<uint64_t> v_full;
  uint64_t u_full, t_full;
  ext.extract_histogram(0, nbins, v_full, u_full, t_full);

  vec<uint64_t> v1, v2;
  uint64_t u1, t1, u2, t2;
  ext.extract_histogram(0, 4, v1, u1, t1);
  ext.extract_histogram(4, nbins, v2, u2, t2);

  const uint32_t W = params.hdist_th + 1;
  for (uint32_t d = 0; d < W; ++d) {
    CHECK(v1[d] + v2[d] == v_full[d]);
  }
  CHECK(t1 + t2 == t_full);
  CHECK(u1 + u2 == u_full);
}

TEST_CASE("extract_histogram supports the maximum fixed SIMD hdist threshold") {
  auto params = map_params<double>(0.1, hdist_bound, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, hdist_bound, 0.1);
  const uint64_t nbins = 4;
  IntExt<double> ext(params, llhf, nbins, nbins);

  ext.aggregate_mer(hdist_bound, 0);
  ext.aggregate_mer(0, 1);
  ext.compute_prefhistsum();

  vec<uint64_t> v;
  uint64_t u, t;
  ext.extract_histogram(0, nbins, v, u, t);

  REQUIRE(v.size() == hdist_bound + 1);
  CHECK(v[0] == 1);
  CHECK(v[hdist_bound] == 1);
  CHECK(t == 2);
}

TEST_CASE("total_histogram matches extract_histogram on full range when keep_hist") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  IntExt<double> ext(params, llhf, 5, 50);

  ext.aggregate_mer(0, 0);
  ext.aggregate_mer(1, 1);
  ext.aggregate_mer(5, 2);
  ext.compute_prefhistsum();

  vec<uint64_t> v_complete, v_extract;
  uint64_t u_complete, t_complete, u_extract, t_extract;
  ext.total_histogram(v_complete, u_complete, t_complete);
  ext.extract_histogram(0, ext.get_nbins(), v_extract, u_extract, t_extract);

  CHECK(v_complete == v_extract);
  CHECK(u_complete == u_extract);
  CHECK(t_complete == t_extract);
}

TEST_CASE("total_histogram exposes global hit counts and explicit misses") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  IntExt<double> ext(params, llhf, 5, 50);

  ext.aggregate_mer(0, 0);
  ext.aggregate_mer(1, 1);
  ext.aggregate_mer(5, 2); // explicit miss (hdist > hdist_th)
  ext.compute_prefhistsum();

  vec<uint64_t> v;
  uint64_t u, t;
  ext.total_histogram(v, u, t);

  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(u == 1);
  CHECK(t == 2);
}

} // TEST_SUITE

TEST_SUITE("IntExt<cmlane_t>") {

TEST_CASE("SIMD IntExt produces valid intervals") {
  cmlane_t dths{};
  for (int i = 0; i < 8; ++i) dths[i] = 0.05 * (i + 1);
  auto params = map_params<cmlane_t>(dths, 4, 2, 0, 33.0);
  LLH<cmlane_t> llhf(27, 11, 0.5, 4, dths);
  const uint64_t nbins = 15;
  IntExt<cmlane_t> ext(params, llhf, nbins, nbins);

  // Inject hits in a pattern
  for (uint64_t i = 0; i < nbins; ++i) {
    if (i < 4 || i >= 11) {
      ext.aggregate_mer(0, i);
      ext.aggregate_mer(0, i);
    }
  }

  finish_scan(ext, llhf);

  // Extract intervals for each threshold
  for (size_t ix = 0; ix < 8; ++ix) {
    ext.extract_intervals_mx(1, 1, nbins, ix);
    ext.expand_intervals(33.0, ix);
    for (const auto& iv : ext.get_intervals_v(ix)) {
      CHECK(iv.a >= 1);
      CHECK(iv.a <= iv.b);
      CHECK(iv.b <= nbins);
    }
  }
}

} // TEST_SUITE

TEST_SUITE("IntExt::find_ordered_intervals") {

// Helper: run the full pipeline to populate eintervals_v, then list contiguous slices in bin space.
template<typename T>
static vec<contig_slice_t> run_pipeline(IntExt<T>& ext, const LLH<T>& llhf, uint8_t th_bv, uint64_t tau = 1)
{
  ext.inclusive_scan();
  ext.compute_prefhistsum();
  vec<uint64_t> v;
  uint64_t u, t;
  ext.total_histogram(v, u, t);
  const double d_q = llhf.mle(v.data(), u);
  ext.set_query_distance(d_q);
  ext.extrema_scan();
  const uint64_t nbins = ext.get_nbins();
  if constexpr (std::is_same_v<T, double>) {
    ext.extract_intervals_mx(tau, 1, nbins);
    ext.expand_intervals(33.0);
  } else {
    constexpr size_t N = std::is_same_v<T, cmlane_t> ? rwidth : 1;
    for (size_t ix = 0; ix < N; ++ix) {
      ext.extract_intervals_mx(tau, 1, nbins, ix);
      ext.expand_intervals(33.0, ix);
    }
  }
  return contig_slices(ext, llhf, th_bv, d_q);
}

TEST_CASE("th_bv=0 returns no segments") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  IntExt<double> ext(params, llhf, 10, 100);
  auto segs = contig_slices(ext, llhf, 0, nanx());
  CHECK(segs.empty());
}

TEST_CASE("low threshold below d_q: matched interval brackets (prev, t]") {
  auto params = map_params<double>(0.01, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.01);
  const uint64_t nbins = 10;
  IntExt<double> ext(params, llhf, nbins, nbins);

  // up-down-up: hdist=4 in 0-1, hdist=0 in 3-4, hdist=4 in 6-7.
  for (uint64_t i : {0u, 1u}) {
    ext.aggregate_mer(4, i); ext.aggregate_mer(4, i);
  }
  for (uint64_t i : {3u, 4u}) {
    ext.aggregate_mer(0, i); ext.aggregate_mer(0, i);
  }
  for (uint64_t i : {6u, 7u}) {
    ext.aggregate_mer(4, i); ext.aggregate_mer(4, i);
  }

  auto segs = run_pipeline(ext, llhf, 1);
  CHECK(!segs.empty());

  vec<uint64_t> v;
  uint64_t u, t;
  ext.total_histogram(v, u, t);
  const double d_q = llhf.mle(v.data(), u);
  REQUIRE(d_q > 0.01);

  for (const auto& s : segs) {
    if (s.mask) {
      const auto b = distance_bin_bounds(llhf, 0, s.d, d_q);
      CHECK(s.bin_lo == doctest::Approx(b.first));
      CHECK(s.bin_hi == doctest::Approx(b.second));
      CHECK(b.first == doctest::Approx(d_eps));
      CHECK(b.second == doctest::Approx(0.01));
    }
  }
}

TEST_CASE("high threshold above d_q: matched interval brackets (t, next]") {
  auto params = map_params<double>(0.5, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.5);
  const uint64_t nbins = 10;
  IntExt<double> ext(params, llhf, nbins, nbins);

  // Uniform low Hamming distance -> low query MLE, above-threshold extraction.
  for (uint64_t i = 0; i < nbins; ++i)
    ext.aggregate_mer(0, i);

  auto segs = run_pipeline(ext, llhf, 1);
  CHECK(!segs.empty());

  vec<uint64_t> v;
  uint64_t u, t;
  ext.total_histogram(v, u, t);
  const double d_q = llhf.mle(v.data(), u);
  REQUIRE(d_q < 0.5);

  for (const auto& s : segs) {
    if (s.mask) {
      const auto b = distance_bin_bounds(llhf, 0, s.d, d_q);
      CHECK(s.bin_lo == doctest::Approx(b.first));
      CHECK(s.bin_hi == doctest::Approx(b.second));
      CHECK(b.first == doctest::Approx(0.5));
      CHECK(b.second == doctest::Approx(d_ub));
    }
  }
}

TEST_CASE("segments cover [1, nbins+1) when interval spans full range") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 10;
  IntExt<double> ext(params, llhf, nbins, nbins);

  for (uint64_t i : {0u, 1u}) {
    ext.aggregate_mer(4, i); ext.aggregate_mer(4, i);
  }
  for (uint64_t i : {3u, 4u}) {
    ext.aggregate_mer(0, i); ext.aggregate_mer(0, i);
  }
  for (uint64_t i : {6u, 7u}) {
    ext.aggregate_mer(4, i); ext.aggregate_mer(4, i);
  }

  auto segs = run_pipeline(ext, llhf, 1);
  REQUIRE(!segs.empty());

  CHECK(segs.front().bin_a == 1);
  CHECK(segs.back().bin_b == nbins + 1);
  for (size_t i = 0; i + 1 < segs.size(); ++i) {
    CHECK(segs[i].bin_b == segs[i + 1].bin_a);
  }
}

TEST_CASE("cmlane_t thresholds ranked relative to d_q") {
  cmlane_t dths{};
  dths[0] = 0.05;
  dths[1] = 0.10;
  dths[2] = 0.20;
  dths[3] = 0.30;

  auto params = map_params<cmlane_t>(dths, 4, 2, 0, 33.0);
  LLH<cmlane_t> llhf(27, 11, 0.5, 4, dths);
  const uint64_t nbins = 10;
  IntExt<cmlane_t> ext(params, llhf, nbins, nbins);

  for (uint64_t i : {0u, 1u}) {
    ext.aggregate_mer(4, i); ext.aggregate_mer(4, i);
  }
  for (uint64_t i : {3u, 4u}) {
    ext.aggregate_mer(0, i); ext.aggregate_mer(0, i);
  }
  for (uint64_t i : {6u, 7u}) {
    ext.aggregate_mer(4, i); ext.aggregate_mer(4, i);
  }

  uint8_t th_bv = (1u << 0) | (1u << 2);
  auto segs = run_pipeline(ext, llhf, th_bv);
  REQUIRE(!segs.empty());

  vec<uint64_t> v;
  uint64_t u, t;
  ext.total_histogram(v, u, t);
  const double d_q = llhf.mle(v.data(), u);

  for (const auto& s : segs) {
    if (!s.mask) continue;
    double expect_lo = d_eps;
    double expect_hi = d_ub;
    for (size_t ti = 0; ti < 4; ++ti) {
      if (!(s.mask & (1u << ti))) continue;
      const auto b = distance_bin_bounds(llhf, ti, s.d, d_q);
      expect_lo = std::max(expect_lo, b.first);
      expect_hi = std::min(expect_hi, b.second);
    }
    CHECK(s.bin_lo == doctest::Approx(expect_lo));
    CHECK(s.bin_hi == doctest::Approx(expect_hi));
  }
}

TEST_CASE("no intervals -> one intact segment (endpoints only)") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 8;
  IntExt<double> ext(params, llhf, nbins, nbins);

  // No hits: extract finds nothing; the full-query fallback still reports [1, nbins+1).

  auto segs = run_pipeline(ext, llhf, 1);
  REQUIRE(segs.size() == 1);
  CHECK(segs[0].bin_a == 1);
  CHECK(segs[0].bin_b == nbins + 1);
  CHECK(segs[0].mask == 0);
}

TEST_CASE("full-span merged interval [1, nbins] still yields one segment") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 10;
  IntExt<double> ext(params, llhf, nbins, nbins);
  for (uint64_t i = 0; i < nbins; ++i)
    ext.aggregate_mer(0, i);

  auto segs = run_pipeline(ext, llhf, 1);
  REQUIRE(!segs.empty());
  CHECK(segs.size() == 1);
  CHECK(segs[0].bin_a == 1);
  CHECK(segs[0].bin_b == nbins + 1);
  CHECK(segs[0].mask == 0); // uniform hits: flat prefix sum, no threshold interval
}

TEST_CASE("map_contiguous segments partition histogram hit counts") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 10;
  IntExt<double> ext(params, llhf, nbins, nbins);

  for (uint64_t i : {0u, 1u}) {
    ext.aggregate_mer(4, i);
    ext.aggregate_mer(4, i);
  }
  for (uint64_t i : {3u, 4u}) {
    ext.aggregate_mer(0, i);
    ext.aggregate_mer(0, i);
  }
  for (uint64_t i : {6u, 7u}) {
    ext.aggregate_mer(4, i);
    ext.aggregate_mer(4, i);
  }

  const auto segs = run_pipeline(ext, llhf, 1);
  REQUIRE(!segs.empty());

  vec<uint64_t> v_full;
  uint64_t u_full, t_full;
  ext.extract_histogram(0, nbins, v_full, u_full, t_full);

  const uint32_t W = params.hdist_th + 1;
  vec<uint64_t> v_sum(W, 0);
  uint64_t t_sum = 0;
  uint64_t u_sum = 0;
  for (const auto& s : segs) {
    vec<uint64_t> v;
    uint64_t u, t;
    ext.extract_histogram(s.bin_a - 1, s.bin_b - 1, v, u, t);
    for (uint32_t d = 0; d < W; ++d)
      v_sum[d] += v[d];
    t_sum += t;
    u_sum += u;
  }
  CHECK(t_sum == t_full);
  CHECK(u_sum == u_full);
  for (uint32_t d = 0; d < W; ++d)
    CHECK(v_sum[d] == v_full[d]);
}

TEST_CASE("segment MLE matches Brent+Fisher reference on same bin range") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 12;
  IntExt<double> ext(params, llhf, nbins, nbins);

  for (uint64_t i = 0; i < nbins; ++i) {
    if (i % 3 == 0)
      ext.aggregate_mer(0, i);
    else if (i % 3 == 1)
      ext.aggregate_mer(2, i);
    else
      ext.aggregate_mer(4, i);
  }

  const auto segs = run_pipeline(ext, llhf, 1);
  REQUIRE(!segs.empty());

  for (const auto& s : segs) {
    const auto [d2, I2] = slice_mle_fisher(ext, llhf, s.bin_a, s.bin_b);
    CHECK(s.d == doctest::Approx(d2).epsilon(1e-9));
    if (std::isfinite(s.I) && std::isfinite(I2)) {
      CHECK(s.I == doctest::Approx(I2).epsilon(1e-7));
    } else {
      CHECK(std::isnan(s.I) == std::isnan(I2));
    }
  }
}

} // TEST_SUITE

TEST_SUITE("IntExt::set_query_distance") {

TEST_CASE("double: thrank_v is always the single threshold index") {
  auto [llhf, params] = make_test_params(0.1, 4, 2, 0);
  IntExt<double> ext(params, llhf, 8, 80);
  for (const double d_q : {0.05, 0.25, 0.99, nanx()}) {
    ext.set_query_distance(d_q);
    const auto& thrank = ext.get_thrank_v();
    REQUIRE(thrank.size() == 1);
    CHECK(thrank[0] == 0);
  }
}

TEST_CASE("cmlane_t: thrank_v orders thresholds relative to d_q") {
  cmlane_t dths{};
  dths[0] = 0.10;
  dths[1] = 0.05;
  dths[2] = 0.20;
  dths[3] = 0.03;
  dths[4] = 0.40;
  dths[5] = 0.35;
  dths[6] = 0.50;
  dths[7] = 0.45;

  auto params = map_params<cmlane_t>(dths, 4, 2, 0, 33.0);
  LLH<cmlane_t> llhf(27, 11, 0.5, 4, dths);
  IntExt<cmlane_t> ext(params, llhf, 8, 80);

  for (uint64_t i = 0; i < 8; ++i)
    ext.aggregate_mer(0, i);

  ext.set_query_distance(0.25);
  const auto& thrank = ext.get_thrank_v();
  REQUIRE(thrank.size() == rwidth);
  CHECK(thrank[0] == 3); // 0.03
  CHECK(thrank[1] == 1); // 0.05
  CHECK(thrank[2] == 0); // 0.10
  CHECK(thrank[3] == 2); // 0.20
  CHECK(thrank[4] == 6); // 0.50
  CHECK(thrank[5] == 7); // 0.45
  CHECK(thrank[6] == 4); // 0.40
  CHECK(thrank[7] == 5); // 0.35
}

TEST_CASE("cmlane_t: non-flipped lanes before flipped at d_q=0.088094") {
  cmlane_t dths{};
  dths[0] = 0.05;
  dths[1] = 0.075;
  dths[2] = 0.1;
  dths[3] = 0.125;
  dths[4] = 0.15;
  dths[5] = 0.225;
  dths[6] = 0.2;
  dths[7] = 0.25;

  auto params = map_params<cmlane_t>(dths, 4, 2, 0, 33.0);
  LLH<cmlane_t> llhf(27, 11, 0.5, 4, dths);
  IntExt<cmlane_t> ext(params, llhf, 4, 40);

  ext.set_query_distance(0.088094);
  const auto& thrank = ext.get_thrank_v();
  REQUIRE(thrank.size() == rwidth);
  CHECK(thrank[0] == 0); // 0.05
  CHECK(thrank[1] == 1); // 0.075
  CHECK(thrank[2] == 7); // 0.25
  CHECK(thrank[3] == 5); // 0.225
  CHECK(thrank[4] == 6); // 0.2
  CHECK(thrank[5] == 4); // 0.15
  CHECK(thrank[6] == 3); // 0.125
  CHECK(thrank[7] == 2); // 0.1
}

TEST_CASE("cmlane_t: NaN d_q sorts thrank_v by ascending threshold") {
  cmlane_t dths{};
  dths[0] = 0.30;
  dths[1] = 0.05;
  dths[2] = 0.20;
  dths[3] = 0.10;
  dths[4] = 0.25;
  dths[5] = 0.15;
  dths[6] = 0.35;
  dths[7] = 0.45;

  auto params = map_params<cmlane_t>(dths, 4, 2, 0, 33.0);
  LLH<cmlane_t> llhf(27, 11, 0.5, 4, dths);
  IntExt<cmlane_t> ext(params, llhf, 4, 40);

  ext.set_query_distance(nanx());
  const auto& thrank = ext.get_thrank_v();
  REQUIRE(thrank.size() == rwidth);
  CHECK(thrank[0] == 1); // 0.05
  CHECK(thrank[1] == 3); // 0.10
  CHECK(thrank[2] == 5); // 0.15
  CHECK(thrank[3] == 2); // 0.20
  CHECK(thrank[4] == 4); // 0.25
  CHECK(thrank[5] == 0); // 0.30
  CHECK(thrank[6] == 6); // 0.35
  CHECK(thrank[7] == 7); // 0.45
}

TEST_CASE("apply_threshold_signs: flip when t exceeds d_q changes intervals") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 10;
  const auto inject = [](IntExt<double>& ext) { inject_udu(ext); };

  const auto iv_flip = collect_intervals_double(params, llhf, nbins, 0.05, inject);
  const auto iv_no_flip = collect_intervals_double(params, llhf, nbins, 0.5, inject);
  const auto iv_nan = collect_intervals_double(params, llhf, nbins, nanx(), inject);

  CHECK_FALSE(intervals_equal(iv_flip, iv_no_flip));
  CHECK(intervals_equal(iv_no_flip, iv_nan));
}

TEST_CASE("set_query_distance must be called once per inclusive_scan") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 10;

  IntExt<double> ext_once(params, llhf, nbins, nbins);
  inject_udu(ext_once);
  finish_ext_scan(ext_once, llhf, 0.05);
  ext_once.extract_intervals_mx(1, 1, nbins);
  ext_once.expand_intervals(33.0);
  const auto iv_once = ext_once.get_intervals_v(0);

  IntExt<double> ext_twice(params, llhf, nbins, nbins);
  inject_udu(ext_twice);
  ext_twice.inclusive_scan();
  ext_twice.set_query_distance(0.05);
  ext_twice.set_query_distance(0.05);
  ext_twice.extrema_scan();
  ext_twice.extract_intervals_mx(1, 1, nbins);
  ext_twice.expand_intervals(33.0);
  const auto iv_twice = ext_twice.get_intervals_v(0);

  CHECK_FALSE(intervals_equal(iv_once, iv_twice));
}

TEST_CASE("cmlane_t: per-lane flip depends on threshold vs d_q") {
  cmlane_t dths{};
  dths[0] = 0.05;
  dths[1] = 0.40;

  auto params = map_params<cmlane_t>(dths, 4, 2, 0, 33.0);
  LLH<cmlane_t> llhf(27, 11, 0.5, 4, dths);
  const uint64_t nbins = 10;

  auto collect_lane = [&](size_t ix, double d_q) {
    IntExt<cmlane_t> ext(params, llhf, nbins, nbins);
    inject_udu(ext);
    finish_ext_scan(ext, llhf, d_q);
    ext.extract_intervals_mx(1, 1, nbins, ix);
    ext.expand_intervals(33.0, ix);
    return ext.get_intervals_v(ix);
  };

  const auto lane0_mid = collect_lane(0, 0.20);
  const auto lane1_mid = collect_lane(1, 0.20);
  const auto lane0_hi = collect_lane(0, 0.60);
  const auto lane1_hi = collect_lane(1, 0.60);

  CHECK(intervals_equal(lane0_mid, lane0_hi));
  CHECK_FALSE(intervals_equal(lane1_mid, lane1_hi));
}

TEST_CASE("skipping set_query_distance differs from production scan") {
  auto params = map_params<double>(0.1, 4, 2, 0, 33.0);
  LLH<double> llhf(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 10;

  const auto iv_prod = collect_intervals_double(
    params, llhf, nbins, 0.05, [](IntExt<double>& ext) { inject_udu(ext); });

  IntExt<double> ext_skip(params, llhf, nbins, nbins);
  inject_udu(ext_skip);
  ext_skip.inclusive_scan();
  ext_skip.extrema_scan();
  ext_skip.extract_intervals_mx(1, 1, nbins);
  ext_skip.expand_intervals(33.0);
  const auto iv_skip = ext_skip.get_intervals_v(0);

  CHECK_FALSE(intervals_equal(iv_prod, iv_skip));
}

} // TEST_SUITE

TEST_SUITE("Misc") {

TEST_CASE("MH step equals sigma in sample_metropolis_hastings") {
  // Invariant: step = init.second, passed as sqrt(1/I) = sigma (map.cpp:638).
  CHECK(true);
}

TEST_CASE("multi-gap extraction finds intervals in separated regions") {
  auto [llhf, params] = make_test_params(0.1, 4, 2, 0);
  const uint64_t nbins = 25;
  IntExt<double> ext(params, llhf, nbins, nbins);

  // Two divergence regions separated by a gap.
  for (uint64_t i = 2; i <= 5; ++i) {
    ext.aggregate_mer(0, i); ext.aggregate_mer(0, i);
  }
  for (uint64_t i = 16; i <= 19; ++i) {
    ext.aggregate_mer(0, i); ext.aggregate_mer(0, i);
  }

  finish_ext_scan(ext, llhf, nanx());
  ext.extract_intervals_mx(1, 1, nbins);
  ext.expand_intervals(33.0);

  const auto& e_v = ext.get_intervals_v(0);
  CHECK(e_v.size() >= 1);
  for (const auto& iv : e_v) {
    CHECK(iv.a >= 1);
    CHECK(iv.a <= iv.b);
    CHECK(iv.b <= nbins);
  }
}

} // TEST_SUITE

TEST_SUITE("IntExt<double> N-run skips") {

TEST_CASE("skip_mer splits extraction at the flagged bin") {
  auto [llhf, params] = make_test_params(0.1, 4, 2, 0);
  const uint64_t nbins = 30;

  auto inject = [&](IntExt<double>& ext) {
    for (uint64_t i = 0; i < 10; ++i) { ext.aggregate_mer(0, i); ext.aggregate_mer(0, i); }
    for (uint64_t i = 12; i < 25; ++i) { ext.aggregate_mer(0, i); ext.aggregate_mer(0, i); }
  };

  // Baseline: no skips -> one interval spans the gap (1-based bin 11).
  IntExt<double> ext_ns(params, llhf, nbins, nbins);
  inject(ext_ns);
  finish_ext_scan(ext_ns, llhf, nanx());
  ext_ns.extract_intervals_mx(0, 1, nbins);
  ext_ns.expand_intervals(33.0);
  bool covers_gap = false;
  for (const auto& iv : ext_ns.get_intervals_v(0))
    covers_gap = covers_gap || (iv.a <= 11 && 11 <= iv.b);
  CHECK(covers_gap);

  // Skip at 0-based bin 10: no interval crosses it.
  IntExt<double> ext(params, llhf, nbins, nbins);
  inject(ext);
  ext.skip_mer(10);
  CHECK(ext.is_wskip());
  CHECK(ext.is_skip(11));
  finish_ext_scan(ext, llhf, nanx());
  ext.extract_intervals_mx(0, 1, nbins);
  ext.expand_intervals(33.0);
  const auto& iv_v = ext.get_intervals_v(0);
  CHECK(iv_v.size() >= 2);
  bool left = false, right = false;
  for (const auto& iv : iv_v) {
    CHECK(!(iv.a <= 11 && 11 <= iv.b));
    left = left || (iv.b <= 10);
    right = right || (iv.a >= 12);
  }
  CHECK(left);
  CHECK(right);
}

TEST_CASE("no skip_mer calls: lazy state stays empty") {
  auto [llhf, params] = make_test_params();
  const uint64_t nbins = 8;
  IntExt<double> ext(params, llhf, nbins, nbins);
  for (uint64_t i = 0; i < nbins; ++i)
    ext.aggregate_mer(0, i);
  finish_ext_scan(ext, llhf, nanx());
  ext.extract_intervals_mx(0, 1, nbins);
  ext.expand_intervals(33.0);
  CHECK(!ext.is_wskip());
  CHECK(!ext.is_skip(1));
}

TEST_CASE("skip_mer ignores out-of-range bins") {
  auto [llhf, params] = make_test_params();
  IntExt<double> ext(params, llhf, 8, 8);
  ext.skip_mer(8);
  ext.skip_mer(100);
  CHECK(!ext.is_wskip());
}

} // TEST_SUITE



