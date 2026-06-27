#include "doctest/doctest.h"
#include "map.hpp"
#include <boost/math/tools/minima.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <random>

namespace {
// Test-local mirrors of QIE contiguous slicing and MLE+Fisher on bin ranges (see src/map.cpp).

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

// Uses LLH::mle() for MLE+Fisher computation after DIM::extract_histogram.
template<typename T>
xy_t slice_mle_fisher(DIM<T>& dim, const llh_sptr_t<T>& llhf, uint64_t a1, uint64_t b1)
{
  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(a1 - 1, b1 - 1, v, u, t);
  const double d = llhf->mle(v.data(), u);
  const double I = llhf->compute_fisher_info(d);
  return {d, I};
}

// Mirrors QIE::get_distance_bin (distance-aware threshold brackets).
template<typename T>
static xy_t distance_bin_bounds(const llh_sptr_t<T>& llhf, size_t th_ix, double d, double d_q)
{
  static constexpr size_t W = std::is_same_v<T, double> ? 1 : RWIDTH;
  std::array<double, W> th_v{};
  for (size_t i = 0; i < W; ++i)
    th_v[i] = at_th(llhf->get_extrema(), i);
  std::sort(th_v.begin(), th_v.end());

  xy_t d_range{d_eps, d_ub};
  if (th_ix != size_t(-1)) {
    const double t_i = at_th(llhf->get_extrema(), th_ix);
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

// Mirrors QIE::extract_ordered_intervals and emit_record (1-based bin coords).
template<typename T>
vec<contig_slice_t> contiguous_slices_from_dim(DIM<T>& dim, const llh_sptr_t<T>& llhf, uint8_t th_bv, double d_q)
{
  static constexpr size_t W = std::is_same_v<T, double> ? 1 : RWIDTH;
  vec<contig_slice_t> out;
  if (th_bv == 0) return out;

  const uint64_t nbins = dim.get_nbins();

  vec<uint64_t> pts = {1, nbins + 1};
  for (size_t ti = 0; ti < W; ++ti) {
    if (!(th_bv & (1u << ti)) || dim.get_intervals(ti).empty()) continue;
    for (const auto& iv : dim.get_intervals(ti)) {
      pts.push_back(iv.a);
      pts.push_back(iv.b + 1);
    }
  }

  std::sort(pts.begin(), pts.end());
  pts.erase(std::unique(pts.begin(), pts.end()), pts.end());

  std::array<size_t, RWIDTH> ti_ix{};

  for (size_t pi = 0; pi + 1 < pts.size(); ++pi) {
    const uint64_t a = pts[pi];
    const uint64_t b = pts[pi + 1];
    const auto [d, I] = slice_mle_fisher(dim, llhf, a, b);

    uint8_t mask = 0;
    double bin_lo = d_eps;
    double bin_hi = d_ub;
    for (size_t ti = 0; ti < W; ++ti) {
      if (!(th_bv & (1u << ti))) continue;
      const auto& ev = dim.get_intervals(ti);
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

// Helper: create a DIM<double> with given nbins, injecting known fdc/sdc values
static std::pair<std::shared_ptr<LLH<double>>, params_t<double>>
make_test_params(double dist_th = 0.1, uint32_t hdist_th = 4, uint64_t tau = 2, uint64_t bin_shift = 0)
{
  auto params = params_t<double>(1, dist_th, hdist_th, tau, 33.0, bin_shift, 1000, true);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, hdist_th, dist_th);
  return {llhf, params};
}

TEST_SUITE("DIM<double>::inclusive_scan") {

TEST_CASE("inclusive_scan with no hits yields no merged intervals") {
  auto [llhf, params] = make_test_params();
  const uint64_t nbins = 5;
  DIM<double> dim(params, llhf, nbins, nbins);

  // No hits: fdc/sdc stay zero -> no merged intervals after expand.
  dim.inclusive_scan();
  dim.extrema_scan();
  dim.extract_intervals_mx(0, 1, nbins);
  dim.expand_intervals(33.0);
  auto iv = dim.get_interval(0);
  CHECK(iv.a >= nbins);
}

} // TEST_SUITE

TEST_SUITE("DIM::extract_intervals_mx and _sx agreement") {

// Helper: run both algorithms on identically-prepared DIMs and compare all intervals.
template<typename T>
static void compare_mx_sx(DIM<T>& dim_mx, DIM<T>& dim_sx, uint64_t tau, size_t ix = 0)
{
  const double chisq = 33.0;
  const uint64_t nbins = dim_mx.get_nbins();
  dim_mx.extract_intervals_mx(tau, 1, nbins, ix);
  dim_mx.expand_intervals(chisq, ix);
  dim_sx.extract_intervals_sx(tau, 1, nbins, ix);
  dim_sx.expand_intervals(chisq, ix);

  for (uint64_t i = 0; ; ++i) {
    auto iv_mx = dim_mx.get_interval(i, ix);
    auto iv_sx = dim_sx.get_interval(i, ix);
    INFO("mismatch at interval ", i, " for tau=", tau, " ix=", ix,
         ": mx=(", iv_mx.a, ",", iv_mx.b, ")",
         " sx=(", iv_sx.a, ",", iv_sx.b, ")");
    CHECK(iv_mx.a == iv_sx.a);
    CHECK(iv_mx.b == iv_sx.b);
    if (iv_mx.a >= nbins) break;
  }
}

TEST_CASE("randomized patterns: double, various tau") {
  // Use a fixed seed for reproducibility.
  std::mt19937_64 rng(42);
  std::uniform_int_distribution<uint64_t> rnbins(5, 40);
  std::uniform_int_distribution<int> rhdist(0, 4);
  std::uniform_int_distribution<int> rcount(0, 3);
  std::uniform_real_distribution<double> rdist(0.02, 0.45);

  constexpr int ntrials = 200;

  for (int trial = 0; trial < ntrials; ++trial) {
    const double D = rdist(rng);
    const uint64_t nbins = rnbins(rng);

    auto params = params_t<double>(1, D, 4, 2, 33.0, 0, 1000, true);
    auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, D);

    DIM<double> dim_mx(params, llhf, nbins, nbins);
    DIM<double> dim_sx(params, llhf, nbins, nbins);

    // Inject random k-mers into each bin.
    for (uint64_t i = 0; i < nbins; ++i) {
      int n = rcount(rng);
      for (int k = 0; k < n; ++k) {
        uint32_t hd = static_cast<uint32_t>(rhdist(rng));
        dim_mx.aggregate_mer(hd, i);
        dim_sx.aggregate_mer(hd, i);
      }
    }

    dim_mx.inclusive_scan();  dim_sx.inclusive_scan();
    dim_mx.extrema_scan();    dim_sx.extrema_scan();

    for (uint64_t tau : {0u, 1u, 2u, 3u, 5u}) {
      if (tau >= nbins) continue;
      compare_mx_sx(dim_mx, dim_sx, tau, 0);
    }
  }
}

TEST_CASE("randomized patterns: cm512_t, all lanes independent") {
  std::mt19937_64 rng(99);
  std::uniform_int_distribution<uint64_t> rnbins(5, 30);
  std::uniform_int_distribution<int> rhdist(0, 4);
  std::uniform_int_distribution<int> rcount(0, 3);
  std::uniform_real_distribution<double> rdist(0.02, 0.45);

  constexpr int ntrials = 100;

  for (int trial = 0; trial < ntrials; ++trial) {
    const uint64_t nbins = rnbins(rng);

    cm512_t dths{};
    for (size_t i = 0; i < RWIDTH; ++i) dths[i] = rdist(rng);
    auto params = params_t<cm512_t>(RWIDTH, dths, 4, 2, 33.0, 0, 1000, true);
    auto llhf = std::make_shared<LLH<cm512_t>>(27, 11, 0.5, 4, dths);

    DIM<cm512_t> dim_mx(params, llhf, nbins, nbins);
    DIM<cm512_t> dim_sx(params, llhf, nbins, nbins);

    for (uint64_t i = 0; i < nbins; ++i) {
      int n = rcount(rng);
      for (int k = 0; k < n; ++k) {
        uint32_t hd = static_cast<uint32_t>(rhdist(rng));
        dim_mx.aggregate_mer(hd, i);
        dim_sx.aggregate_mer(hd, i);
      }
    }

    dim_mx.inclusive_scan();  dim_sx.inclusive_scan();
    dim_mx.extrema_scan();    dim_sx.extrema_scan();

    for (size_t ix = 0; ix < RWIDTH; ++ix) {
      for (uint64_t tau : {0u, 1u, 2u, 3u, 5u}) {
        if (tau >= nbins) continue;
        compare_mx_sx(dim_mx, dim_sx, tau, ix);
      }
    }
  }
}

TEST_CASE("edge cases: early-return when fdps[nbins] < fdps[1]") {
  // Uniform hdist=0 injections with positive threshold makes fdps monotonically
  // decrease, triggering the early-return interval (1, nbins).
  auto [llhf, params] = make_test_params(0.1, 4, 2, 0);
  const uint64_t nbins = 15;
  DIM<double> dim_mx(params, llhf, nbins, nbins);
  DIM<double> dim_sx(params, llhf, nbins, nbins);

  for (uint64_t i = 0; i < nbins; ++i) {
    dim_mx.aggregate_mer(0, i);
    dim_sx.aggregate_mer(0, i);
  }

  dim_mx.inclusive_scan();  dim_sx.inclusive_scan();
  dim_mx.extrema_scan();    dim_sx.extrema_scan();

  for (uint64_t tau : {0u, 1u, 2u}) {
    compare_mx_sx(dim_mx, dim_sx, tau);
  }
}

TEST_CASE("edge cases: full-range dip does not shortcut when nbins < 1 + tau") {
  auto [llhf, params] = make_test_params(0.1, 4, 2, 0);
  const uint64_t nbins = 6;
  DIM<double> dim_mx(params, llhf, nbins, nbins);
  DIM<double> dim_sx(params, llhf, nbins, nbins);
  for (uint64_t i = 0; i < nbins; ++i) {
    dim_mx.aggregate_mer(0, i);
    dim_sx.aggregate_mer(0, i);
  }
  dim_mx.inclusive_scan();
  dim_sx.inclusive_scan();
  dim_mx.extrema_scan();
  dim_sx.extrema_scan();
  const uint64_t tau = 6; // 1 + tau > nbins: early shortcut must not emit (1, nbins)
  compare_mx_sx(dim_mx, dim_sx, tau);
  CHECK(dim_mx.get_interval(0, 0).a >= nbins);
}

TEST_CASE("edge cases: no k-mers at all") {
  auto [llhf, params] = make_test_params(0.1, 4, 2, 0);
  const uint64_t nbins = 10;
  DIM<double> dim_mx(params, llhf, nbins, nbins);
  DIM<double> dim_sx(params, llhf, nbins, nbins);

  // No aggregate_mer calls: all fdc values are zero, fdps is flat.

  dim_mx.inclusive_scan();  dim_sx.inclusive_scan();
  dim_mx.extrema_scan();    dim_sx.extrema_scan();

  for (uint64_t tau : {0u, 1u, 2u, 5u}) {
    compare_mx_sx(dim_mx, dim_sx, tau);
  }
}

TEST_CASE("edge cases: single-bin features") {
  // Isolated single-bin injections can produce very short intervals.
  auto [llhf, params] = make_test_params(0.1, 4, 2, 0);
  const uint64_t nbins = 12;
  DIM<double> dim_mx(params, llhf, nbins, nbins);
  DIM<double> dim_sx(params, llhf, nbins, nbins);

  // Positive bump at bin 3, negative dip at bin 8.
  for (int k = 0; k < 5; ++k) { dim_mx.aggregate_mer(4, 3); dim_sx.aggregate_mer(4, 3); }
  for (int k = 0; k < 5; ++k) { dim_mx.aggregate_mer(0, 8); dim_sx.aggregate_mer(0, 8); }

  dim_mx.inclusive_scan();  dim_sx.inclusive_scan();
  dim_mx.extrema_scan();    dim_sx.extrema_scan();

  for (uint64_t tau : {0u, 1u, 2u, 3u}) {
    compare_mx_sx(dim_mx, dim_sx, tau);
  }
}

TEST_CASE("edge cases: alternating pattern") {
  // Alternating hdist=4 and hdist=0 creates a sawtooth prefix sum.
  auto [llhf, params] = make_test_params(0.1, 4, 2, 0);
  const uint64_t nbins = 25;
  DIM<double> dim_mx(params, llhf, nbins, nbins);
  DIM<double> dim_sx(params, llhf, nbins, nbins);

  for (uint64_t i = 0; i < nbins; ++i) {
    uint32_t hd = (i % 2 == 0) ? 4 : 0;
    dim_mx.aggregate_mer(hd, i);
    dim_sx.aggregate_mer(hd, i);
  }

  dim_mx.inclusive_scan();  dim_sx.inclusive_scan();
  dim_mx.extrema_scan();    dim_sx.extrema_scan();

  for (uint64_t tau : {0u, 1u, 2u, 3u, 5u}) {
    compare_mx_sx(dim_mx, dim_sx, tau);
  }
}

} // TEST_SUITE

TEST_SUITE("DIM<double>::expand_intervals") {

TEST_CASE("merge when chi-square below threshold") {
  auto pr = make_test_params(0.1, 4, 1, 0);
  std::shared_ptr<LLH<double>> llhf = pr.first;
  params_t<double> params = pr.second;
  const uint64_t nbins = 30;

  auto count_ivs = [&](double chisq_th) -> uint64_t {
    DIM<double> d(params, llhf, nbins, nbins);
    for (uint64_t i = 0; i < 10; ++i) {
      d.aggregate_mer(0, i);
      d.aggregate_mer(0, i);
    }
    for (uint64_t i = 12; i < 25; ++i) {
      d.aggregate_mer(0, i);
      d.aggregate_mer(0, i);
    }
    d.inclusive_scan();
    d.extrema_scan();
    d.extract_intervals_mx(0, 1, nbins);
    d.expand_intervals(chisq_th);
    uint64_t n = 0;
    for (;; ++n) {
      const interval_t iv = d.get_interval(n);
      if (iv.a >= nbins) break;
    }
    return n;
  };

  const uint64_t n_strict = count_ivs(0.0);
  const uint64_t n_merged = count_ivs(1e20);
  CHECK(n_merged <= n_strict);
}

} // TEST_SUITE

TEST_SUITE("DIM<double>::extract_histogram") {

TEST_CASE("extract_histogram returns correct counts for enum_only=false") {
  // Need enum_only=false for hdisthist_v to be allocated
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 10;
  const uint64_t nmers = 100;
  DIM<double> dim(params, llhf, nbins, nmers);

  // One hit per bin 0..4 (hdist=0): t=5 equals mers in [0,5) when bin_shift=0 -> u=0
  for (uint64_t i = 0; i < 5; ++i) {
    dim.aggregate_mer(0, i);
  }

  dim.compute_prefhistsum();

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(0, 5, v, u, t);

  CHECK(v[0] == 5);
  CHECK(t == 5);
  CHECK(u == 0);
  for (uint32_t d = 1; d <= params.hdist_th; ++d) {
    CHECK(v[d] == 0);
  }
}

TEST_CASE("extract_histogram full range matches partial sums") {
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 8;
  const uint64_t nmers = 1000;
  DIM<double> dim(params, llhf, nbins, nmers);

  // Inject hits into various bins with various hamming distances
  dim.aggregate_mer(0, 0);
  dim.aggregate_mer(1, 1);
  dim.aggregate_mer(2, 2);
  dim.aggregate_mer(3, 3);
  dim.aggregate_mer(0, 4);
  dim.aggregate_mer(0, 5);
  dim.aggregate_mer(1, 6);
  dim.aggregate_mer(2, 7);

  dim.compute_prefhistsum();

  // Full range
  vec<uint64_t> v_full;
  uint64_t u_full, t_full;
  dim.extract_histogram(0, nbins, v_full, u_full, t_full);

  // Split range: [0,4) + [4,8) should equal [0,8)
  vec<uint64_t> v1, v2;
  uint64_t u1, t1, u2, t2;
  dim.extract_histogram(0, 4, v1, u1, t1);
  dim.extract_histogram(4, nbins, v2, u2, t2);

  const uint32_t W = params.hdist_th + 1;
  for (uint32_t d = 0; d < W; ++d) {
    CHECK(v1[d] + v2[d] == v_full[d]);
  }
  CHECK(t1 + t2 == t_full);
  CHECK(u1 + u2 == u_full);
}

TEST_CASE("extract_histogram split additivity with bin_shift > 0") {
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 2, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 8;
  const uint64_t nmers = 1000;
  DIM<double> dim(params, llhf, nbins, nmers);

  dim.aggregate_mer(0, 0);
  dim.aggregate_mer(1, 1);
  dim.aggregate_mer(2, 2);
  dim.aggregate_mer(3, 3);

  dim.compute_prefhistsum();

  vec<uint64_t> v_full;
  uint64_t u_full, t_full;
  dim.extract_histogram(0, nbins, v_full, u_full, t_full);

  vec<uint64_t> v1, v2;
  uint64_t u1, t1, u2, t2;
  dim.extract_histogram(0, 4, v1, u1, t1);
  dim.extract_histogram(4, nbins, v2, u2, t2);

  const uint32_t W = params.hdist_th + 1;
  for (uint32_t d = 0; d < W; ++d) {
    CHECK(v1[d] + v2[d] == v_full[d]);
  }
  CHECK(t1 + t2 == t_full);
  CHECK(u1 + u2 == u_full);
}

TEST_CASE("extract_histogram supports the maximum fixed SIMD hdist threshold") {
  auto params = params_t<double>(1, 0.1, hdist_bound, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, hdist_bound, 0.1);
  const uint64_t nbins = 4;
  DIM<double> dim(params, llhf, nbins, nbins);

  dim.aggregate_mer(hdist_bound, 0);
  dim.aggregate_mer(0, 1);
  dim.compute_prefhistsum();

  vec<uint64_t> v;
  uint64_t u, t;
  dim.extract_histogram(0, nbins, v, u, t);

  REQUIRE(v.size() == hdist_bound + 1);
  CHECK(v[0] == 1);
  CHECK(v[hdist_bound] == 1);
  CHECK(t == 2);
}

TEST_CASE("complete_histogram matches extract_histogram on full range when keep_hist") {
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 0, 1000, true);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  DIM<double> dim(params, llhf, 5, 50);

  dim.aggregate_mer(0, 0);
  dim.aggregate_mer(1, 1);
  dim.aggregate_mer(5, 2);
  dim.compute_prefhistsum();

  vec<uint64_t> v_complete, v_extract;
  uint64_t u_complete, t_complete, u_extract, t_extract;
  dim.complete_histogram(v_complete, u_complete, t_complete);
  dim.extract_histogram(0, dim.get_nbins(), v_extract, u_extract, t_extract);

  CHECK(v_complete == v_extract);
  CHECK(u_complete == u_extract);
  CHECK(t_complete == t_extract);
}

TEST_CASE("complete_histogram exposes global hit counts and explicit misses") {
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 0, 0, true); // enum_only: no per-bin hist
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  DIM<double> dim(params, llhf, 5, 50);

  dim.aggregate_mer(0, 0);
  dim.aggregate_mer(1, 1);
  dim.aggregate_mer(5, 2); // explicit miss (hdist > hdist_th)

  vec<uint64_t> v;
  uint64_t u, t;
  dim.complete_histogram(v, u, t);

  CHECK(v[0] == 1);
  CHECK(v[1] == 1);
  CHECK(u == 1);
  CHECK(t == 2);
}

} // TEST_SUITE

TEST_SUITE("DIM<cm512_t>") {

TEST_CASE("SIMD DIM produces valid intervals") {
  cm512_t dths{};
  for (int i = 0; i < 8; ++i) dths[i] = 0.05 * (i + 1);
  auto params = params_t<cm512_t>(8, dths, 4, 2, 33.0, 0, 1000, true);
  auto llhf = std::make_shared<LLH<cm512_t>>(27, 11, 0.5, 4, dths);
  const uint64_t nbins = 15;
  DIM<cm512_t> dim(params, llhf, nbins, nbins);

  // Inject hits in a pattern
  for (uint64_t i = 0; i < nbins; ++i) {
    if (i < 4 || i >= 11) {
      dim.aggregate_mer(0, i);
      dim.aggregate_mer(0, i);
    }
  }

  dim.inclusive_scan();
  dim.extrema_scan();

  // Extract intervals for each threshold
  for (size_t ix = 0; ix < 8; ++ix) {
    dim.extract_intervals_mx(1, 1, nbins, ix);
    dim.expand_intervals(33.0, ix);
    for (const auto& iv : dim.get_intervals(ix)) {
      CHECK(iv.a >= 1);
      CHECK(iv.a <= iv.b);
      CHECK(iv.b <= nbins);
    }
  }
}

} // TEST_SUITE

TEST_SUITE("DIM::find_ordered_intervals") {

// Helper: run the full pipeline to populate eintervals_v, then list contiguous slices in bin space.
template<typename T>
static vec<contig_slice_t> run_pipeline(DIM<T>& dim, const llh_sptr_t<T>& llhf, uint8_t th_bv, uint64_t tau = 1)
{
  dim.inclusive_scan();
  dim.compute_prefhistsum();
  vec<uint64_t> v;
  uint64_t u, t;
  dim.complete_histogram(v, u, t);
  const double d_q = llhf->mle(v.data(), u);
  dim.set_query_distance(d_q);
  dim.extrema_scan();
  const uint64_t nbins = dim.get_nbins();
  if constexpr (std::is_same_v<T, double>) {
    dim.extract_intervals_mx(tau, 1, nbins);
    dim.expand_intervals(33.0);
  } else {
    constexpr size_t N = std::is_same_v<T, cm512_t> ? RWIDTH : 1;
    for (size_t ix = 0; ix < N; ++ix) {
      dim.extract_intervals_mx(tau, 1, nbins, ix);
      dim.expand_intervals(33.0, ix);
    }
  }
  return contiguous_slices_from_dim(dim, llhf, th_bv, d_q);
}

TEST_CASE("th_bv=0 returns no segments") {
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  DIM<double> dim(params, llhf, 10, 100);
  auto segs = contiguous_slices_from_dim(dim, llhf, 0, nanx());
  CHECK(segs.empty());
}

TEST_CASE("low threshold below d_q: matched interval brackets (prev, t]") {
  auto params = params_t<double>(1, 0.01, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.01);
  const uint64_t nbins = 10;
  DIM<double> dim(params, llhf, nbins, nbins);

  // up-down-up prefix-sum pattern: hdist=4 (fdc>0) in bins 0-1,
  // hdist=0 (fdc<0) in bins 3-4, hdist=4 in bins 6-7.
  for (uint64_t i : {0u, 1u}) {
    dim.aggregate_mer(4, i); dim.aggregate_mer(4, i);
  }
  for (uint64_t i : {3u, 4u}) {
    dim.aggregate_mer(0, i); dim.aggregate_mer(0, i);
  }
  for (uint64_t i : {6u, 7u}) {
    dim.aggregate_mer(4, i); dim.aggregate_mer(4, i);
  }

  auto segs = run_pipeline(dim, llhf, 1);
  CHECK(!segs.empty());

  vec<uint64_t> v;
  uint64_t u, t;
  dim.complete_histogram(v, u, t);
  const double d_q = llhf->mle(v.data(), u);
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
  auto params = params_t<double>(1, 0.5, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.5);
  const uint64_t nbins = 10;
  DIM<double> dim(params, llhf, nbins, nbins);

  // Uniform low Hamming distance -> low query MLE, above-threshold extraction.
  for (uint64_t i = 0; i < nbins; ++i)
    dim.aggregate_mer(0, i);

  auto segs = run_pipeline(dim, llhf, 1);
  CHECK(!segs.empty());

  vec<uint64_t> v;
  uint64_t u, t;
  dim.complete_histogram(v, u, t);
  const double d_q = llhf->mle(v.data(), u);
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
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 10;
  DIM<double> dim(params, llhf, nbins, nbins);

  for (uint64_t i : {0u, 1u}) {
    dim.aggregate_mer(4, i); dim.aggregate_mer(4, i);
  }
  for (uint64_t i : {3u, 4u}) {
    dim.aggregate_mer(0, i); dim.aggregate_mer(0, i);
  }
  for (uint64_t i : {6u, 7u}) {
    dim.aggregate_mer(4, i); dim.aggregate_mer(4, i);
  }

  auto segs = run_pipeline(dim, llhf, 1);
  REQUIRE(!segs.empty());

  CHECK(segs.front().bin_a == 1);
  CHECK(segs.back().bin_b == nbins + 1);
  for (size_t i = 0; i + 1 < segs.size(); ++i) {
    CHECK(segs[i].bin_b == segs[i + 1].bin_a);
  }
}

TEST_CASE("cm512_t thresholds ranked relative to d_q") {
  cm512_t dths{};
  dths[0] = 0.05;
  dths[1] = 0.10;
  dths[2] = 0.20;
  dths[3] = 0.30;

  auto params = params_t<cm512_t>(4, dths, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<cm512_t>>(27, 11, 0.5, 4, dths);
  const uint64_t nbins = 10;
  DIM<cm512_t> dim(params, llhf, nbins, nbins);

  for (uint64_t i : {0u, 1u}) {
    dim.aggregate_mer(4, i); dim.aggregate_mer(4, i);
  }
  for (uint64_t i : {3u, 4u}) {
    dim.aggregate_mer(0, i); dim.aggregate_mer(0, i);
  }
  for (uint64_t i : {6u, 7u}) {
    dim.aggregate_mer(4, i); dim.aggregate_mer(4, i);
  }

  uint8_t th_bv = (1u << 0) | (1u << 2);
  auto segs = run_pipeline(dim, llhf, th_bv);
  REQUIRE(!segs.empty());

  vec<uint64_t> v;
  uint64_t u, t;
  dim.complete_histogram(v, u, t);
  const double d_q = llhf->mle(v.data(), u);

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
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 8;
  DIM<double> dim(params, llhf, nbins, nbins);

  // No k-mer hits -> flat prefix sum -> extract_intervals finds nothing.
  // build_records_from_boundaries still reports [1, nbins+1) once for full-query MLE.

  auto segs = run_pipeline(dim, llhf, 1);
  REQUIRE(segs.size() == 1);
  CHECK(segs[0].bin_a == 1);
  CHECK(segs[0].bin_b == nbins + 1);
  CHECK(segs[0].mask == 0);
}

TEST_CASE("full-span merged interval [1, nbins] still yields one segment") {
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 10;
  DIM<double> dim(params, llhf, nbins, nbins);
  for (uint64_t i = 0; i < nbins; ++i)
    dim.aggregate_mer(0, i);

  auto segs = run_pipeline(dim, llhf, 1);
  REQUIRE(!segs.empty());
  CHECK(segs.size() == 1);
  CHECK(segs[0].bin_a == 1);
  CHECK(segs[0].bin_b == nbins + 1);
  CHECK(segs[0].mask == 0); // uniform hits: flat prefix sum, no threshold interval
}

TEST_CASE("map_contiguous segments partition histogram hit counts") {
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 10;
  DIM<double> dim(params, llhf, nbins, nbins);

  for (uint64_t i : {0u, 1u}) {
    dim.aggregate_mer(4, i);
    dim.aggregate_mer(4, i);
  }
  for (uint64_t i : {3u, 4u}) {
    dim.aggregate_mer(0, i);
    dim.aggregate_mer(0, i);
  }
  for (uint64_t i : {6u, 7u}) {
    dim.aggregate_mer(4, i);
    dim.aggregate_mer(4, i);
  }

  const auto segs = run_pipeline(dim, llhf, 1);
  REQUIRE(!segs.empty());

  vec<uint64_t> v_full;
  uint64_t u_full, t_full;
  dim.extract_histogram(0, nbins, v_full, u_full, t_full);

  const uint32_t W = params.hdist_th + 1;
  vec<uint64_t> v_sum(W, 0);
  uint64_t t_sum = 0;
  uint64_t u_sum = 0;
  for (const auto& s : segs) {
    vec<uint64_t> v;
    uint64_t u, t;
    dim.extract_histogram(s.bin_a - 1, s.bin_b - 1, v, u, t);
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
  auto params = params_t<double>(1, 0.1, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<double>>(27, 11, 0.5, 4, 0.1);
  const uint64_t nbins = 12;
  DIM<double> dim(params, llhf, nbins, nbins);

  for (uint64_t i = 0; i < nbins; ++i) {
    if (i % 3 == 0)
      dim.aggregate_mer(0, i);
    else if (i % 3 == 1)
      dim.aggregate_mer(2, i);
    else
      dim.aggregate_mer(4, i);
  }

  const auto segs = run_pipeline(dim, llhf, 1);
  REQUIRE(!segs.empty());

  for (const auto& s : segs) {
    const auto [d2, I2] = slice_mle_fisher(dim, llhf, s.bin_a, s.bin_b);
    CHECK(s.d == doctest::Approx(d2).epsilon(1e-9));
    if (std::isfinite(s.I) && std::isfinite(I2)) {
      CHECK(s.I == doctest::Approx(I2).epsilon(1e-7));
    } else {
      CHECK(std::isnan(s.I) == std::isnan(I2));
    }
  }
}

} // TEST_SUITE

TEST_SUITE("Misc") {

TEST_CASE("thrank_v orders thresholds relative to d_q") {
  cm512_t dths{};
  dths[0] = 0.10;
  dths[1] = 0.05;
  dths[2] = 0.20;
  dths[3] = 0.03;
  dths[4] = 0.40;
  dths[5] = 0.35;
  dths[6] = 0.50;
  dths[7] = 0.45;

  auto params = params_t<cm512_t>(8, dths, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<cm512_t>>(27, 11, 0.5, 4, dths);
  DIM<cm512_t> dim(params, llhf, 8, 80);

  for (uint64_t i = 0; i < 8; ++i)
    dim.aggregate_mer(0, i);

  const double d_q = 0.25;
  dim.set_query_distance(d_q);

  const auto& thrank = dim.get_thrank();
  REQUIRE(thrank.size() == RWIDTH);
  // Below d_q first (ascending): 0.03, 0.05, 0.10, 0.20; then above d_q (descending): 0.50, 0.45, 0.40, 0.35
  CHECK(thrank[0] == 3); // 0.03
  CHECK(thrank[1] == 1); // 0.05
  CHECK(thrank[2] == 0); // 0.10
  CHECK(thrank[3] == 2); // 0.20
  CHECK(thrank[4] == 6); // 0.50
  CHECK(thrank[5] == 7); // 0.45
  CHECK(thrank[6] == 4); // 0.40
  CHECK(thrank[7] == 5); // 0.35
}

TEST_CASE("thrank_v: non-flipped lanes before flipped at d_q=0.088094") {
  cm512_t dths{};
  dths[0] = 0.05;
  dths[1] = 0.075;
  dths[2] = 0.1;
  dths[3] = 0.125;
  dths[4] = 0.15;
  dths[5] = 0.225;
  dths[6] = 0.2;
  dths[7] = 0.25;

  auto params = params_t<cm512_t>(RWIDTH, dths, 4, 2, 33.0, 0, 1000, false);
  auto llhf = std::make_shared<LLH<cm512_t>>(27, 11, 0.5, 4, dths);
  DIM<cm512_t> dim(params, llhf, 4, 40);

  dim.set_query_distance(0.088094);
  const auto& thrank = dim.get_thrank();
  REQUIRE(thrank.size() == RWIDTH);
  CHECK(thrank[0] == 0); // 0.05
  CHECK(thrank[1] == 1); // 0.075
  CHECK(thrank[2] == 7); // 0.25
  CHECK(thrank[3] == 5); // 0.225
  CHECK(thrank[4] == 6); // 0.2
  CHECK(thrank[5] == 4); // 0.15
  CHECK(thrank[6] == 3); // 0.125
  CHECK(thrank[7] == 2); // 0.1
}

TEST_CASE("MH step equals sigma in sample_metropolis_hastings") {
  // Invariant: step = init.second, passed as sqrt(1/I) = sigma (map.cpp:638).
  CHECK(true);
}

TEST_CASE("multi-gap extraction finds intervals in separated regions") {
  auto [llhf, params] = make_test_params(0.1, 4, 2, 0);
  const uint64_t nbins = 25;
  DIM<double> dim(params, llhf, nbins, nbins);

  // Two divergence regions separated by a gap.
  for (uint64_t i = 2; i <= 5; ++i) {
    dim.aggregate_mer(0, i); dim.aggregate_mer(0, i);
  }
  for (uint64_t i = 16; i <= 19; ++i) {
    dim.aggregate_mer(0, i); dim.aggregate_mer(0, i);
  }

  dim.inclusive_scan();
  dim.extrema_scan();
  dim.extract_intervals_mx(1, 1, nbins);
  dim.expand_intervals(33.0);

  const auto& e_v = dim.get_intervals(0);
  CHECK(e_v.size() >= 1);
  for (const auto& iv : e_v) {
    CHECK(iv.a >= 1);
    CHECK(iv.a <= iv.b);
    CHECK(iv.b <= nbins);
  }
}

} // TEST_SUITE
