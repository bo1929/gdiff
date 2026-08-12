#include "dim.hpp"
#include "random.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <limits>
#include <numeric>
#include <random>

namespace {
  struct pv_t
  {
    uint64_t pos;
    double val;
  };
} // namespace

HDHist::HDHist(const uint64_t nbins, const uint32_t hdist_th, const uint64_t bin_shift)
  : nbins(nbins)
  , hdist_th(hdist_th)
  , bin_shift(bin_shift)
  , hist_v((nbins + 1) * (hdist_th + 1), 0)
  , miss_v(nbins + 1, 0)
{
  assert(this->bin_shift <= 16);
}

template<bool Atomic>
void HDHist::aggregate_mer(const uint32_t hdist_min, const uint64_t i)
{
  if (i >= nbins) return;
  if (hdist_min <= hdist_th) {
    if constexpr (Atomic)
      __atomic_add_fetch(&hist_v[((i + 1) * (hdist_th + 1)) + hdist_min], 1, __ATOMIC_RELAXED);
    else
      ++hist_v[((i + 1) * (hdist_th + 1)) + hdist_min];
  } else {
    if constexpr (Atomic)
      __atomic_add_fetch(&miss_v[i + 1], 1, __ATOMIC_RELAXED);
    else
      ++miss_v[i + 1];
  }
}

template void HDHist::aggregate_mer<false>(uint32_t, uint64_t);
template void HDHist::aggregate_mer<true>(uint32_t, uint64_t);

void HDHist::compute_prefhistsum()
{
  const uint32_t W = hdist_th + 1;
  for (uint64_t i = 0; i < nbins; ++i) {
    for (uint32_t d = 0; d < W; ++d) {
      hist_v[((i + 1) * W) + d] += hist_v[(i * W) + d];
    }
    miss_v[i + 1] += miss_v[i];
  }
}

void HDHist::compute_prefhistsum_parallel(ThreadPool& pool, uint32_t nchunks)
{
  const uint32_t W = hdist_th + 1;
  const uint64_t nrows = nbins + 1;
  if (nchunks <= 1 || nrows < (uint64_t(1) << 16)) {
    compute_prefhistsum();
    return;
  }
  nchunks = std::min<uint32_t>(nchunks, static_cast<uint32_t>((nrows + 4095) / 4096));
  if (nchunks <= 1) {
    compute_prefhistsum();
    return;
  }
  const uint64_t rows_per = (nrows + nchunks - 1) / nchunks;
  uint64_t* h = hist_v.data();
  // Phase A: exclusive-local prefix sums inside each chunk.
  pool.parallel_for(nchunks, 1, [&](uint64_t c) {
    const uint64_t c0 = c * rows_per;
    const uint64_t c1 = std::min(c0 + rows_per, nrows);
    for (uint64_t r = c0 + 1; r < c1; ++r) {
      for (uint32_t d = 0; d < W; ++d)
        h[r * W + d] += h[(r - 1) * W + d];
    }
  });
  // Serial combine of chunk bases (few chunks).
  vec<uint64_t> base(static_cast<uint64_t>(nchunks) * W, 0);
  for (uint32_t c = 1; c < nchunks; ++c) {
    const uint64_t last = std::min((c * rows_per), nrows) - 1;
    for (uint32_t d = 0; d < W; ++d)
      base[c * W + d] = base[(c - 1) * W + d] + h[last * W + d];
  }
  // Phase B: add each chunk's base to all of its rows.
  pool.parallel_for(nchunks - 1, 1, [&](uint64_t ci) {
    const uint64_t c = ci + 1;
    const uint64_t c0 = c * rows_per;
    const uint64_t c1 = std::min(c0 + rows_per, nrows);
    const uint64_t* b = &base[c * W];
    for (uint64_t r = c0; r < c1; ++r) {
      for (uint32_t d = 0; d < W; ++d)
        h[r * W + d] += b[d];
    }
  });
  for (uint64_t i = 0; i < nbins; ++i)
    miss_v[i + 1] += miss_v[i];
}

void HDHist::extract_histogram(const uint64_t a, const uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const
{
  assert(a <= b && b <= nbins);
  const uint32_t W = hdist_th + 1;
  assert(W <= RWIDTH && W <= hdist_bound + 1);
  // Copy into 8-lane scratch first. A direct mm512_maskz_loadu on &hist_v[row*W]
  // is unsafe when W < 8: SIMDe's non-native path (and some native masked loads
  // near a page boundary) still touch a full 64-byte vector and can SEGV on the
  // last prefix-sum rows.
  alignas(64) uint64_t hb[RWIDTH] = {};
  alignas(64) uint64_t ha[RWIDTH] = {};
  std::memcpy(hb, &hist_v[b * W], W * sizeof(uint64_t));
  std::memcpy(ha, &hist_v[a * W], W * sizeof(uint64_t));

  v.resize(hdist_bound + 1);
  const simde__mmask8 mask = static_cast<simde__mmask8>((1u << W) - 1);
  const simde__m512i vb = simde_mm512_maskz_loadu_epi64(mask, hb);
  const simde__m512i va = simde_mm512_maskz_loadu_epi64(mask, ha);
  const simde__m512i vd = simde_mm512_sub_epi64(vb, va);
  simde_mm512_storeu_si512(v.data(), vd);
  const simde__m256i lend = simde_mm512_castsi512_si256(vd);
  const simde__m256i rend = simde_mm512_extracti64x4_epi64(vd, 1);
  const simde__m256i s4 = simde_mm256_add_epi64(lend, rend);
  const simde__m128i s4_lend = simde_mm256_castsi256_si128(s4);
  const simde__m128i s4_rend = simde_mm256_extracti128_si256(s4, 1);
  const simde__m128i s2 = simde_mm_add_epi64(s4_lend, s4_rend);
  t = static_cast<uint64_t>(simde_mm_extract_epi64(s2, 0) + simde_mm_extract_epi64(s2, 1));
  u = miss_v[b] - miss_v[a];
}

void HDHist::extract_histogram(uint64_t a, uint64_t b, window_counts_t& wc) const
{
  uint64_t t = 0;
  extract_histogram(a, b, wc.hist_v, wc.u, t);
}

void HDHist::extract_histogram(uint64_t a, uint64_t b, swindow_counts_t& wc, bool is_rc) const
{
  uint64_t t = 0;
  if (is_rc)
    extract_histogram(a, b, wc.hist_rc_v, wc.u_rc, t);
  else
    extract_histogram(a, b, wc.hist_fw_v, wc.u_fw, t);
}

template<typename T>
DIM<T>::DIM(const params_t<T>& params, const llh_sptr_t<T>& llhf, uint64_t nbins, uint64_t nmers)
  : params(params)
  , llhf(llhf)
  , nbins(nbins)
  , nmers(nmers)
  , keep_hist(!params.enum_only || params.sample_size > 0)
{
  fdc_v.resize(nbins);
  sdc_v.resize(nbins);
  thneg_v.fill(false);
  thrank_v.resize(WIDTH);
  for (size_t i = 0; i < WIDTH; ++i)
    thrank_v[i] = i;
  if (keep_hist) {
    hdhist = HDHist(nbins, params.hdist_th, params.bin_shift);
  } else {
    hist_v.assign(hdist_bound + 1, 0);
  }
}

template<typename T>
void DIM<T>::set_query_distance(const double d_q)
{
  const bool is_valid = is_valid_distance(d_q);

  if constexpr (std::is_same_v<T, double>) {
    const double t = params.dist_th;
    thneg_v.front() = is_valid && t > d_q;
    thrank_v = {0};
  } else {
    arr<std::pair<double, size_t>, WIDTH> tp;
    for (size_t i = 0; i < WIDTH; ++i) {
      const double t = at(params.dist_th, i);
      thneg_v[i] = is_valid && t > d_q;
      tp[i] = {t, i};
    }

    if (!is_valid) {
      std::sort(tp.begin(), tp.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    } else {
      std::sort(tp.begin(), tp.end(), [&](const auto& a, const auto& b) {
        const double sa = a.first <= d_q ? a.first : (2.0 * d_ub - a.first);
        const double sb = b.first <= d_q ? b.first : (2.0 * d_ub - b.first);
        return sa < sb;
      });
    }

    thrank_v.resize(WIDTH);
    for (size_t i = 0; i < WIDTH; ++i)
      thrank_v[i] = tp[i].second;
  }

  apply_threshold_signs();
}

template<typename T>
void DIM<T>::apply_threshold_signs()
{
  // Flip the sign of fdps_v lanes whose threshold is on the opposite side of d_q.
  if constexpr (std::is_same_v<T, double>) {
    if (!thneg_v.front()) return;
    for (uint64_t i = 1; i < fdps_v.size(); ++i) {
      fdps_v[i] = -fdps_v[i];
    }
  } else {
    alignas(64) double xor_mask[WIDTH];
    bool flip_sign = false;
    for (size_t i = 0; i < WIDTH; ++i) {
      xor_mask[i] = thneg_v[i] ? -0.0 : 0.0;
      if (thneg_v[i]) {
        flip_sign = true;
      }
    }
    if (!flip_sign) return;
    const simde__m512d s = simde_mm512_loadu_pd(xor_mask);
    for (uint64_t i = 1; i < fdps_v.size(); ++i) {
      simde__m512d v = simde_mm512_loadu_pd(fdps_v[i].data());
      v = simde_mm512_xor_pd(v, s);
      simde_mm512_storeu_pd(fdps_v[i].data(), v);
    }
  }
}

template<typename T>
void DIM<T>::aggregate_mer(uint32_t hdist_min, uint64_t i)
{
  if (hdist_min <= params.hdist_th) {
    t_q++;
    if (keep_hist)
      hdhist.aggregate_mer(hdist_min, i);
    else
      hist_v[hdist_min]++;
    add_to(sdc_v[i], llhf->get_sdc(hdist_min));
    add_to(fdc_v[i], llhf->get_fdc(hdist_min));
  } else {
    u_q++;
    if (keep_hist) hdhist.aggregate_mer(hdist_min, i);
    add_to(sdc_v[i], llhf->get_sdc());
    add_to(fdc_v[i], llhf->get_fdc());
  }
}

template<typename T>
void DIM<T>::skip_mer(const uint64_t i)
{
  if (i >= nbins) return;
  if (skip_v.empty()) skip_v.assign(nbins, 0);
  skip_v[i] = 1;
  has_skips = true;
}

template<typename T>
void DIM<T>::inclusive_scan()
{
  assert(nbins > 0);
  const uint64_t s = nbins + 1;

  fdps_v.resize(s);
  sdps_v.resize(s);

  if constexpr (std::is_same_v<T, double>) {
    fdps_v[0] = 0.0;
    sdps_v[0] = 0.0;
    for (uint64_t i = 1; i < s; ++i) {
      fdps_v[i] = fdps_v[i - 1] + fdc_v[i - 1];
      sdps_v[i] = sdps_v[i - 1] + sdc_v[i - 1];
    }
  } else {
    fdps_v[0].fill(0.0);
    sdps_v[0].fill(0.0);
    simde__m512d fdps_acc = simde_mm512_setzero_pd();
    simde__m512d sdps_acc = simde_mm512_setzero_pd();
    for (uint64_t i = 1; i < s; ++i) {
      const simde__m512d fdc = simde_mm512_loadu_pd(fdc_v[i - 1].data());
      const simde__m512d sdc = simde_mm512_loadu_pd(sdc_v[i - 1].data());
      fdps_acc = simde_mm512_add_pd(fdps_acc, fdc);
      sdps_acc = simde_mm512_add_pd(sdps_acc, sdc);
      simde_mm512_storeu_pd(fdps_v[i].data(), fdps_acc);
      simde_mm512_storeu_pd(sdps_v[i].data(), sdps_acc);
    }
  }
}

template<typename T>
void DIM<T>::extrema_scan()
{
  const uint64_t s = nbins + 1;
  fdpmax_v.resize(s + 1);
  fdsmin_v.resize(s + 1);

  if constexpr (std::is_same_v<T, double>) {
    fdpmax_v[0] = ninf();
    fdsmin_v[0] = pinf();
    std::inclusive_scan(
      fdps_v.begin() + 1, fdps_v.end(), fdpmax_v.begin() + 1, [](double a, double b) { return std::max(a, b); });
    std::inclusive_scan(
      fdps_v.rbegin(), fdps_v.rend() - 1, fdsmin_v.rbegin() + 1, [](double a, double b) { return std::min(a, b); });
    fdpmax_v[s] = pinf();
    fdsmin_v[s] = ninf();
  } else {
    fdpmax_v[0].fill(ninf());
    fdsmin_v[0].fill(pinf());
    simde__m512d fdpmax_acc = simde_mm512_loadu_pd(fdpmax_v[0].data());
    simde__m512d fdsmin_acc = simde_mm512_loadu_pd(fdsmin_v[0].data());
    for (uint64_t i = 1; i < s; ++i) {
      const simde__m512d fdps_front = simde_mm512_loadu_pd(fdps_v[i].data());
      const simde__m512d fdps_back = simde_mm512_loadu_pd(fdps_v[s - i].data());
      fdpmax_acc = simde_mm512_max_pd(fdpmax_acc, fdps_front);
      fdsmin_acc = simde_mm512_min_pd(fdsmin_acc, fdps_back);
      simde_mm512_storeu_pd(fdpmax_v[i].data(), fdpmax_acc);
      simde_mm512_storeu_pd(fdsmin_v[s - i].data(), fdsmin_acc);
    }
    fdpmax_v[s].fill(pinf());
    fdsmin_v[s].fill(ninf());
  }
}

// Find maximal intervals [a, b] within [lix, rix] where the prefix sum drops below a prior maximum.
//
// Conditions for a valid interval (a, b):
//   1. a is a strict prefix maximum:  fdps[a] > fdpmax[a-1]
//   2. b is right-maximal:            fdsmin[b] < fdps[a] but fdsmin[b+1] >= fdps[a]
//   3. Minimum length:                b >= a + tau
//   4. Negative sum:                  fdps[b] < fdps[a]
//   5. Left-maximal (non-redundant):  fdps[b] >= fdpmax[a-1]
//   6. b not claimed by earlier a:    b != b_prev
//
// Early return: if prefix sum ends below where it started (fdps[rix] < fdps[lix]), the entire [lix, rix] is one interval.
template<typename T>
void DIM<T>::extract_intervals_mx(const uint64_t tau, const uint64_t lix, const uint64_t rix, const size_t ix)
{
  if (!has_skips) {
    extract_mx(tau, lix, rix, ix);
    return;
  }
  // Split [lix, rix] into maximal skip-free subranges (1-based bins, skip_v is 0-based).
  uint64_t a = lix;
  for (uint64_t x = lix; x <= rix; ++x) {
    if (skip_v[x - 1]) {
      if (x > a) extract_mx(tau, a, x - 1, ix);
      a = x + 1;
    }
  }
  if (rix >= a) extract_mx(tau, a, rix, ix);
}

template<typename T>
void DIM<T>::extract_mx(const uint64_t tau, const uint64_t lix, const uint64_t rix, const size_t ix)
{
  uint64_t b_curr = lix;
  uint64_t b_prev = std::numeric_limits<uint64_t>::max();

  if (rix >= lix + tau && at(fdps_v[rix], ix) < at(fdps_v[lix], ix)) {
    intervals_v[ix].emplace_back(lix, rix);
    return;
  }

  for (uint64_t a = lix; a <= rix; ++a) {
    const double fdpmax_a = at(fdpmax_v[a - 1], ix);
    const double fdps_a = at(fdps_v[a], ix);

    if (fdpmax_a >= fdps_a) continue; // Condition 1

    while ((b_curr + 1) <= rix && (at(fdsmin_v[b_curr + 1], ix) < fdps_a))
      ++b_curr; // Condition 2

    const uint64_t b_star = b_curr;
    if (b_star < (a + tau)) continue;               // Condition 3
    if (at(fdps_v[b_star], ix) >= fdps_a) continue; // Condition 4
    if (b_star == b_prev) continue;                 // Condition 6

    if (at(fdps_v[b_star], ix) >= fdpmax_a) { // Condition 5
      intervals_v[ix].emplace_back(a, b_star);
      b_prev = b_star;
    }
  }
}

template<typename T>
void DIM<T>::extract_intervals_sx(const uint64_t tau, const uint64_t lix, const uint64_t rix, const size_t ix)
{
  if (!has_skips) {
    extract_sx(tau, lix, rix, ix);
    return;
  }
  // Split [lix, rix] into maximal skip-free subranges (1-based bins, skip_v is 0-based).
  uint64_t a = lix;
  for (uint64_t x = lix; x <= rix; ++x) {
    if (skip_v[x - 1]) {
      if (x > a) extract_sx(tau, a, x - 1, ix);
      a = x + 1;
    }
  }
  if (rix >= a) extract_sx(tau, a, rix, ix);
}

template<typename T>
void DIM<T>::extract_sx(const uint64_t tau, const uint64_t lix, const uint64_t rix, const size_t ix)
{
  // Every valid right endpoint b* is a suffix minimum of fdps_v
  // Suffix minimum values are strictly increasing left-to-right,
  // Hence, the pointer into the list is monotone across record highs which is O(k) total
  const uint64_t gap_len = rix - lix + 1;
  if (gap_len >= 1 + tau && at(fdps_v[rix], ix) < at(fdps_v[lix], ix)) {
    intervals_v[ix].emplace_back(lix, rix);
    return;
  }

  vec<pv_t> pv_v;
  {
    double y_min = pinf();
    for (uint64_t j = rix; j >= lix; --j) {
      const double v = at(fdps_v[j], ix);
      if (v < y_min) {
        y_min = v;
        pv_v.push_back({j, v});
      }
    }
    std::reverse(pv_v.begin(), pv_v.end());
  }
  if (pv_v.empty()) return;

  size_t yix_min = 0;
  double running_max = ninf();
  uint64_t b_prev = std::numeric_limits<uint64_t>::max();

  for (uint64_t a = lix; a <= rix; ++a) {
    const double fdps_a = at(fdps_v[a], ix);
    const double fdpmax_a = running_max;
    if (fdps_a > running_max) running_max = fdps_a;

    if (fdpmax_a >= fdps_a) continue; // Skip if a is not a record high

    // Advance yix_min to the last suffix minimum with val < fdps_a
    while (yix_min + 1 < pv_v.size() && pv_v[yix_min + 1].val < fdps_a) {
      ++yix_min;
    }
    // Skip if no valid right endpoint with val < fdps_a
    if (pv_v[yix_min].val >= fdps_a) continue;

    const uint64_t b_star = pv_v[yix_min].pos;
    const double fdps_bstar = pv_v[yix_min].val;

    if (b_star < a + tau) continue; // Skip if no valid right endpoint in [a+tau, nbins]
    if (b_star == b_prev) continue; // Skip if b* was already claimed

    if (fdps_bstar >= fdpmax_a) {              // Left maximal
      intervals_v[ix].emplace_back(a, b_star); // 1-based inclusive coordinates
      b_prev = b_star;
    }
  }
}

template<typename T>
void DIM<T>::expand_intervals(const double chisq_th, const size_t ix)
{
  auto& iv_ix = intervals_v[ix];
  if (iv_ix.empty()) return;

  double fdiff, sdiff, chisq_val;
  uint64_t a, ap, b, bp;
  size_t w = 0;
  ap = iv_ix[0].a;
  bp = iv_ix[0].b;

  for (size_t i = 1; i < iv_ix.size(); ++i) {
    a = iv_ix[i].a;
    b = iv_ix[i].b;
    fdiff = at(fdps_v[b], ix) - at(fdps_v[ap], ix);
    sdiff = at(sdps_v[ap], ix) - at(sdps_v[b], ix);
    chisq_val = ((fdiff * fdiff) + eps) / (sdiff + eps);

    // Never merge across an N-run break: a skip bin in the gap (bp, a) keeps them apart.
    bool skip_gap = false;
    if (has_skips) {
      for (uint64_t x = bp + 1; x < a; ++x) {
        if (skip_v[x - 1]) {
          skip_gap = true;
          break;
        }
      }
    }

    if (!skip_gap && (chisq_val < chisq_th) && (a < bp)) {
      a = ap;
    } else {
      iv_ix[w++] = {ap, bp}; // 1-based inclusive coordinates
    }

    ap = a;
    bp = b;
  }

  iv_ix[w++] = {ap, bp}; // 1-based inclusive coordinates
  iv_ix.resize(w);
}

template<typename T>
void DIM<T>::compute_prefhistsum()
{
  if (!keep_hist) return;
  hdhist.compute_prefhistsum();
}

template<typename T>
void DIM<T>::total_histogram(vec<uint64_t>& v, uint64_t& u, uint64_t& t) const
{
  if (keep_hist) {
    hdhist.extract_histogram(0, nbins, v, u, t);
  } else {
    v = hist_v;
    u = u_q;
    t = t_q;
  }
}

template<typename T>
void DIM<T>::extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const
{
  hdhist.extract_histogram(a, b, v, u, t);
}

template<typename T>
void DIM<T>::extract_histogram(uint64_t a, uint64_t b, window_counts_t& wc) const
{
  uint64_t t = 0;
  extract_histogram(a, b, wc.hist_v, wc.u, t);
}

template<typename T>
void DIM<T>::extract_histogram(uint64_t a, uint64_t b, swindow_counts_t& wc, bool is_rc) const
{
  uint64_t t = 0;
  if (is_rc)
    extract_histogram(a, b, wc.hist_rc_v, wc.u_rc, t);
  else
    extract_histogram(a, b, wc.hist_fw_v, wc.u_fw, t);
}

template<typename T>
vec<sample_t> DIM<T>::sample_random_intervals(const uint64_t nwin_bins, const uint64_t bix) const
{
  vec<sample_t> out_v;
  if (nwin_bins == 0 || params.sample_size == 0 || nwin_bins > nbins) return out_v;

  const uint64_t npos = nbins - nwin_bins + 1;
  const size_t goal = static_cast<size_t>(std::min(params.sample_size, npos));
  out_v.reserve(goal);
  window_counts_t wc(llhf->hdist_th);

  vec<uint64_t> starts_v(npos);
  std::iota(starts_v.begin(), starts_v.end(), uint64_t(0));
  std::shuffle(starts_v.begin(), starts_v.end(), gen);

  for (const uint64_t x : starts_v) {
    if (out_v.size() >= goal) break;
    if (has_skips) {
      const auto first = skip_v.begin() + static_cast<ptrdiff_t>(x);
      const auto last = first + static_cast<ptrdiff_t>(nwin_bins);
      if (std::find(first, last, uint8_t(1)) != last) continue;
    }

    const uint64_t a_bin = x + 1;
    const uint64_t b_bin = x + nwin_bins + 1;
    wc.clear();
    extract_histogram(a_bin - 1, b_bin - 1, wc);
    const double d = llhf->mle(wc.hist(), wc.u);
    if (!is_valid_distance(d)) continue;
    const double I = llhf->compute_fisher_info(wc.hist(), wc.u, d);
    if (!std::isfinite(I) || I <= 0.0) continue;
    out_v.push_back({d, I, bix, {a_bin, b_bin}});
  }
  return out_v;
}

bool filter_background_samples(const vec<sample_t>& in_v, const record_t& r, uint64_t sample_size, vec<sample_t>& out_v)
{
  bool excluded = false;
  out_v.clear();
  const size_t goal = std::min(in_v.size(), static_cast<size_t>(sample_size));
  out_v.reserve(goal);
  size_t neligible = 0;

  for (const auto& s : in_v) {
    if (s.bin_iv.b - s.bin_iv.a != r.nbins) continue;
    if (s.bix == r.bix && overlaps_half_open(s.bin_iv, r.bin_iv)) {
      excluded = true;
      continue;
    }
    if (!std::isfinite(s.I) || s.I <= 0.0) continue;
    if (!is_valid_distance(s.d)) continue;

    ++neligible;
    if (out_v.size() < goal) {
      out_v.push_back(s);
    } else if (goal > 0) {
      const size_t j = std::uniform_int_distribution<size_t>(0, neligible - 1)(gen);
      if (j < goal) out_v[j] = s;
    }
  }
  return excluded;
}

template class DIM<double>;
template class DIM<cm512_t>;
