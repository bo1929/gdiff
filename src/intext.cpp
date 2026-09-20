#include "intext.hpp"
#include "scan.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <limits>
#include <numeric>

namespace {
  // Split [lix, rix] into maximal skip-free subranges (1-based bins, skip_v is 0-based).
  template<typename Fn>
  void for_each_contiguous(const vec<uint8_t>& skip_v, uint64_t lix, uint64_t rix, Fn&& fn)
  {
    uint64_t a = lix;
    for (uint64_t x = lix; x <= rix; ++x) {
      if (skip_v[x - 1]) {
        if (x > a) fn(a, x - 1);
        a = x + 1;
      }
    }
    if (rix >= a) fn(a, rix);
  }
} // namespace

Histogram::Histogram(uint64_t nbins, uint32_t hdist_th, uint64_t bin_shift)
  : nbins(nbins)
  , hdist_th(hdist_th)
  , hist_v((nbins + 1) * (hdist_th + 1), 0)
  , miss_v(nbins + 1, 0)
{
  assert(bin_shift <= 16);
}

void Histogram::aggregate_mer(uint32_t hdist_min, uint64_t i)
{
  if (i >= nbins) return;
  if (hdist_min <= hdist_th) {
    ++hist_v[((i + 1) * (hdist_th + 1)) + hdist_min];
  } else {
    ++miss_v[i + 1];
  }
}

void Histogram::compute_prefhistsum()
{
  const uint32_t W = hdist_th + 1;
  for (uint64_t i = 0; i < nbins; ++i) {
    for (uint32_t hd = 0; hd < W; ++hd) {
      hist_v[((i + 1) * W) + hd] += hist_v[(i * W) + hd];
    }
    miss_v[i + 1] += miss_v[i];
  }
}

void Histogram::extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const
{
  assert(a <= b && b <= nbins);
  const uint32_t W = hdist_th + 1;
  assert(W <= rwidth && W <= hdist_bound + 1);
  // Copy into 8-lane scratch; a masked load of W < 8 can fault near a page end.
  alignas(64) uint64_t hb[rwidth] = {};
  alignas(64) uint64_t ha[rwidth] = {};
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

void Histogram::extract_histogram(uint64_t a, uint64_t b, window_counts_t& wc) const
{
  uint64_t t = 0;
  extract_histogram(a, b, wc.hist_v, wc.u, t);
}

void Histogram::extract_histogram(uint64_t a, uint64_t b, swindow_counts_t& wc, bool is_rc) const
{
  uint64_t t = 0;
  if (is_rc)
    extract_histogram(a, b, wc.hist_rc_v, wc.u_rc, t);
  else
    extract_histogram(a, b, wc.hist_fw_v, wc.u_fw, t);
}

template<typename T>
IntExt<T>::IntExt(const map_params<T>& params, const LLH<T>& llhf, uint64_t nbins, uint64_t nmers)
  : params(params)
  , llhf(llhf)
  , nbins(nbins)
  , nmers(nmers)
  , hist(nbins, params.hdist_th, params.bin_shift)
{
  // Extra slot so inclusive_scan can prefix-sum in place. Aggregation uses [0, nbins).
  fdc_v.resize(nbins + 1);
  sdc_v.resize(nbins + 1);
  thneg_v.fill(false);
  thrank_v.resize(WIDTH);
  for (size_t i = 0; i < WIDTH; ++i)
    thrank_v[i] = i;
}

template<typename T>
void IntExt<T>::set_query_distance(double d_q)
{
  const bool is_valid = is_valid_distance(d_q);

  if constexpr (std::is_same_v<T, double>) {
    const double t = params.dist_th;
    thneg_v.front() = is_valid && t > d_q;
    thrank_v = {0};
  } else {
    arr<std::pair<double, size_t>, WIDTH> tp;
    for (size_t i = 0; i < WIDTH; ++i) {
      const double t = lane_at(params.dist_th, i);
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
void IntExt<T>::apply_threshold_signs()
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
void IntExt<T>::aggregate_mer(uint32_t hdist_min, uint64_t i)
{
  if (hdist_min <= params.hdist_th) {
    hist.aggregate_mer(hdist_min, i);
    add_to(sdc_v[i], llhf.get_sdc(hdist_min));
    add_to(fdc_v[i], llhf.get_fdc(hdist_min));
  } else {
    hist.aggregate_mer(hdist_min, i);
    add_to(sdc_v[i], llhf.get_sdc());
    add_to(fdc_v[i], llhf.get_fdc());
  }
}

template<typename T>
void IntExt<T>::skip_mer(uint64_t i)
{
  if (i >= nbins) return;
  if (skip_v.empty()) skip_v.assign(nbins, 0);
  skip_v[i] = 1;
  wskip = true;
}

// In-place prefix sum: C[0]=0, C[i]=C[i-1]+c[i-1].
template<typename T>
void IntExt<T>::inclusive_scan()
{
  assert(nbins > 0);
  const uint64_t s = nbins + 1;
  assert(fdc_v.size() == s && sdc_v.size() == s);

  std::move_backward(fdc_v.begin(), fdc_v.end() - 1, fdc_v.end());
  std::move_backward(sdc_v.begin(), sdc_v.end() - 1, sdc_v.end());

  if constexpr (std::is_same_v<T, double>) {
    fdc_v[0] = 0.0;
    sdc_v[0] = 0.0;
    for (uint64_t i = 1; i < s; ++i) {
      fdc_v[i] += fdc_v[i - 1];
      sdc_v[i] += sdc_v[i - 1];
    }
  } else {
    fdc_v[0].fill(0.0);
    sdc_v[0].fill(0.0);
    simde__m512d fdps_acc = simde_mm512_setzero_pd();
    simde__m512d sdps_acc = simde_mm512_setzero_pd();
    for (uint64_t i = 1; i < s; ++i) {
      fdps_acc = simde_mm512_add_pd(fdps_acc, simde_mm512_loadu_pd(fdc_v[i].data()));
      sdps_acc = simde_mm512_add_pd(sdps_acc, simde_mm512_loadu_pd(sdc_v[i].data()));
      simde_mm512_storeu_pd(fdc_v[i].data(), fdps_acc);
      simde_mm512_storeu_pd(sdc_v[i].data(), sdps_acc);
    }
  }
  fdps_v = std::move(fdc_v);
  sdps_v = std::move(sdc_v);
  vec<T>().swap(fdc_v);
  vec<T>().swap(sdc_v);
}

template<typename T>
void IntExt<T>::extrema_scan()
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

// Maximal [a,b] in [lix,rix] where the prefix sum drops below a prior maximum.
template<typename T>
void IntExt<T>::extract_intervals_mx(uint64_t tau, uint64_t lix, uint64_t rix, size_t ix)
{
  if (!wskip) {
    extract_mx(tau, lix, rix, ix);
    return;
  }
  for_each_contiguous(skip_v, lix, rix, [&](uint64_t a, uint64_t b) { extract_mx(tau, a, b, ix); });
}

template<typename T>
void IntExt<T>::extract_mx(uint64_t tau, uint64_t lix, uint64_t rix, size_t ix)
{
  uint64_t b_curr = lix;
  uint64_t b_prev = std::numeric_limits<uint64_t>::max();

  if (rix >= lix + tau && lane_at(fdps_v[rix], ix) < lane_at(fdps_v[lix], ix)) {
    intervals_v[ix].emplace_back(lix, rix);
    return;
  }

  for (uint64_t a = lix; a <= rix; ++a) {
    const double fdpmax_a = lane_at(fdpmax_v[a - 1], ix);
    const double fdps_a = lane_at(fdps_v[a], ix);

    if (fdpmax_a >= fdps_a) continue; // Condition 1

    while ((b_curr + 1) <= rix && (lane_at(fdsmin_v[b_curr + 1], ix) < fdps_a))
      ++b_curr; // Condition 2

    const uint64_t b_star = b_curr;
    if (b_star < (a + tau)) continue;                    // Condition 3
    if (lane_at(fdps_v[b_star], ix) >= fdps_a) continue; // Condition 4
    if (b_star == b_prev) continue;                      // Condition 6

    if (lane_at(fdps_v[b_star], ix) >= fdpmax_a) { // Condition 5
      intervals_v[ix].emplace_back(a, b_star);
      b_prev = b_star;
    }
  }
}

template<typename T>
void IntExt<T>::expand_intervals(double chisq_th, size_t ix)
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
    fdiff = lane_at(fdps_v[b], ix) - lane_at(fdps_v[ap], ix);
    sdiff = lane_at(sdps_v[ap], ix) - lane_at(sdps_v[b], ix);
    chisq_val = ((fdiff * fdiff) + eps) / (sdiff + eps);

    // Never merge across an N-run break: a skip bin in the gap (bp, a) keeps them apart.
    bool skip_gap = false;
    if (wskip) {
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
void IntExt<T>::compute_prefhistsum()
{
  hist.compute_prefhistsum();
}

template<typename T>
void IntExt<T>::total_histogram(vec<uint64_t>& v, uint64_t& u, uint64_t& t) const
{
  hist.extract_histogram(0, nbins, v, u, t);
}

template<typename T>
void IntExt<T>::extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const
{
  hist.extract_histogram(a, b, v, u, t);
}

template<typename T>
void IntExt<T>::extract_histogram(uint64_t a, uint64_t b, window_counts_t& wc) const
{
  uint64_t t = 0;
  extract_histogram(a, b, wc.hist_v, wc.u, t);
}

template<typename T>
void IntExt<T>::extract_histogram(uint64_t a, uint64_t b, swindow_counts_t& wc, bool is_rc) const
{
  uint64_t t = 0;
  if (is_rc)
    extract_histogram(a, b, wc.hist_rc_v, wc.u_rc, t);
  else
    extract_histogram(a, b, wc.hist_fw_v, wc.u_fw, t);
}

template class IntExt<double>;
template class IntExt<cmlane_t>;
