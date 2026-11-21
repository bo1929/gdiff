#ifndef _RECORDS_HPP
#define _RECORDS_HPP

// The per-interval output record and the helpers that turn bin ranges into reported coordinates
// and adjust the interval p-values.

#include "distance.hpp"
#include "types.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>

struct record_t
{
  uint64_t bix;      // batch index of the source query
  uint64_t L;        // effective query length (enmers + k - 1)
  interval_t seq_iv; // 1-based inclusive coordinates on query
  bool is_rc;        // Is the source query on the reverse-complement strand?
  double d;          // MLE distance for this interval; NaN without hits (see distance.hpp for the range)
  size_t th_ix;      // Index of the threshold lane that claimed the interval
  // Per-sequence fields, filled once both strands' distances are known
  double d_seq = nanx();  // MLE distance for the source strand (its own background)
  double d_diff = nanx(); // strand difference encoding; see strand_diff()
  // Significance
  double p_interval = nanx(); // chi-square p-value of the interval vs the sampled null center
  double q_interval = nanx(); // BH-adjusted p_interval, per strand
  double p_gamma = nanx();    // two-sided p-value from the fitted hurdle gamma

  record_t(uint64_t bix, uint64_t L, interval_t seq_iv, bool is_rc, double d, size_t th_ix)
    : bix(bix)
    , L(L)
    , seq_iv(seq_iv)
    , is_rc(is_rc)
    , d(d)
    , th_ix(th_ix)
  {
  }
};

// 1-based inclusive bp coordinates of the bin range [bin_iv.a, bin_iv.b); the bin coordinates
// are 1-based, so both must be at least 1 or the shift underflows.
inline interval_t get_coordinates(const interval_t& bin_iv, uint64_t bin_shift, uint64_t enmers, uint32_t k)
{
  assert(bin_iv.a >= 1 && bin_iv.b >= 1);
  const uint64_t a = ((bin_iv.a - 1) << bin_shift) + 1;
  const uint64_t b = std::min((bin_iv.b - 1) << bin_shift, enmers) + k - 1;
  return {a, b};
}

// Per-strand Benjamini-Hochberg adjustment of the interval-level p-values.
inline void benjamini_hochberg_intervals(vec<record_t>& records_v)
{
  for (bool is_rc : {false, true}) {
    vec<size_t> idx_v;
    for (size_t i = 0; i < records_v.size(); ++i)
      if (records_v[i].is_rc == is_rc && std::isfinite(records_v[i].p_interval)) idx_v.push_back(i);
    if (idx_v.empty()) continue;

    std::sort(
      idx_v.begin(), idx_v.end(), [&](size_t a, size_t b) { return records_v[a].p_interval < records_v[b].p_interval; });

    const double m = static_cast<double>(idx_v.size());
    double q_min = 1.0;
    for (size_t i = idx_v.size(); i-- > 0;) {
      record_t& r = records_v[idx_v[i]];
      q_min = std::min(q_min, std::min(1.0, r.p_interval * m / static_cast<double>(i + 1)));
      r.q_interval = q_min;
    }
  }
}

#endif
