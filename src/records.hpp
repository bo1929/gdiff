#ifndef _RECORDS_HPP
#define _RECORDS_HPP

// The per-interval output records, the per-window background-sample row, and the
// interval/coordinate helpers that turn bin ranges into reported coordinates.

#include "distance.hpp"
#include "types.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

// Sentinel for `record_t::th_ix`: no threshold (a background/unreported row).
inline constexpr size_t no_threshold = std::numeric_limits<size_t>::max();

struct record_t
{
  uint64_t bix;      // batch index of the source query
  uint64_t L;        // effective query length (enmers + k - 1)
  interval_t seq_iv; // 1-based inclusive coordinates on query
  uint64_t nbins;    // number of bins covered: bin_iv.b - bin_iv.a
  interval_t bin_iv; // 1-based inclusive bin start, 1-based exclusive bin end
  bool is_rc;        // Is the source query on the reverse-complement strand?
  double d;          // MLE distance for this interval in [d_eps, d_ub]
  double I;          // Observed Fisher information I(d)
  size_t th_ix;      // Threshold index; no_threshold marks background gaps (not reported anymore)
  // Per-query fields, filled once both strands' distances are known
  double d_q = nanx();        // MLE distance for the source strand
  double d_diff = nanx();     // strand difference encoding; see strand_diff()
  double fold = nanx();       // fold change: d / median(background samples)
  double percentile = nanx(); // two-sided percentile for the closer strand (reference), otherwise cdf
  double qvalue = nanx();     // Benjamini-Hochberg adjusted percentile
  double lr_bg;               // likelihood-ratio statistic vs the background distance
  double lr_ub;               // likelihood-ratio statistic vs max estimable distance (extreme match counts)

  bool is_intact() const { return seq_iv.a == 1 && seq_iv.b == L; }
  interval_t get_interval() const { return {bin_iv.a - 1, bin_iv.b - 1}; } // 0-based half-open bin-boundary

  record_t(uint64_t bix,
           uint64_t L,
           interval_t seq_iv,
           interval_t bin_iv,
           bool is_rc,
           double d,
           double I,
           size_t th_ix,
           double lr_bg = nanx(),
           double lr_ub = nanx())
    : bix(bix)
    , L(L)
    , seq_iv(seq_iv)
    , nbins(bin_iv.b - bin_iv.a)
    , bin_iv(bin_iv)
    , is_rc(is_rc)
    , d(d)
    , I(I)
    , th_ix(th_ix)
    , lr_bg(lr_bg)
    , lr_ub(lr_ub)
  {
  }
};

// One sample pool entry: a window sampled from a query sequence (bix)
struct sample_t
{
  double d;          // MLE distance
  double I;          // Fisher information
  uint64_t bix;      // query batch index (for overlap filtering)
  interval_t bin_iv; // 1-based half-open bin coordinates
};

// 1-based half-open bin range [a_bin, b_bin) plus the covering threshold.
struct bp_t
{
  uint64_t a_bin;
  uint64_t b_bin;
  size_t ix;
};

// The thresholds bracketing a distance, for the d_bin output column; the full distance
// domain when d is not valid.
inline xy_t bracket_distance(double d, const vec<double>& th_v)
{
  xy_t d_range{d_eps, d_ub};
  if (!is_valid_distance(d)) return d_range;
  const auto it = std::lower_bound(th_v.begin(), th_v.end(), d);
  if (it != th_v.begin()) d_range.first = *(it - 1);
  if (it != th_v.end()) d_range.second = *it;
  return d_range;
}

// 1-based inclusive bp coordinates of the bin range [bin_iv.a, bin_iv.b).
inline interval_t get_coordinates(const interval_t& bin_iv, uint64_t bin_shift, uint64_t enmers, uint32_t k)
{
  const uint64_t a = ((bin_iv.a - 1) << bin_shift) + 1;
  const uint64_t b = std::min((bin_iv.b - 1) << bin_shift, enmers) + k - 1;
  return {a, b};
}

// 1-based half-open interval convention: inclusive start, exclusive end.
inline bool overlaps_half_open(const interval_t& lhs, const interval_t& rhs) { return lhs.a < rhs.b && rhs.a < lhs.b; }

// Per-strand Benjamini-Hochberg adjustment of record percentiles into qvalues.
inline void benjamini_hochberg_correction(vec<record_t>& records_v)
{
  for (bool is_rc : {false, true}) {
    vec<size_t> idx_v;
    for (size_t i = 0; i < records_v.size(); ++i)
      if (records_v[i].is_rc == is_rc && std::isfinite(records_v[i].percentile)) idx_v.push_back(i);
    if (idx_v.empty()) continue;

    std::sort(
      idx_v.begin(), idx_v.end(), [&](size_t a, size_t b) { return records_v[a].percentile < records_v[b].percentile; });

    const double m = static_cast<double>(idx_v.size());
    double q_min = 1.0;
    for (size_t rank = idx_v.size(); rank >= 1; --rank) {
      record_t& r = records_v[idx_v[rank - 1]];
      q_min = std::min(q_min, std::min(1.0, r.percentile * m / static_cast<double>(rank)));
      r.qvalue = q_min;
    }
  }
}

#endif
