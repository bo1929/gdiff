#ifndef _RECORDS_HPP
#define _RECORDS_HPP

// The per-interval output records, the per-window background-sample line, and the
// interval/coordinate helpers that turn bin ranges into reported coordinates.

#include "distance.hpp"
#include "types.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ostream>

// Tab-separated output helpers: the variadic line writer and the "(low, high)" distance
// bracket that inherits the stream's flags and precision.
template<typename... Args>
inline std::ostream& write_tsv(std::ostream& os, const Args&... args)
{
  size_t n = 0;
  ((os << (n++ ? "\t" : "") << args), ...);
  return os;
}

struct bracket_t
{
  double lo;
  double hi;
};

inline std::ostream& operator<<(std::ostream& os, const bracket_t& b) { return os << '(' << b.lo << ", " << b.hi << ')'; }

// Sentinel for `record_t::th_ix`: no threshold (a background/unreported line).
inline constexpr size_t no_threshold = std::numeric_limits<size_t>::max();

// One bit per threshold lane (bit i = lane i); widened past the eight lanes it needs so it
// streams as a number rather than as an `unsigned char`.
using mask_t = uint32_t;

struct record_t
{
  uint64_t bix;       // batch index of the source query
  uint64_t L;         // effective query length (enmers + k - 1)
  interval_t seq_iv;  // 1-based inclusive coordinates on query
  bool is_rc;         // Is the source query on the reverse-complement strand?
  double d;           // MLE distance for this interval in [d_eps, d_ub]
  double fisher_info; // Observed Fisher information at d
  size_t th_ix;       // Threshold index; no_threshold marks background gaps (not reported anymore)
  // Per-sequence fields, filled once both strands' distances are known
  double d_seq = nanx();  // MLE distance for the source strand (its own background)
  double d_diff = nanx(); // strand difference encoding; see strand_diff()
  double fold = nanx();   // fold change: d / median(background samples)
  double pvalue = nanx(); // empirical two-sided p-value against the background pool
  double qvalue = nanx(); // Benjamini-Hochberg adjusted p-value
  // Likelihood-ratio statistics.
  double lr_ub; // this interval vs the max estimable distance (extreme match counts)
  double lr_lb; // this interval vs perfect identity (all k-mers matching)
  double lr_null = nanx();     // this interval vs the sampled null center
  double p_interval = nanx();  // chi-square p-value of lr_null
  double q_interval = nanx();  // BH-adjusted p_interval, per strand
  double p_gamma = nanx();     // two-sided p-value from the fitted hurdle gamma

  record_t(uint64_t bix,
           uint64_t L,
           interval_t seq_iv,
           bool is_rc,
           double d,
           double fisher_info,
           size_t th_ix,
           double lr_ub = nanx(),
           double lr_lb = nanx())
    : bix(bix)
    , L(L)
    , seq_iv(seq_iv)
    , is_rc(is_rc)
    , d(d)
    , fisher_info(fisher_info)
    , th_ix(th_ix)
    , lr_ub(lr_ub)
    , lr_lb(lr_lb)
  {
  }
};

// The thresholds bracketing a distance, for the d_bin output column; the full distance
// domain when d is not valid. `th_v` must be sorted ascending.
inline xy_t bracket_distance(double d, const vec<double>& th_v)
{
  assert(std::is_sorted(th_v.begin(), th_v.end()));
  xy_t d_range{d_eps, d_ub};
  if (!is_valid_distance(d)) return d_range;
  const auto it = std::lower_bound(th_v.begin(), th_v.end(), d);
  if (it != th_v.begin()) d_range.first = *(it - 1);
  if (it != th_v.end()) d_range.second = *it;
  return d_range;
}

// 1-based inclusive bp coordinates of the bin range [bin_iv.a, bin_iv.b); the bin coordinates
// are 1-based, so both must be at least 1 or the shift underflows.
inline interval_t get_coordinates(const interval_t& bin_iv, uint64_t bin_shift, uint64_t enmers, uint32_t k)
{
  assert(bin_iv.a >= 1 && bin_iv.b >= 1);
  const uint64_t a = ((bin_iv.a - 1) << bin_shift) + 1;
  const uint64_t b = std::min((bin_iv.b - 1) << bin_shift, enmers) + k - 1;
  return {a, b};
}

// Minimum background windows for a p-value to mean anything. `map` raises this to a tenth of its
// `--sample-size`; 8 is the floor for callers that do not carry a sample size.
inline constexpr size_t min_null_samples = 8;

// Empirical two-sided p-value of r.d within a sorted background pool; sets pvalue and fold. Both
// tails are treated alike, so a record that is unusually *similar* and one that is unusually
// *divergent* are scored on the same scale - the scale `benjamini_hochberg_correction` consumes,
// which is why the reported column and the q-value can never disagree.
inline void apply_empirical_significance(record_t& r, const vec<double>& pool_v, size_t min_samples = min_null_samples)
{
  if (!is_valid_distance(r.d) || pool_v.size() < min_samples) return;
  assert(std::is_sorted(pool_v.begin(), pool_v.end()));
  const size_t le = static_cast<size_t>(std::upper_bound(pool_v.begin(), pool_v.end(), r.d) - pool_v.begin());
  const double prob = static_cast<double>(le) / static_cast<double>(pool_v.size());
  r.pvalue = std::clamp(2.0 * std::min(prob, 1.0 - prob), 0.0, 1.0);
  const double median = linear_quantile(pool_v, 0.5);
  if (median > eps) r.fold = r.d / median;
}

// Per-strand Benjamini-Hochberg adjustment of the interval-level p-values.
inline void benjamini_hochberg_intervals(vec<record_t>& records_v)
{
  for (bool is_rc : {false, true}) {
    vec<size_t> idx_v;
    for (size_t i = 0; i < records_v.size(); ++i)
      if (records_v[i].is_rc == is_rc && std::isfinite(records_v[i].p_interval)) idx_v.push_back(i);
    if (idx_v.empty()) continue;

    std::sort(idx_v.begin(), idx_v.end(), [&](size_t a, size_t b) { return records_v[a].p_interval < records_v[b].p_interval; });

    const double m = static_cast<double>(idx_v.size());
    double q_min = 1.0;
    for (size_t i = idx_v.size(); i-- > 0;) {
      record_t& r = records_v[idx_v[i]];
      q_min = std::min(q_min, std::min(1.0, r.p_interval * m / static_cast<double>(i + 1)));
      r.q_interval = q_min;
    }
  }
}

// Per-strand Benjamini-Hochberg adjustment of the p-values into qvalues.
inline void benjamini_hochberg_correction(vec<record_t>& records_v)
{
  for (bool is_rc : {false, true}) {
    vec<size_t> idx_v;
    for (size_t i = 0; i < records_v.size(); ++i)
      if (records_v[i].is_rc == is_rc && std::isfinite(records_v[i].pvalue)) idx_v.push_back(i);
    if (idx_v.empty()) continue;

    std::sort(idx_v.begin(), idx_v.end(), [&](size_t a, size_t b) { return records_v[a].pvalue < records_v[b].pvalue; });

    const double m = static_cast<double>(idx_v.size());
    double q_min = 1.0;
    // Walk from the largest p down; `rank` counts from 1, so `i + 1` is the rank.
    for (size_t i = idx_v.size(); i-- > 0;) {
      record_t& r = records_v[idx_v[i]];
      q_min = std::min(q_min, std::min(1.0, r.pvalue * m / static_cast<double>(i + 1)));
      r.qvalue = q_min;
    }
  }
}

#endif
