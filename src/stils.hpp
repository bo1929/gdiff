#ifndef _STILS_HPP
#define _STILS_HPP

#include "types.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>

static constexpr uint32_t hdist_bound = 7;

static constexpr double d_ub = 1.0;
static constexpr double d_lb = 0.0;
static constexpr double d_eps = 0.00001;

// MLE search domain (Brent bounds); UB is the no-homology plateau end.
static constexpr double UB = 0.99;
static constexpr double LB = 0.0;

static constexpr double eps = 1e-7;

inline double nanx() noexcept { return std::numeric_limits<double>::quiet_NaN(); }
inline double pinf() noexcept { return std::numeric_limits<double>::infinity(); }
inline double ninf() noexcept { return -std::numeric_limits<double>::infinity(); }

inline double strand_diff(const double d_q_fw, const double d_q_rc) noexcept
{
  const bool fw_valid = std::isfinite(d_q_fw);
  const bool rc_valid = std::isfinite(d_q_rc);
  if (!fw_valid && !rc_valid) return nanx(); // NaN = neither finite
  if (fw_valid && !rc_valid) return ninf();  // -inf = only fw is finite
  if (!fw_valid && rc_valid) return pinf();  // inf = only rc is finite
  return d_q_fw - d_q_rc;                    // both are finite
}

inline char report_strand(const bool is_rc, const double d_diff) noexcept
{
  if (std::isnan(d_diff)) return '.';
  if (std::isinf(d_diff)) {
    if (d_diff < 0.0) return is_rc ? '.' : '+';
    return is_rc ? '+' : '.';
  }
  return (is_rc == (d_diff > 0.0)) ? '+' : '-';
}

inline bool is_valid_distance(const double d) noexcept { return std::isfinite(d) && d >= d_lb && d < d_ub - eps; }

inline double validate_distance(const double d) noexcept { return is_valid_distance(d) ? d : nanx(); }

// Per-window HD counts (single strand / canonical): scan aggregator and
// DIM/HDHist extract sink. No fw/rc split.
struct window_counts_t
{
  vec<uint64_t> hist_v;
  uint64_t u = 0;
  uint32_t hdist_th = 0;

  window_counts_t() = default;

  explicit window_counts_t(uint32_t hdist_th)
    : hdist_th(hdist_th)
  {
    hist_v.assign(hdist_bound + 1, 0);
  }

  void clear() noexcept
  {
    std::fill(hist_v.begin(), hist_v.end(), 0);
    u = 0;
  }

  uint64_t* hist() noexcept { return hist_v.data(); }
  const uint64_t* hist() const noexcept { return hist_v.data(); }

  uint64_t t() const noexcept
  {
    uint64_t sum = 0;
    for (uint64_t c : hist_v)
      sum += c;
    return sum;
  }

  // scan_mers_range aggregator (is_rc ignored; canonical scan never sets it)
  inline void operator()(uint64_t /*bin*/, uint32_t hdist, bool /*is_rc*/) noexcept
  {
    if (hdist <= hdist_th)
      ++hist_v[hdist];
    else
      ++u;
  }

  inline void skip_mer(uint64_t /*bin*/) const noexcept {}
};

// Per-window HD counts with explicit fw/rc halves for strand-aware scans.
struct swindow_counts_t
{
  vec<uint64_t> hist_fw_v;
  vec<uint64_t> hist_rc_v;
  uint64_t u_fw = 0;
  uint64_t u_rc = 0;
  uint32_t hdist_th = 0;

  swindow_counts_t() = default;

  explicit swindow_counts_t(uint32_t hdist_th)
    : hdist_th(hdist_th)
  {
    hist_fw_v.assign(hdist_bound + 1, 0);
    hist_rc_v.assign(hdist_bound + 1, 0);
  }

  void clear() noexcept
  {
    std::fill(hist_fw_v.begin(), hist_fw_v.end(), 0);
    std::fill(hist_rc_v.begin(), hist_rc_v.end(), 0);
    u_fw = 0;
    u_rc = 0;
  }

  uint64_t* hist_fw() noexcept { return hist_fw_v.data(); }
  uint64_t* hist_rc() noexcept { return hist_rc_v.data(); }
  const uint64_t* hist_fw() const noexcept { return hist_fw_v.data(); }
  const uint64_t* hist_rc() const noexcept { return hist_rc_v.data(); }

  uint64_t t_fw() const noexcept
  {
    uint64_t sum = 0;
    for (uint64_t c : hist_fw_v)
      sum += c;
    return sum;
  }

  uint64_t t_rc() const noexcept
  {
    uint64_t sum = 0;
    for (uint64_t c : hist_rc_v)
      sum += c;
    return sum;
  }

  // scan_mers_range aggregator
  inline void operator()(uint64_t /*bin*/, uint32_t hdist, bool is_rc) noexcept
  {
    if (hdist <= hdist_th)
      ++(is_rc ? hist_rc_v[hdist] : hist_fw_v[hdist]);
    else
      ++(is_rc ? u_rc : u_fw);
  }

  inline void skip_mer(uint64_t /*bin*/) const noexcept {}
};

// 1-based half-open interval convention: inclusive start, exclusive end.
inline bool overlaps_half_open(const interval_t& lhs, const interval_t& rhs) { return lhs.a < rhs.b && rhs.a < lhs.b; }

// One sample pool entry: a window sampled from a query sequence (bix)
struct sample_t
{
  double d;          // MLE distance
  double I;          // Fisher information
  uint64_t bix;      // query batch index (for overlap filtering)
  interval_t bin_iv; // 1-based half-open bin coordinates
};

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
  size_t th_ix;      // Threshold index; size_t(-1) reserved for background gaps (not reported anymore)
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

// Breakpoints of an extracted interval
// 1-based half-open bin range [a_bin, b_bin) with the threshold that covers it.
struct bp_t
{
  uint64_t a_bin;
  uint64_t b_bin;
  size_t ix;
};

// 1-based inclusive bp coordinates of the bin range [bin_iv.a, bin_iv.b).
// The end always covers the full span of the last k-mer: bins [a, b) cover mer
// starts [(a-1)<<bin_shift, (b-1)<<bin_shift), so the last covered mer starts
// at ((b-1)<<bin_shift)-1 and its k-mer ends at min((b-1)<<bin_shift, enmers)+k-1.
inline interval_t get_coordinates(const interval_t& bin_iv, uint64_t bin_shift, uint64_t enmers, uint32_t k)
{
  const uint64_t a = ((bin_iv.a - 1) << bin_shift) + 1;
  const uint64_t b = std::min((bin_iv.b - 1) << bin_shift, enmers) + k - 1;
  return {a, b};
}

// Per-strand Benjamini-Hochberg adjustment of record percentiles into qvalues.
inline void benjamini_hochberg_correction(vec<record_t>& records_v)
{
  for (bool is_rc : {false, true}) {
    vec<size_t> idx;
    for (size_t i = 0; i < records_v.size(); ++i)
      if (records_v[i].is_rc == is_rc && std::isfinite(records_v[i].percentile)) idx.push_back(i);
    if (idx.empty()) continue;

    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return records_v[a].percentile < records_v[b].percentile; });

    const double m = static_cast<double>(idx.size());
    double q_min = 1.0;
    for (size_t rank = idx.size(); rank >= 1; --rank) {
      record_t& r = records_v[idx[rank - 1]];
      q_min = std::min(q_min, std::min(1.0, r.percentile * m / static_cast<double>(rank)));
      r.qvalue = q_min;
    }
  }
}

template<typename... Args>
inline std::ostream& write_tsv(std::ostream& os, const Args&... args)
{
  size_t n = 0;
  ((os << (n++ ? "\t" : "") << args), ...);
  return os;
}

#endif
