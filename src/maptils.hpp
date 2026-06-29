#ifndef _MAPTILS_HPP
#define _MAPTILS_HPP

#include "types.hpp"
#include <limits>
#include <sstream>
#include <utility>

static constexpr uint32_t hdist_bound = 7; // Hamming distance bound for the SIMD alignment

static constexpr double d_ub = 1.0;
static constexpr double d_lb = 0.0;
static constexpr double d_eps = 0.00001;

static constexpr double eps = 1e-7; // Tolerance for floating-point comparisons

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

inline double validate_distance(const double d)
{
  if (d >= d_ub - eps || std::isnan(d) || (d < d_lb)) {
    return nanx();
  }
  return d;
}

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

struct p_t
{
  double d; // MLE distance
  double I; // Fisher information
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
  double fold = nanx();       // fold change: d / median(null samples)
  double percentile = nanx(); // two-sided percentile for the closer strand (reference), otherwise cdf
  double qvalue = nanx();     // Benjamini-Hochberg adjusted percentile

  bool is_intact() const { return seq_iv.a == 1 && seq_iv.b == L; }
  interval_t get_interval() const { return {bin_iv.a - 1, bin_iv.b - 1}; } // 0-based half-open bin-boundary

  record_t(uint64_t bix, uint64_t L, interval_t seq_iv, interval_t bin_iv, bool is_rc, double d, double I, size_t th_ix)
    : bix(bix)
    , L(L)
    , seq_iv(seq_iv)
    , nbins(bin_iv.b - bin_iv.a)
    , bin_iv(bin_iv)
    , is_rc(is_rc)
    , d(d)
    , I(I)
    , th_ix(th_ix)
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

inline interval_t get_coordinates(const interval_t& bin_iv, uint64_t bin_shift, uint64_t enmers, uint32_t k, bool is_last)
{
  const uint64_t a = ((bin_iv.a - 1) << bin_shift) + 1;
  const uint64_t b = is_last ? (enmers + k - 1) : std::min(((bin_iv.b - 1) << bin_shift) + 1, enmers + 1);
  return {a, b};
}

template<typename... Args>
inline std::ostream& write_tsv(std::ostream& os, const Args&... args)
{
  size_t n = 0;
  ((os << (n++ ? "\t" : "") << args), ...);
  return os;
}

#endif
