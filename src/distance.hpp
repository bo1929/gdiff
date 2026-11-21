#ifndef _DISTANCE_HPP
#define _DISTANCE_HPP

// The distance domain, the strand conventions built on it, and small numeric helpers.

#include "types.hpp"
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

// --- Numeric helpers -----------------------------------------------------------------------------

inline double nanx() noexcept { return std::numeric_limits<double>::quiet_NaN(); }
inline double pinf() noexcept { return std::numeric_limits<double>::infinity(); }
inline double ninf() noexcept { return -std::numeric_limits<double>::infinity(); }

// Linear-interpolation quantile of a sorted vector; NaN when empty.
inline double linear_quantile(const vec<double>& v, double p)
{
  if (v.empty()) return nanx();
  if (v.size() == 1) return v.front();
  const double ix = p * static_cast<double>(v.size() - 1);
  const size_t lo = static_cast<size_t>(std::floor(ix));
  const size_t hi = static_cast<size_t>(std::ceil(ix));
  return v[lo] + (ix - static_cast<double>(lo)) * (v[hi] - v[lo]);
}

// --- Distance domain -----------------------------------------------------------------------------

static constexpr uint32_t hdist_bound = 7;

// Reachable distances: Brent's search interval for the MLE, and the range is_valid_distance accepts.
static constexpr double LB = 0.0 + 1e-7;
static constexpr double UB = 0.5 - 1e-7;

// Padding for unused SIMD threshold lanes.
static constexpr double d_eps = 1e-7;

inline bool is_valid_distance(double d) noexcept { return std::isfinite(d) && d >= LB && d <= UB; }

inline double validate_distance(double d) noexcept { return is_valid_distance(d) ? d : nanx(); }

inline bool is_high_side(double d1, double d0) noexcept { return is_valid_distance(d0) && d1 > d0; }

// --- Strand conventions --------------------------------------------------------------------------

// Encodes which strand is closer: NaN = neither distance is finite.
// -inf = only fw is finite; +inf = only rc is finite, otherwise fw - rc.
inline double strand_diff(double d_fw, double d_rc) noexcept
{
  const bool fw_valid = std::isfinite(d_fw);
  const bool rc_valid = std::isfinite(d_rc);
  if (!fw_valid && !rc_valid) return nanx();
  if (fw_valid && !rc_valid) return ninf();
  if (!fw_valid && rc_valid) return pinf();
  return d_fw - d_rc;
}

// The reported strand symbol of a record on strand `is_rc`, given the sequence's strand_diff.
inline char report_strand(bool is_rc, double d_diff) noexcept
{
  if (std::isnan(d_diff)) {
    return '.';
  }
  if (std::isinf(d_diff)) {
    if (d_diff < 0.0) {
      return is_rc ? '.' : '+';
    } else {
      return is_rc ? '+' : '.';
    }
  }
  return (is_rc == (d_diff > 0.0)) ? '+' : '-';
}

// The closer of the two strand distances and its symbol; fw wins ties, '.' when neither is valid.
inline std::pair<double, char> select_strand_distance(double d_fw, double d_rc)
{
  const bool fw_valid = is_valid_distance(d_fw);
  const bool rc_valid = is_valid_distance(d_rc);
  if (!fw_valid && !rc_valid) return {nanx(), '.'};
  if (!rc_valid || (fw_valid && d_fw <= d_rc)) return {d_fw, '+'};
  return {d_rc, '-'};
}

#endif
