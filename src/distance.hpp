#ifndef _DISTANCE_HPP
#define _DISTANCE_HPP

// Distance model constants and predicates.

#include "types.hpp"
#include <cmath>
#include <cstddef>

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

inline double strand_diff(double d_q_fw, double d_q_rc) noexcept
{
  const bool fw_valid = std::isfinite(d_q_fw);
  const bool rc_valid = std::isfinite(d_q_rc);
  if (!fw_valid && !rc_valid) return nanx(); // NaN = neither finite
  if (fw_valid && !rc_valid) return ninf();  // -inf = only fw is finite
  if (!fw_valid && rc_valid) return pinf();  // inf = only rc is finite
  return d_q_fw - d_q_rc;                    // both are finite
}

inline char report_strand(bool is_rc, double d_diff) noexcept
{
  if (std::isnan(d_diff)) return '.';
  if (std::isinf(d_diff)) {
    if (d_diff < 0.0) return is_rc ? '.' : '+';
    return is_rc ? '+' : '.';
  }
  return (is_rc == (d_diff > 0.0)) ? '+' : '-';
}

inline bool is_valid_distance(double d) noexcept { return std::isfinite(d) && d >= d_lb && d < d_ub - eps; }

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

inline double validate_distance(double d) noexcept { return is_valid_distance(d) ? d : nanx(); }

inline std::pair<double, char> select_strand_distance(double d_fw, double d_rc)
{
  const bool fw_valid = is_valid_distance(d_fw);
  const bool rc_valid = is_valid_distance(d_rc);
  if (!fw_valid && !rc_valid) return {nanx(), '.'};
  if (!rc_valid || (fw_valid && d_fw <= d_rc)) return {d_fw, '+'};
  return {d_rc, '-'};
}

#endif
