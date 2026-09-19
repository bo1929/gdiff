#ifndef _TSV_HPP
#define _TSV_HPP

// Tab-separated output helpers: the variadic row writer and the "(low, high)"
// distance bracket that inherits the stream's flags and precision.

#include <cstddef>
#include <ostream>

template<typename... Args>
inline std::ostream& write_tsv(std::ostream& os, const Args&... args)
{
  size_t n = 0;
  ((os << (n++ ? "\t" : "") << args), ...);
  return os;
}

// A distance bracket rendered as "(low, high)".
struct bracket_t
{
  double lo;
  double hi;
};

inline std::ostream& operator<<(std::ostream& os, const bracket_t& b) { return os << '(' << b.lo << ", " << b.hi << ')'; }

#endif
