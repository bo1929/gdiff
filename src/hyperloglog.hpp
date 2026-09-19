#ifndef _HYPERLOGLOG_HPP
#define _HYPERLOGLOG_HPP

// HyperLogLog cardinality estimator (Hideaki Ohno, 2013).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace hll {

  class HyperLogLog
  {
  public:
    static constexpr uint8_t min_precision = 8;
    static constexpr uint8_t max_precision = 18;

    explicit HyperLogLog(uint8_t precision)
      : p(precision)
      , q(static_cast<uint8_t>(64 - precision))
      , m(uint32_t(1) << precision)
      , registers(size_t(1) << precision, 0)
    {
    }

    void add(uint64_t hash) noexcept
    {
      const uint32_t ix = static_cast<uint32_t>(hash >> q);
      const uint64_t suffix = hash << p;
      // clzll(0) is undefined; an all-zero suffix means the maximal rank.
      const uint8_t rank = suffix ? static_cast<uint8_t>(__builtin_clzll(suffix) + 1) : static_cast<uint8_t>(q + 1);
      if (rank > registers[ix]) registers[ix] = rank;
    }

    [[nodiscard]] double estimate() const
    {
      // Register-value multiplicities; values are in [0, q+1] by construction.
      std::vector<uint32_t> mult(static_cast<size_t>(q) + 2, 0);
      for (const uint8_t r : registers)
        ++mult[r];

      const double md = static_cast<double>(m);
      double z = md * tau(1.0 - static_cast<double>(mult[q + 1]) / md);
      for (uint32_t k = q; k >= 1; --k)
        z = 0.5 * (z + static_cast<double>(mult[k]));
      z += md * sigma(static_cast<double>(mult[0]) / md);
      if (!(z > 0.0)) return 0.0;
      return alpha_inf * md * md / z;
    }

  private:
    // 1 / (2 ln 2): the m -> infinity limit of the HLL bias constant.
    static constexpr double alpha_inf = 0.7213475204444817;

    static double sigma(double x)
    {
      if (x == 1.0) return std::numeric_limits<double>::infinity();
      double y = 1.0;
      double z = x;
      double z_prev;
      do {
        x *= x;
        z_prev = z;
        z += x * y;
        y *= 2.0;
      } while (z != z_prev);
      return z;
    }

    static double tau(double x)
    {
      if (x == 0.0 || x == 1.0) return 0.0;
      double y = 1.0;
      double z = 1.0 - x;
      double z_prev;
      do {
        x = std::sqrt(x);
        z_prev = z;
        y *= 0.5;
        const double d = 1.0 - x;
        z -= d * d * y;
      } while (z != z_prev);
      return z / 3.0;
    }

    uint8_t p;
    uint8_t q;
    uint32_t m;
    std::vector<uint8_t> registers;
  };

} // namespace hll

#endif