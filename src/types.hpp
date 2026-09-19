#ifndef _TYPES_HPP
#define _TYPES_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

inline constexpr size_t rwidth = 8; // AVX-512 double lanes: the maximum threshold width

using cmlane_t = std::array<double, rwidth>; // multi-threshold SIMD lanes

class LSHF;
class Buckets;
class Sketch;
class Container;

using enc_t = uint32_t;
using str = std::string;
using strstream = std::stringstream;

struct interval_t
{
  uint64_t a;
  uint64_t b;

  interval_t() = default;
  interval_t(uint64_t a, uint64_t b)
    : a(a)
    , b(b)
  {
  }
};

using xy_t = std::pair<double, double>;
using lshf_sptr_t = std::shared_ptr<LSHF>;

template<typename T, size_t WIDTH>
using arr = std::array<T, WIDTH>;

template<typename T>
using vec = std::vector<T>;

template<typename T>
using vvec = std::vector<std::vector<T>>;

#endif
