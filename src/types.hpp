#ifndef _TYPES_HPP
#define _TYPES_HPP

#include <array>
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
// #include "btree.h"
// #include "phmap.h"

#define RWIDTH 8

class RSeq;
class QSeq;
class LSHF;
class SDHM;
class SFHM;
class Sketch;

// using inc_t = uint32_t; // This might be just OK
using inc_t = uint64_t;
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
using rseq_sptr_t = std::shared_ptr<RSeq>;
using qseq_sptr_t = std::shared_ptr<QSeq>;
using lshf_sptr_t = std::shared_ptr<LSHF>;
using sdhm_sptr_t = std::shared_ptr<SDHM>;
using sfhm_sptr_t = std::shared_ptr<SFHM>;
using sketch_sptr_t = std::shared_ptr<Sketch>;

template<typename T, size_t WIDTH>
using arr = std::array<T, WIDTH>;

template<typename T>
using vec = std::vector<T>;

template<typename T>
using vvec = std::vector<std::vector<T>>;

using cm512_t = std::array<double, RWIDTH>;

template<typename T>
struct params_t
{
  T dist_th;            // Distance threshold used for detection across varying scales
  uint32_t hdist_th;    // Hamming distance threshold used for k-mer search
  uint64_t tau_bin;     // The minimum length threshold in number of bins instead of sites
  double chisq;         // Chi-square threshold in the statistical test for interval merging
  uint64_t bin_shift;   // Shift value for fast bin index calculation
  uint64_t bin_size;    // Bin size in sites, equals to pow(2, bin_shift)
  uint64_t sample_size; // Number of background samples for significance test (0 = skip)
  bool canonical;       // Strand-agnostic sketch mode and canonical k-mers
  bool enum_only;       // Only enumerate intervals, no iterative interval removal

  params_t(T dist_th,
           uint32_t hdist_th,
           uint64_t tau,
           double chisq,
           uint64_t bin_shift,
           uint64_t sample_size,
           bool canonical,
           bool enum_only)
    : dist_th(dist_th)
    , hdist_th(hdist_th)
    , tau_bin((tau + (uint64_t(1) << bin_shift) - 1) >> bin_shift)
    , chisq(chisq)
    , bin_shift(bin_shift)
    , bin_size(uint64_t(1) << bin_shift)
    , sample_size(sample_size)
    , canonical(canonical)
    , enum_only(enum_only)
  {
  }
};

#endif
