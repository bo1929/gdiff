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

template<typename T>
class LLH;
template<typename T>
class DIM;
template<typename T>
class QIE;
class RSeq;
class QSeq;
class LSHF;
class SDHM;
class SFHM;
class Sketch;

// using inc_t = uint32_t; // This might be just OK...
using inc_t = uint64_t;
using enc_t = uint32_t;
using str = std::string;
using strstream = std::stringstream;
using interval_t = std::pair<uint64_t, uint64_t>;
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

template<typename T>
using llh_sptr_t = std::shared_ptr<LLH<T>>;

using cm512_t = std::array<double, RWIDTH>;

struct hmer_t
{
  uint64_t x, y, z;
};

template<typename T>
struct params_t
{
  size_t n;           // Number of distance thresholds given, also equals to WIDTH later
  T dist_th;          // Distance threshold used for detection across varying scales
  uint32_t hdist_th;  // Hamming distance threshold used for k-mer search
  double chisq;       // Chi-square threshold in the statistical test for interval merging
  uint64_t bin_shift; // Shift value for fast bin index calculation
  uint64_t bin_size;  // Bin size in sites, equals to pow(2, bin_shift)
  uint64_t tau;       // The minimum length threshold in sites
  uint64_t tau_bin;   // The minimum length threshold in number of bins instead of sites
  bool enum_only;

  params_t(size_t n, T dist_th, uint32_t hdist_th, uint64_t tau, double chisq, uint64_t bin_shift, bool enum_only)
    : n(n)
    , dist_th(dist_th)
    , hdist_th(hdist_th)
    , tau(tau)
    , chisq(chisq)
    , bin_shift(bin_shift)
    , bin_size(uint64_t(1) << bin_shift)
    , tau_bin((tau + (uint64_t(1) << bin_shift) - 1) >> bin_shift)
    , enum_only(enum_only)
  {
  }
};

// struct alignas(64) cm512_t
// {
//   arr<double, RWIDTH> v{};
// };

// #define EXTRAARGS                                                                                                           \
//   phmap::priv::hash_default_hash<K>, phmap::priv::hash_default_eq<K>, std::allocator<std::pair<const K, V>>, 4

// template<class K, class V>
// using parallel_flat_phmap = phmap::parallel_flat_hash_map<K, V, EXTRAARGS, std::mutex>;

// template<class K, class V>
// using parallel_node_phmap = phmap::parallel_node_hash_map<K, V, EXTRAARGS, std::mutex>;

// template<class K, class V>
// using fparallel_flat_phmap = phmap::parallel_flat_hash_map<K, V, EXTRAARGS>;

// template<class K, class V>
// using fparallel_node_phmap = phmap::parallel_node_hash_map<K, V, EXTRAARGS>;

// template<class K, class V>
// using flat_phmap = phmap::flat_hash_map<K, V>;

// template<class K, class V>
// using btree_phmap = phmap::btree_map<K, V>;

// template<class K, class V>
// using node_phmap = phmap::node_hash_map<K, V>;

#endif
