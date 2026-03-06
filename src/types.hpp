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

typedef uint64_t inc_t;
typedef uint32_t enc_t;
typedef std::string str;
typedef std::stringstream strstream;
typedef std::pair<uint64_t, uint64_t> interval_t;
typedef std::shared_ptr<RSeq> rseq_sptr_t;
typedef std::shared_ptr<QSeq> qseq_sptr_t;
typedef std::shared_ptr<LSHF> lshf_sptr_t;
typedef std::shared_ptr<SDHM> sdhm_sptr_t;
typedef std::shared_ptr<SFHM> sfhm_sptr_t;
typedef std::shared_ptr<Sketch> sketch_sptr_t;

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
  size_t n;
  T dist_th;
  uint32_t hdist_th;
  uint64_t min_length;
  double chisq;
  uint64_t bin_shift;
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
