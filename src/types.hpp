#ifndef _TYPES_HPP
#define _TYPES_HPP

#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include "btree.h"
#include "phmap.h"

class SBatch; // TODO
class RSeq;
class QSeq;
class LSHF;
class SDHM;
class SFHM;
class Sketch;

typedef uint64_t inc_t;
typedef uint32_t enc_t;

typedef std::stringstream strstream;

template<typename T>
using vvec = std::vector<std::vector<T>>;
template<typename T>
using vec = std::vector<T>;

typedef std::shared_ptr<SBatch> sbatch_sptr_t; // TODO
typedef std::shared_ptr<RSeq> rseq_sptr_t;
typedef std::shared_ptr<QSeq> qseq_sptr_t;
typedef std::shared_ptr<LSHF> lshf_sptr_t;
typedef std::shared_ptr<SDHM> sdhm_sptr_t;
typedef std::shared_ptr<SFHM> sfhm_sptr_t;
typedef std::shared_ptr<Sketch> sketch_sptr_t;

struct hmer_t
{
  uint64_t x, y, z;
};

#define EXTRAARGS                                                                                                           \
  phmap::priv::hash_default_hash<K>, phmap::priv::hash_default_eq<K>, std::allocator<std::pair<const K, V>>, 4

template<class K, class V>
using parallel_flat_phmap = phmap::parallel_flat_hash_map<K, V, EXTRAARGS, std::mutex>;

template<class K, class V>
using parallel_node_phmap = phmap::parallel_node_hash_map<K, V, EXTRAARGS, std::mutex>;

template<class K, class V>
using fparallel_flat_phmap = phmap::parallel_flat_hash_map<K, V, EXTRAARGS>;

template<class K, class V>
using fparallel_node_phmap = phmap::parallel_node_hash_map<K, V, EXTRAARGS>;

template<class K, class V>
using flat_phmap = phmap::flat_hash_map<K, V>;

template<class K, class V>
using btree_phmap = phmap::btree_map<K, V>;

template<class K, class V>
using node_phmap = phmap::node_hash_map<K, V>;

#endif
