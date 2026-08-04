#ifndef _DISTSTAT_HPP
#define _DISTSTAT_HPP

#include "dim.hpp"
#include "gamma.hpp"
#include "msg.hpp"
#include "llh.hpp"
#include "maptils.hpp"
#include "random.hpp"
#include "types.hpp"
#include <algorithm>
#include <cmath>
#include <random>

template<typename T>
class DistanceStat
{
public:
  DistanceStat(const params_t<T>& params, const llh_sptr_t<T>& llhf);

  void clear_samples();
  void sample_null_pool(const DIM<T>& dim, uint64_t tau_eff, uint64_t bix);
  bool test_significance(record_t& r, uint64_t sample_size, const str& qid);
  void benjamini_hochberg_correction(vec<record_t>& records);

  const vec<sample_t>& samples() const { return samples_v; }
  uint64_t get_n_unmapped() const { return n_unmapped; }

private:
  // Fills p_v with the pool samples usable for r (excludes same-query windows
  // overlapping r, downsamples to sample_size). Returns true when any pool
  // sample was excluded for overlap (i.e. the set is record-specific).
  bool filter_sample(const record_t& r, vec<p_t>& p_v, uint64_t sample_size) const;

  const params_t<T>& params;
  const llh_sptr_t<T> llhf;
  vec<sample_t> samples_v;
  uint64_t n_unmapped = 0; // sampled null windows with no k-mer hits (t == 0)

  // Gamma fit cache: records of one query (bix) are contiguous and usually
  // share the same filtered null set, so the Nelder-Mead fit and the latent
  // median are computed once per query unless an overlap exclusion applied.
  uint64_t fit_cache_bix = std::numeric_limits<uint64_t>::max();
  GammaModel::params_t fit_cache_params{1.0, 1.0};
  double fit_cache_median = nanx();
  bool fit_cache_ok = false;
};

#endif
