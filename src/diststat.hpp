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

// When defined, Metropolis-Hastings MCMC adds posterior draws.
// Otherwise, only the samples are from MLE based on mesh sampling.
// #define MCMC

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

private:
  void filter_sample(const record_t& r, vec<p_t>& p_v, uint64_t sample_size) const;
  bool sample_metropolis_hastings(const p_t& init, vec<p_t>& p_v);
  static double sample_box_muller(std::mt19937& rng);

  static constexpr uint64_t S = 20;  // TODO: Revisit: MH MCMC parameters.
  static constexpr uint64_t B = 100; // TODO: Revisit: MH MCMC parameters.

  const params_t<T>& params;
  const llh_sptr_t<T> llhf;
  vec<sample_t> samples_v;
};

#endif
