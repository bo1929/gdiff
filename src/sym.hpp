#ifndef _SYM_HPP
#define _SYM_HPP

#include "distance.hpp"
#include "types.hpp"

struct ds_t
{
  double d = nanx(); // the MLE distance estimate
  double s = nanx(); // likelihood-ratio statistic for the upper bound
};

struct summary_t
{
  double d = nanx();         // reported estimate
  double d_median = nanx();  // median over the distances behind `d`
  double d_mean = nanx();    // plain mean over the symmetrized distances, before filtering
  double d_highest = nanx(); // max d over every mapped window
  double d_upper = nanx();   // max d over the kept windows only
  vec<double> d_v;           // distances after symmetrization
  uint64_t n_na = 0;         // # of windows without a distance estimate
  uint64_t n_ub = 0;         // # of windows whose likelihood-ratio statistic is exactly zero
  uint64_t n_filtered = 0;   // # of windows rejected by the likelihood-ratio filter
};

summary_t summarize_symmetric(vec<ds_t> ab_v, vec<ds_t> ba_v, double lr_th, double min_portion);

constexpr double lr_th_default = 3.841;
constexpr double min_portion_default = 0.66;

#endif
