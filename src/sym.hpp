#ifndef _SYM_HPP
#define _SYM_HPP

#include "stils.hpp"
#include "types.hpp"

struct dpoint_t
{
  double d = nanx();
  double lr_ub = nanx();
};

struct sym_merge_t
{
  vec<dpoint_t> rows; // ascending by d, all-NaN ranks removed
  uint64_t n_na = 0;
};

sym_merge_t sym_merge(vec<dpoint_t> ab, vec<dpoint_t> ba);

struct sym_est_t
{
  double distance = nanx();         // reported estimate (filtered or unfiltered mean)
  double median = nanx();           // median over the rows behind `distance`
  uint64_t num_filtered = 0;        // rows rejected by the lr_ub filter
  double alternative_mean = nanx(); // the mean that was not reported
  uint64_t num_na = 0;
  double max_unfiltered = nanx(); // max d over all non-NA rows
  double max_distance = nanx();   // max d over kept rows only
  uint64_t n_lr_zero = 0;
  uint64_t n_total = 0; // non-NA reconciled rows
  uint64_t n_kept = 0;  // rows passing the filter
  bool used_filtered = false;
  vec<double> null_d_v; // distances behind `distance`, ascending; detect's null
};

sym_est_t sym_estimate(const sym_merge_t& rc, double lr_th, double min_portion);

constexpr double lr_th_default = 3.841;
constexpr double min_portion_default = 0.66;

#endif
