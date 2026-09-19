#include "sym.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace {

  bool order_finite_first(const ds_t& lhs, const ds_t& rhs)
  {
    const bool lnan = std::isnan(lhs.d);
    const bool rnan = std::isnan(rhs.d);
    if (lnan || rnan) return !lnan && rnan;
    return lhs.d < rhs.d;
  }

  std::pair<vec<ds_t>, uint64_t> merge_directions(vec<ds_t>& ab_v, vec<ds_t>& ba_v)
  {
    std::sort(ab_v.begin(), ab_v.end(), order_finite_first);
    std::sort(ba_v.begin(), ba_v.end(), order_finite_first);

    uint64_t n_na = 0;
    const size_t nranks = std::max(ab_v.size(), ba_v.size());
    vec<ds_t> ds_v;
    ds_v.reserve(nranks);
    for (size_t i = 0; i < nranks; ++i) {
      const bool has_ab = i < ab_v.size() && !std::isnan(ab_v[i].d);
      const bool has_ba = i < ba_v.size() && !std::isnan(ba_v[i].d);
      if (has_ab && has_ba) {
        ds_v.push_back(ab_v[i].d <= ba_v[i].d ? ab_v[i] : ba_v[i]);
      } else if (has_ab) {
        ds_v.push_back(ab_v[i]);
      } else if (has_ba) {
        ds_v.push_back(ba_v[i]);
      } else {
        ++n_na;
      }
    }
    return {std::move(ds_v), n_na};
  }

} // namespace

summary_t summarize_symmetric(vec<ds_t> ab_v, vec<ds_t> ba_v, double lr_th, double min_portion)
{
  summary_t s;
  const auto [ds_v, n_na] = merge_directions(ab_v, ba_v);
  s.n_na = n_na;

  double bf_sum = 0.0, af_sum = 0.0;
  vec<double> bf_d_v, af_d_v;
  bf_d_v.reserve(ds_v.size());
  af_d_v.reserve(ds_v.size());

  for (const ds_t& ds : ds_v) {
    if (ds.s == 0.0) ++s.n_ub;

    bf_sum += ds.d;
    bf_d_v.push_back(ds.d);
    if (std::isnan(s.d_highest) || ds.d > s.d_highest) s.d_highest = ds.d;

    // A missing likelihood-ratio bound carries no evidence, so it is rejected.
    if (!std::isnan(ds.s) && ds.s > lr_th) {
      af_sum += ds.d;
      af_d_v.push_back(ds.d);
      if (std::isnan(s.d_upper) || ds.d > s.d_upper) s.d_upper = ds.d;
    } else {
      ++s.n_filtered;
    }
  }

  const size_t n_total = bf_d_v.size();
  const size_t n_kept = af_d_v.size();
  const double portion = n_total ? static_cast<double>(n_kept) / static_cast<double>(n_total) : 0.0;
  const double bf_mean = n_total ? bf_sum / static_cast<double>(n_total) : nanx();
  const double af_mean = n_kept ? af_sum / static_cast<double>(n_kept) : nanx();

  // The plain mean is always reported, so a caller can compare it with the filtered estimate.
  s.d_mean = bf_mean;
  if (portion > min_portion) {
    s.d = af_mean;
    s.d_v = std::move(af_d_v);
  } else {
    s.d = bf_mean;
    s.d_v = std::move(bf_d_v);
  }
  s.d_median = linear_quantile(s.d_v, 0.5);
  return s;
}
