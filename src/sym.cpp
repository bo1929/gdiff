#include "sym.hpp"

#include <algorithm>
#include <cmath>

namespace
{

  // Ascending by d, NaN last.
  bool d_nan_last(const dpoint_t& lhs, const dpoint_t& rhs)
  {
    const bool lnan = std::isnan(lhs.d);
    const bool rnan = std::isnan(rhs.d);
    if (lnan || rnan) return !lnan && rnan;
    return lhs.d < rhs.d;
  }

} // namespace

sym_merge_t sym_merge(vec<dpoint_t> ab, vec<dpoint_t> ba)
{
  std::sort(ab.begin(), ab.end(), d_nan_last);
  std::sort(ba.begin(), ba.end(), d_nan_last);

  sym_merge_t rc;
  const size_t nmax = std::max(ab.size(), ba.size());
  rc.rows.reserve(nmax);
  for (size_t i = 0; i < nmax; ++i) {
    const bool has_ab = i < ab.size() && !std::isnan(ab[i].d);
    const bool has_ba = i < ba.size() && !std::isnan(ba[i].d);
    if (has_ab && has_ba) {
      rc.rows.push_back(ab[i].d <= ba[i].d ? ab[i] : ba[i]);
    } else if (has_ab) {
      rc.rows.push_back(ab[i]);
    } else if (has_ba) {
      rc.rows.push_back(ba[i]);
    } else {
      ++rc.n_na;
    }
  }
  return rc;
}

sym_est_t sym_estimate(const sym_merge_t& rc, const double lr_th, const double min_portion)
{
  sym_est_t est;
  est.num_na = rc.n_na;

  double unf_sum = 0.0, fil_sum = 0.0;
  vec<double> unf_d, fil_d;
  unf_d.reserve(rc.rows.size());
  fil_d.reserve(rc.rows.size());

  for (const dpoint_t& row : rc.rows) {
    if (std::isnan(row.d)) {
      ++est.num_na;
      continue;
    }
    if (row.lr_ub == 0.0) ++est.n_lr_zero;

    unf_sum += row.d;
    unf_d.push_back(row.d);
    if (std::isnan(est.max_unfiltered) || row.d > est.max_unfiltered) est.max_unfiltered = row.d;

    // Missing lr_ub is rejected.
    if (!std::isnan(row.lr_ub) && row.lr_ub > lr_th) {
      fil_sum += row.d;
      fil_d.push_back(row.d);
      if (std::isnan(est.max_distance) || row.d > est.max_distance) est.max_distance = row.d;
    } else {
      ++est.num_filtered;
    }
  }

  est.n_total = unf_d.size();
  est.n_kept = fil_d.size();
  const double portion =
    est.n_total ? static_cast<double>(est.n_kept) / static_cast<double>(est.n_total) : 0.0;
  const double unf_mean = est.n_total ? unf_sum / static_cast<double>(est.n_total) : nanx();
  const double fil_mean = est.n_kept ? fil_sum / static_cast<double>(est.n_kept) : nanx();

  est.used_filtered = portion > min_portion;
  if (est.used_filtered) {
    est.distance = fil_mean;
    est.alternative_mean = unf_mean;
    est.null_d_v = std::move(fil_d);
  } else {
    est.distance = unf_mean;
    est.alternative_mean = fil_mean;
    est.null_d_v = std::move(unf_d);
  }
  // rc.rows is already ascending, so the retained subsequence is too.
  est.median = linear_quantile(est.null_d_v, 0.5);
  return est;
}
