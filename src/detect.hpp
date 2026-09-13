#ifndef _DETECT_HPP
#define _DETECT_HPP

#include "CLI11.hpp"
#include "gamma.hpp"
#include "llh.hpp"
#include "sym.hpp"
#include "stils.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include "tpool.hpp"
#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

class DistanceSampler;

struct clvl_t
{
  double alpha;
  double t_low;
  double t_high;
  bool high = true; // false when high-side detection is disabled for this level
  bool low = true;  // false when low-side detection is disabled for this level
};

// Result of fitting Gamma(shape, scale) to background window distances.
struct bggamma_t
{
  GammaModel::params_t params{1.0, 1.0};
  double objective = nanx();
  int niter = 0;
  uint64_t nsamples = 0;
  uint64_t nfloored = 0;
  uint64_t ndropped = 0;
  bool ok = false;
};

struct thcfg_t
{
  vec<clvl_t> levels; // alpha descending
  arr<double, RWIDTH> extrema{};
  vec<double> high_v; // ascending t_high
  vec<double> low_v;  // ascending t_low

  // Background window closest to d_median; used for the high-side screen and lr_bg.
  arr<uint64_t, hdist_bound + 1> med_hist{};
  uint64_t med_u = 0;
  double d_median = nanx();
  bool med_valid = false;

  [[nodiscard]] size_t nlevels() const { return levels.size(); }
  [[nodiscard]] size_t nlanes() const { return 2 * levels.size(); }
  [[nodiscard]] bool empty() const { return levels.empty(); }

  [[nodiscard]] bool is_high_side(size_t ix) const { return ix < nlevels(); }
  [[nodiscard]] double alpha(size_t ix) const
  {
    assert(ix < nlanes());
    return levels[is_high_side(ix) ? ix : ix - nlevels()].alpha;
  }

  [[nodiscard]] bool high_enabled(size_t level_ix) const
  {
    assert(level_ix < nlevels());
    return levels[level_ix].high;
  }

  [[nodiscard]] bool low_enabled(size_t level_ix) const
  {
    assert(level_ix < nlevels());
    return levels[level_ix].low;
  }

  [[nodiscard]] bool lane_enabled(size_t ix) const
  {
    assert(ix < nlanes());
    return is_high_side(ix) ? high_enabled(ix) : low_enabled(ix - nlevels());
  }

  void set_median_window(const uint64_t* hist, uint64_t u, double d_med)
  {
    d_median = d_med;
    med_valid = hist != nullptr && u > 0 && is_valid_distance(d_med);
    med_u = med_valid ? u : 0;
    if (med_valid)
      std::copy(hist, hist + hdist_bound + 1, med_hist.begin());
    else
      med_hist.fill(0);
  }

  void pack()
  {
    assert(nlanes() <= RWIDTH);
    extrema.fill(d_eps);
    high_v.clear();
    low_v.clear();
    high_v.reserve(nlevels());
    low_v.reserve(nlevels());
    for (size_t j = 0; j < nlevels(); ++j) {
      // Disabled lanes stay at d_eps and are skipped during extraction.
      if (levels[j].high) extrema[j] = -levels[j].t_high;
      if (levels[j].low) extrema[nlevels() + j] = levels[j].t_low;
      high_v.push_back(levels[j].t_high);
      low_v.push_back(levels[j].t_low);
    }
    std::reverse(low_v.begin(), low_v.end());
  }
};

struct lvlstat_t
{
  uint64_t nintervals = 0;
  uint64_t bp_covered = 0;
};

// Per-sketch worker: sample background, fit thresholds, extract outliers.
class Detector
{
public:
  Detector(const Sketch& sketch,
           const vec<qseq_t>& batch_v,
           uint64_t tau,
           uint64_t bin_shift,
           uint32_t hdist_th,
           double chisq,
           uint64_t sample_size,
           const vec<double>& levels,
           const vec<double>& fit_quantiles,
           bool per_sequence,
           bool low_memory,
           uint32_t verbosity);

  // Empty rev_rows: one-directional null. Otherwise reconcile both directions.
  void run(std::ostream& out, ThreadPool& pool, const vec<dpoint_t>& rev_rows, double lr_th, double min_portion);

private:
  bggamma_t fit(const vec<double>& d_v) const;
  thcfg_t thresholds_for(const bggamma_t& fit,
                         const LLH<double>& llhf,
                         const uint64_t* win_hist,
                         uint64_t win_u,
                         double d_median,
                         bool disable_high,
                         bool disable_low) const;
  // fit_v is the gamma sample (pooled or per-sequence); may differ from d_per_seq.
  vec<thcfg_t> plan(const DistanceSampler& sampler,
                    const vvec<double>& d_per_seq,
                    const vvec<double>& fit_v,
                    const LLH<double>& llhf) const;
  void extract_batch(const vec<thcfg_t>& sets,
                     bool per_sequence,
                     const LLH<double>& llhf,
                     vec<lvlstat_t>& stats,
                     uint64_t& unmapped_iv,
                     uint64_t& unmapped_bp,
                     strstream& sout,
                     ThreadPool& pool) const;
  void
  report_fit(const bggamma_t& fit, const thcfg_t& thresholds, double mean, double sd, uint64_t nunmapped, uint64_t nmapped)
    const;
  void
  report_stats(const thcfg_t& thresholds, const vec<lvlstat_t>& stats, uint64_t unmapped_iv, uint64_t unmapped_bp) const;

  const Sketch& sketch;
  const vec<qseq_t>& batch_v;
  const uint64_t tau;
  const uint64_t bin_shift;
  const uint32_t hdist_th;
  const double chisq;
  const uint64_t sample_size;
  const vec<double>& levels;
  const vec<double>& fit_quantiles;
  const bool per_sequence;
  const bool low_memory;
  const uint32_t verbosity;
  double sym_distance = nanx(); // reconciled point estimate, reported with the fit
};

class DetectSC
{
public:
  explicit DetectSC(CLI::App& sc);
  void detect();
  bool validate_configuration();

private:
  str target_path;
  std::filesystem::path sketch_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  uint64_t tau = 0;
  uint64_t sample_size = 1000;
  uint64_t bin_shift = 0;
  uint32_t hdist_th = 4;
  double chisq = 33.00051; // chi-square(1) survival ~1e-8
  vec<double> levels;
  vec<double> fit_quantiles;
  bool per_sequence = false;
  bool symmetric = true;
  double lr_th = lr_th_default;
  double min_portion = min_portion_default;
  uint64_t batch_bases = 0; // 0 = load the whole query, preserving the pooled fit
  bool low_memory = false;
  uint32_t verbosity = 1;
};

#endif
