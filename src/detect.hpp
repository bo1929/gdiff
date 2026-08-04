#ifndef _DETECT_HPP
#define _DETECT_HPP

#include "CLI11.hpp"
#include "dim.hpp"
#include "dist.hpp"
#include "gamma.hpp"
#include "llh.hpp"
#include "maptils.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include "tpool.hpp"
#include "types.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

struct detect_threshold_t
{
  double alpha; // two-sided confidence level
  double t_lo;  // Q(alpha/2): extract regions with local distance < t_lo
  double t_hi;  // Q(1-alpha/2): extract regions with local distance > t_hi
};

struct gamma_fit_t
{
  GammaModel::params_t params{1.0, 1.0}; // shape, scale
  double objective = nanx();             // final SSE over the fit quantiles
  int niter = 0;                         // Nelder-Mead iterations used
  uint64_t n_samples = 0;                // finite samples used in the fit
  uint64_t n_zeros = 0;                  // samples floored to d_eps
  uint64_t n_dropped = 0;                // non-finite samples dropped
  bool ok = false;
};

// SIMD lane layout for detection: [-t_hi per level (alpha descending),
// +t_lo per level (alpha descending)]. The LLH contributions are the
// log-likelihood derivative f'(t), so prefix-sum drops correspond to local
// d < t with positive sign and to local d > t with negative sign: negative
// extrema mark high-side (too diverged) lanes, positive extrema low-side
// (too similar) lanes. Unused lanes are filled with 0.5 and never extracted.
struct lane_config_t
{
  size_t n_levels = 0;
  size_t n_lanes = 0;
  arr<double, RWIDTH> extrema{};
  arr<double, RWIDTH> alpha_v{};
  vec<detect_threshold_t> thresholds;
};

struct detect_level_stats_t
{
  uint64_t n_intervals = 0;
  uint64_t bp_covered = 0;
};

class DetectSC
{
public:
  explicit DetectSC(CLI::App& sc);
  void detect();
  bool validate_configuration();

private:
  // Pass 1: sample tau-mer windows across the query file and compute their
  // MLE distances; d_per_seq is indexed by query batch index. Windows with no
  // k-mer hits (t == 0) are unmapped: excluded from d_per_seq, counted in
  // n_unmapped.
  void sample_null(const sketch_sptr_t& sketch,
                   const vec<str>& seq_batch,
                   vec<vec<double>>& d_per_seq,
                   uint64_t& n_unmapped,
                   ThreadPool& pool) const;

  gamma_fit_t fit_null(const vec<double>& d_v) const;
  bool build_lanes(const gamma_fit_t& fit, lane_config_t& lanes) const;

  // Pass 2: scan every query and extract maximal outlier intervals per lane.
  // Intervals with no k-mer hits are emitted as side="unmapped" (d = NaN) and
  // counted in unmapped_iv / unmapped_bp instead of the per-lane stats.
  void detect_queries(const sketch_sptr_t& sketch,
                      const vec<str>& seq_batch,
                      const vec<str>& qid_batch,
                      const lane_config_t* lanes_pooled,
                      const vec<lane_config_t>* lanes_per_seq,
                      vec<detect_level_stats_t>& stats,
                      uint64_t& unmapped_iv,
                      uint64_t& unmapped_bp,
                      strstream& sout,
                      ThreadPool& pool) const;

  void report_fit(const str& rid,
                  const dist_summary_t& summary,
                  const gamma_fit_t& fit,
                  const lane_config_t& lanes,
                  uint64_t n_unmapped) const;
  void report_stats(const str& rid,
                    const lane_config_t& lanes,
                    const vec<detect_level_stats_t>& stats,
                    uint64_t unmapped_iv,
                    uint64_t unmapped_bp) const;

  str query_path;
  std::filesystem::path sketch_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  uint64_t tau = 0;
  uint64_t sample_size = 1000;
  uint64_t bin_shift = 0;
  uint32_t hdist_th = 4;
  double chisq = 33.00051; // 1e-10
  vec<double> levels;
  vec<double> fit_quantiles;
  str fit_scope = "per-sketch";
  uint32_t verbosity = 1;
};

#endif
