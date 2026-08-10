#ifndef _DETECT_HPP
#define _DETECT_HPP

#include "CLI11.hpp"
#include "gamma.hpp"
#include "llh.hpp"
#include "stils.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include "tpool.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>

class DetectSC
{
public:
  explicit DetectSC(CLI::App& sc);
  void detect();
  bool validate_configuration();

private:
  struct detect_threshold_t
  {
    double alpha;
    double t_lo;
    double t_hi;
  };

  struct background_fit_t
  {
    GammaModel::params_t params{1.0, 1.0};
    double objective = nanx();
    int niter = 0;
    uint64_t nsamples = 0;
    uint64_t nfloored = 0;
    uint64_t ndropped = 0;
    bool ok = false;
  };

  struct lane_config_t
  {
    size_t nlevels = 0;
    size_t nlanes = 0;
    arr<double, RWIDTH> extrema{};
    arr<double, RWIDTH> alpha_v{};
    vec<detect_threshold_t> thresholds;
  };

  struct detect_level_stats_t
  {
    uint64_t nintervals = 0;
    uint64_t bp_covered = 0;
  };

  background_fit_t fit_background(const vec<double>& d_v) const;

  // Builds per-level lanes from the fit. When win_hist/win_u describe a
  // representative background window (nullptr/0 disables the screen), levels
  // whose thresholds a typical window cannot distinguish from the fitted median
  // d_med (chi-square(1) LR test at the level's own alpha) are dropped.
  bool build_lanes(const background_fit_t& fit,
                   lane_config_t& lanes,
                   const LLH<double>& llhf,
                   const uint64_t* win_hist,
                   uint64_t win_u,
                   double d_med) const;

  void detect_queries(const sketch_sptr_t& sketch,
                      const vec<qseq_t>& batch_v,
                      const std::optional<lane_config_t>& lanes_pooled,
                      const vec<lane_config_t>* lanes_per_seq,
                      vec<detect_level_stats_t>& stats,
                      uint64_t& unmapped_iv,
                      uint64_t& unmapped_bp,
                      strstream& sout,
                      ThreadPool& pool) const;

  void report_fit(const str& rname,
                  const background_fit_t& fit,
                  const lane_config_t& lanes,
                  double mean,
                  double sd,
                  uint64_t nunmapped,
                  uint64_t nmapped) const;
  void report_stats(const str& rname,
                    const lane_config_t& lanes,
                    const vec<detect_level_stats_t>& stats,
                    uint64_t unmapped_iv,
                    uint64_t unmapped_bp) const;

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
  uint32_t verbosity = 1;
};

#endif
