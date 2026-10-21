#ifndef _MAP_HPP
#define _MAP_HPP

#include "CLI11.hpp"
#include "intext.hpp"
#include "msg.hpp"
#include "records.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include "tpool.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

// Reject a bin shift that cannot be represented, or that quantises wider than the window.
inline bool validate_binning(uint64_t bin_shift, uint64_t tau)
{
  bool is_invalid = false;
  if (bin_shift > 16) {
    is_invalid = true;
    cerr_msg("--bin-shift must be less than or equal to 16; got ", bin_shift);
  }
  const uint64_t bin_size = (bin_shift <= 16) ? (uint64_t(1) << bin_shift) : 0;
  if (tau && bin_size > tau) {
    is_invalid = true;
    cerr_msg("--bin-shift gives bin_size=", bin_size, ", which exceeds -l=", tau);
  }
  return !is_invalid;
}

// One map pass: thresholds, background sampling and the detection mode.
struct map_opts
{
  uint32_t hdist_th = 3;
  uint64_t tau = 1;
  uint64_t bin_shift = 0;
  double chisq = 33.00051;
  uint64_t sample_size = 1000;
  vec<double> levels_v; // one to three two-sided tail probabilities (two thresholds each)
  // Deprecated compatibility fields for direct IntMap callers; the CLI no longer sets them.
  vec<double> thresholds_v;
  bool iterative_removal = false;
  // Low and high thresholds are always independent. By default, levels of the same side are
  // claimed in order of extremity and later levels only fill the gaps left by earlier ones.
  // The independent mode lets every level enumerate its own intervals.
  bool independent_levels = false;
  uint32_t verbosity = 1;
};

// A background pool smaller than this cannot support a significance verdict; `--sample-size` is
// either 0 (skip the test) or at least this large, and a record is tested against at least a tenth
// of the requested pool.
inline constexpr uint64_t min_sample_size = 100;
inline constexpr uint64_t sample_size_reqf = 10; // denominator for the positive-tail fit minimum

// `--levels` budgets the lanes: two thresholds per level, so at most rwidth / 2 levels fit the
// register; the CLI caps the list at max_levels.
inline constexpr size_t max_levels = (rwidth / 2) - 1;

// One scored background window: its source sequence, start and MLE distance.
struct sample_point_t
{
  uint64_t bix = 0;
  uint64_t start = 0;
  double d = nanx();
};

// Samples tau-long windows from the query sequences and scores them against one reference,
// producing that reference's empirical background null.
class BackgroundSampler
{
public:
  BackgroundSampler(const Sketch& sketch, const vec<qseq_t>& batch_v, uint64_t tau, uint64_t bin_shift, uint32_t hdist_th);

  // Draw at most `sample_size` windows across the whole batch and score them against the
  // reference. Returned in source order. Call once: the points are moved out.
  vec<sample_point_t> sample(uint64_t sample_size, ThreadPool& pool);

private:
  // A contiguous run of points drawn from one source sequence.
  struct batch_t
  {
    uint64_t bix = 0;
    uint64_t enmers = 0;
    uint64_t point0 = 0;
    uint64_t nsamples = 0;
  };

  void plan(uint64_t sample_size);
  void evaluate(ThreadPool& pool);

  const Sketch& sketch;
  const vec<qseq_t>& batch_v;
  bool sketch_canonical;
  uint32_t hdist_th;
  uint64_t tau;
  uint64_t bin_shift;
  uint64_t nwinmers;
  LLH<double> llhf;
  vec<sample_point_t> points_v;
  vec<batch_t> batches_v;
};

// A 1-based half-open bin range plus the threshold that claimed it.
struct bp_t
{
  uint64_t a_bin = 0;
  uint64_t b_bin = 0;
  size_t ix = 0;
};

// The background distance distribution: a hurdle of a point mass p0 at the observed
// minimum d_min (windows with the same histogram get the same MLE, so that spike is exact)
// above it. `via_ml` says whether the gamma came from the likelihood or the moment fallback.
struct gamma_fit_t
{
  double shape = nanx();
  double scale = nanx();
  double p0 = 0.0;
  double d_min = nanx();
  bool via_ml = false;
};

// Method-of-moments gamma fit; false when the sample is smaller than `min_samples` or degenerate
// (zero variance). It is the start and the fallback of `fit_background_gamma`.
bool fit_gamma_mom(const vec<double>& sample_v, gamma_fit_t& fit, size_t min_samples);

// plus a gamma above it.
// by ML from the moment start. False when fewer than `min_samples` windows sit above the spike or
// the continuous part is degenerate.
bool fit_background_gamma(const vec<double>& sample_v, gamma_fit_t& fit, size_t min_samples,
                          double floor = nanx());

// The --levels threshold set from a sorted null sample, including the alpha that placed each
// distinct threshold. Thresholds at the estimable floor are skipped on either side.
bool level_thresholds(const vec<double>& sample_v,
                      const vec<double>& levels_v,
                      vec<double>& th_v,
                      vec<double>& alpha_v,
                      size_t min_samples,
                      gamma_fit_t* fit_out,
                      double floor = nanx());

// Two-sided thresholds for `levels_v` from the background pool itself: one low and one high
// threshold per level, sorted ascending. `min_samples` is the smallest pool that may place a
// quantile; equal thresholds are left to `compose_lanes` to collapse into one lane.
bool empirical_thresholds_from_levels(const vec<double>& sample_v,
                                      const vec<double>& levels_v,
                                      vec<double>& out_v,
                                      size_t min_samples);

// The `--levels` lane set: the two bounds plus the model's own thresholds, sorted and reduced to the
// distinct values a lane can carry. Nothing is clamped - a threshold outside the bounds stays a lane
// of its own - so the only threshold a lane can do without is an exact duplicate of another, which
// is one band. Returns the number of lanes; `th_v` is sorted ascending and free of duplicates.
size_t compose_lanes(vec<double>& th_v, const a_range_t& ar);

// Two-sided p-value and fold from the fitted background gamma: the model-based counterpart of
// `apply_empirical_significance`, used when `--bg-model gamma` placed the thresholds. Reads the
// distance off the fitted CDF and the fold off the fitted median, so thresholds and significance
// come from one model. Leaves the row untested when the distance is invalid or the fit is unusable.
void apply_gamma_significance(record_t& r, const gamma_fit_t& fit);

// Two-sided thresholds for `levels_v` as quantiles of the gamma fitted to `sample_v`, which must be
// sorted. A tail quantile that reaches the estimable floor is lifted to the nearest observed
// distance. When the gamma is rejected (too few samples, no variance, or quantiles that cannot be
// placed) this falls back to the empirical quantiles of the same sample; it returns false only when
// that fails too, in which case the caller skips the sketch. `fit_out` receives the fitted
// parameters on the gamma path and NaN on the empirical one, so the caller can report which was used.
bool thresholds_from_levels(const vec<double>& sample_v,
                            const vec<double>& levels_v,
                            vec<double>& out_v,
                            size_t min_samples,
                            gamma_fit_t* fit_out = nullptr,
                            double floor = nanx());

// Per-reference detection pass: null, thresholds, scan, intervals, significance.
template<typename T>
class IntMap
{
  static constexpr size_t WIDTH = std::is_same_v<T, double> ? 1 : rwidth;

public:
  IntMap(const map_opts& opts, const Sketch& sketch, const vec<qseq_t>& batch_v);
  void map_sequences(std::ostream& sout, const str& rname, ThreadPool& pool);
  uint64_t get_nunmapped() const { return nunmapped; }

private:
  void build_null(ThreadPool& pool);
  void scan_sequence(const map_params<T>& params, const LLH<T>& llhf, size_t bix);
  void extract_lane_intervals(IntExt<T>& ext, const LLH<T>& llhf, bool is_rc, uint64_t tau_eff, size_t bix);
  void extract_iterative_intervals(IntExt<T>& ext, const LLH<T>& llhf, bool is_rc, uint64_t tau_eff, size_t bix, bool low_side);
  // `a_bin`/`b_bin` are the half-open 0-based bin range the extractor scored: bins [a_bin, b_bin).
  void emit_record(IntExt<T>& ext, const LLH<T>& llhf, size_t bix, uint64_t a_bin, uint64_t b_bin, size_t th_ix, bool is_rc);
  void report_null(const str& rname) const;
  void report_contiguous(std::ostream& sout, const str& rname) const;

  const map_opts& opts;
  const Sketch& sketch;
  const vec<qseq_t>& batch_v;
  const uint32_t k;
  const uint32_t hpos;
  const bool sketch_canonical;
  vec<double> th_v;          // thresholds, sorted ascending and distinct; lanes are [0, n_active)
  vec<double> th_alpha_v;    // the level alpha that placed th_v[i]
  size_t n_active = 0;       // lanes carrying a threshold
  vec<double> null_v;        // sampled background distances, sorted
  gamma_fit_t gamma_fit;     // hurdle fit behind th_v (NaN on the empirical fallback)
  double d_null = nanx();   // median of the sampled null
  double d_min = nanx();    // minimum sampled null distance
  double d_max = nanx();    // maximum sampled null distance
  double d_floor = nanx();  // minimum estimable distance
  double d_ceiling = nanx();// maximum estimable distance
  uint64_t null_floored = 0; // background samples below d_eps
  vec<uint64_t> scratch_v;
  vec<uint64_t> acc_v;
  uint64_t u_acc = 0;
  vec<bp_t> bp_v;
  vec<record_t> records_v;
  double d_acc = nanx();
  uint64_t nunmapped = 0;
};

class MapSC
{
public:
  explicit MapSC(CLI::App& sc);
  void map();
  bool validate_configuration();
  uint64_t get_total_qseq() const { return total_qseq; }

private:
  str query_path;
  std::filesystem::path sketch_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  map_opts params;
  uint64_t total_qseq = 0;
};

#endif
