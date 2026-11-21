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
    error_msg("--bin-shift must be in [0, 16]; got ", bin_shift);
  }
  const uint64_t bin_size = (bin_shift <= 16) ? (uint64_t(1) << bin_shift) : 0;
  if (tau && bin_size > tau) {
    is_invalid = true;
    error_msg("--bin-shift gives a bin size of ", bin_size, ", which exceeds -l ", tau);
  }
  return !is_invalid;
}

// One map pass: thresholds, background sampling and the detection mode.
struct map_opts
{
  uint32_t hdist_th = 3;
  uint64_t tau = 1;
  uint64_t bin_shift = 1;
  double chisq = 33.00051;
  uint64_t sample_size = 1000;
  vec<double> levels_v; // one to three tail probabilities pi (a lower and an upper threshold each)
  // --thresholds: explicit distance thresholds, one lane each. When non-empty, map_sequences
  // uses them as-is and skips the background detection that places thresholds for levels_v.
  vec<double> thresholds_v;
  // Low and high thresholds are always independent. By default, levels of the same side are
  // claimed in order of extremity and later levels only fill the gaps left by earlier ones.
  // The independent mode lets every level enumerate its own intervals.
  bool independent_levels = false;
};

// A background pool smaller than this cannot support a significance verdict; `--sample-size` is
// either 0 (skip the test) or at least this large, and a record is tested against at least a tenth
// of the requested pool.
inline constexpr uint64_t min_sample_size = 100;
inline constexpr uint64_t sample_size_reqf = 10; // denominator for the positive-tail fit minimum

// The background fit needs at least min_fit_windows windows above the floor spike whatever the
// pool size (the reference is skipped below it), and thresholds from fewer than warn_fit_windows
// are reported as unstable.
inline constexpr size_t min_fit_windows = 50;
inline constexpr size_t warn_fit_windows = 100;

// `--levels` budgets the lanes: two thresholds per level, so at most rwidth / 2 levels fit the
// register; the CLI caps the list at max_levels.
inline constexpr size_t max_levels = (rwidth / 2) - 1;

// `--thresholds` places one lane per value, so at most rwidth thresholds fit the register.
inline constexpr size_t max_thresholds = rwidth;

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
  BackgroundSampler(const Sketch& sketch,
                    const vec<qseq_t>& batch_v,
                    uint64_t tau,
                    uint64_t bin_shift,
                    uint32_t hdist_th,
                    const vec<double>& spurious_f_v = {});

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

// The background distance distribution: a point mass p0 at `d_floor` plus a gamma above it.
struct gamma_fit_t
{
  double shape = nanx();
  double scale = nanx();
  double p0 = 0.0;
  double d_floor = nanx();
  bool via_ml = false;
};

// Method-of-moments gamma fit; false when the sample is smaller than `min_samples` or degenerate
// (zero variance). It is the start and the fallback of `fit_background_gamma`.
bool fit_gamma_mom(const vec<double>& sample_v, gamma_fit_t& fit, size_t min_samples);

// Hurdle fit of a sorted background sample: the spike at the floor plus a gamma above it, fitted
// by ML from the moment start. False when fewer than `min_samples` windows sit above the spike or
// the continuous part is degenerate.
bool fit_background_gamma(const vec<double>& sample_v, gamma_fit_t& fit, size_t min_samples, double d_floor = nanx());

// A threshold that was asked for and not placed: the quantile that fell out of range (or the
// user's distance), the level pi that asked for it (NaN under --thresholds), its tail, and why.
struct lane_drop_t
{
  double distance = nanx();
  double pi = nanx();
  bool low = false;
  str reason;
};

// What `level_thresholds` has to say beyond its verdict: why it returned false, and every
// threshold it dropped on the way.
struct level_report_t
{
  str skip;
  vec<lane_drop_t> drops_v;
};

// The --levels threshold set from a sorted null sample. Each level pi gives a lower threshold at
// the pi quantile and an upper one at the 1 - pi quantile of the hurdle fit, with the pi that
// placed each distinct threshold (a repeated threshold keeps the smallest pi). Thresholds inside
// the floor spike, at the estimable floor or outside the distance domain are dropped and listed
// in `report`; false (with `report->skip` set) when no threshold survives or the fit is rejected.
bool level_thresholds(const vec<double>& sample_v,
                      const vec<double>& levels_v,
                      vec<double>& th_v,
                      vec<double>& pi_v,
                      size_t min_samples,
                      gamma_fit_t* fit_out,
                      double d_floor = nanx(),
                      level_report_t* report = nullptr);

// Two-sided p-value and fold from the fitted background gamma. Reads the distance off the fitted
// CDF and the fold off the fitted median, so thresholds and significance come from one model.
// Leaves the row untested when the distance is invalid or the fit is unusable.
void apply_gamma_significance(record_t& r, const gamma_fit_t& fit);

// What the sampled windows say about one reference: the null sample, the hurdle fit and the
// thresholds placed from it. Filled once per reference, then read-only during the scan.
struct background_t
{
  vec<double> null_v;     // sampled background distances, sorted
  gamma_fit_t gamma_fit;  // hurdle fit behind th_v (NaN only when --thresholds mode rejects it)
  vec<double> th_v;       // thresholds, sorted ascending and distinct; lane i carries th_v[i]
  vec<double> th_pi_v;    // the level pi that placed th_v[i] (NaN under --thresholds)
  double d_null = nanx(); // median of the null sample (the background for delta and p_interval)
  double d_min = nanx();  // smallest null window (the d_min output column)
  double d_max = nanx();  // largest null window (the d_max output column)
  // MLE of an all-matching window: the search boundary, not an interior optimum.
  double d_floor = nanx();
  // MLE of the weakest detectable window: one hit at `hdist_th`, the rest misses.
  double d_ceiling = nanx();
};

// Per-reference detection pass: null, thresholds, scan, intervals, significance.
class IntMap
{
  using params_t = map_params<cmlane_t>;
  using llh_t = LLH<cmlane_t>;
  using ext_t = IntExt<cmlane_t>;

public:
  IntMap(const map_opts& opts, const Sketch& sketch, const vec<qseq_t>& batch_v);
  void map_sequences(std::ostream& sout, const str& rname, ThreadPool& pool);
  uint64_t get_nunmapped() const { return nunmapped; }

private:
  void build_null(ThreadPool& pool);
  bool place_thresholds(std::ostream& sout, const str& rname);
  void scan_sequence(const params_t& params, const llh_t& llhf, size_t bix);
  void extract_intervals(ext_t& ext, const llh_t& llhf, bool is_rc, uint64_t tau_eff, size_t bix);
  void extract_lane_intervals(ext_t& ext, const llh_t& llhf, bool is_rc, uint64_t tau_eff, size_t bix);
  void
  extract_iterative_intervals(ext_t& ext, const llh_t& llhf, bool is_rc, uint64_t tau_eff, size_t bix, bool low_side);
  // `a_bin`/`b_bin` are the half-open 0-based bin range the extractor scored: bins [a_bin, b_bin).
  void emit_record(ext_t& ext, const llh_t& llhf, size_t bix, uint64_t a_bin, uint64_t b_bin, size_t th_ix, bool is_rc);
  void report_skip(std::ostream& sout, const str& rname, const str& why) const;
  void report_lanes(std::ostream& sout, const str& rname, const vec<lane_drop_t>& drops_v) const;
  void warn_unresolvable_lanes(const str& rname) const;
  void report_null(const str& rname) const;
  void report_calls(std::ostream& sout, const str& rname) const;
  void write_call(std::ostream& sout, const str& rname, const record_t& r) const;

  const map_opts& opts;
  const Sketch& sketch;
  const vec<qseq_t>& batch_v;
  const uint32_t k;
  const uint32_t hpos;
  const bool sketch_canonical;
  vec<double> spurious_f_v; // experimental non-self hit distribution
  background_t bg;
  vec<uint64_t> scratch_v;
  vec<bp_t> bp_v;
  vec<record_t> records_v;
  uint64_t nunmapped = 0;
};

class MapSC
{
public:
  explicit MapSC(CLI::App& sc);
  void map();
  bool validate_configuration();

private:
  str query_path;
  std::filesystem::path sketch_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  map_opts params;
};

#endif
