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
  uint64_t sample_size = 500;
  vec<double> thresholds_v; // -d: exactly 1 or rwidth absolute distance thresholds
  vec<double> levels_v;     // --levels: exactly 4 two-sided tail probabilities
  bool enum_only = false;
  bool per_sequence = false;
  uint32_t verbosity = 1;
  // TODO: --lr-th / --min-portion may return once a reverse null is available (see dist.hpp).
};

// One scored background window: its source sequence, start, and distance.
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

  // Draw at most `sample_size` windows (per query sequence when `per_sequence`) and score
  // them against the reference. Returned in source order. Call once: the points are moved out.
  vec<sample_point_t> sample(uint64_t sample_size, bool per_sequence, ThreadPool& pool);

private:
  // A contiguous run of points drawn from one source sequence.
  struct batch_t
  {
    uint64_t bix = 0;
    uint64_t enmers = 0;
    uint64_t point0 = 0;
    uint64_t nsamples = 0;
  };

  void plan(uint64_t sample_size, bool per_sequence);
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

// Empirical two-sided thresholds for `levels_v` over a sorted background pool.
bool thresholds_from_levels(const vec<double>& pool_v, const vec<double>& levels_v, vec<double>& out_v);

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
  void extract_simple_intervals(IntExt<T>& ext, const LLH<T>& llhf, bool is_rc, uint64_t tau_eff, size_t bix);
  void extract_ordered_intervals(IntExt<T>& ext, const LLH<T>& llhf, bool is_rc, uint64_t tau_eff, size_t bix);
  void emit_record(IntExt<T>& ext, const LLH<T>& llhf, size_t bix, uint64_t a_bin, uint64_t b_bin, size_t th_ix, bool is_rc);
  xy_t get_distance_bin(const record_t& r, const LLH<T>& llhf) const;
  void report_null(const str& rname) const;
  void report_contiguous(std::ostream& sout, const str& rname, const LLH<T>& llhf) const;

  const map_opts& opts;
  const Sketch& sketch;
  const vec<qseq_t>& batch_v;
  const uint32_t k;
  const uint32_t hpos;
  const bool sketch_canonical;
  vec<double> th_v;          // resolved thresholds, sorted ascending
  vec<double> null_v;        // pooled background distances, sorted
  uint64_t null_floored = 0; // background samples below d_eps
  vvec<double> null_seq_v;   // per-sequence nulls (--per-sequence)
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
