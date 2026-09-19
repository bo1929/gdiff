#ifndef _MAP_HPP
#define _MAP_HPP

#include "CLI11.hpp"
#include "dim.hpp"
#include "records.hpp"
#include "windows.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

template<typename T>
class QIE
{
  static constexpr size_t WIDTH = std::is_same_v<T, double> ? 1 : rwidth;

public:
  QIE(const dim_params<T>& params, const Sketch& sketch, const vec<qseq_t>& batch_v);
  void map_sequences(std::ostream& sout, const str& rname);

  uint64_t get_nunmapped() const { return nunmapped; }

private:
  void sample_background(DIM<T>& dim, size_t first_record);
  void extract_ordered_intervals(DIM<T>& dim, bool is_rc, uint64_t tau_eff, double d_q_bg);
  void extract_simple_intervals(DIM<T>& dim, bool is_rc, uint64_t tau_eff, double d_q_bg);
  void emit_record(DIM<T>& dim, uint64_t a_bin, uint64_t b_bin, size_t th_ix, bool is_rc, double d_q_bg);
  xy_t get_distance_bin(const record_t& r, const vec<double>& th_sorted) const;
  void report_contiguous(std::ostream& sout, const str& rname) const;

  const dim_params<T>& params;
  const Sketch& sketch;
  const vec<qseq_t>& batch_v;
  const uint32_t k;
  const uint32_t hpos;
  LLH<T> llhf;
  bool skip_test;
  bool keep_hist;
  bool enum_only;
  bool coordinates_only;
  uint64_t enmers;
  uint64_t nbins;
  uint64_t bix;
  vec<bp_t> bp_v;
  vec<sample_t> samples_v;
  vec<record_t> records_v;
  vec<uint64_t> scratch_v;
  vec<uint64_t> acc_v;
  uint64_t u_acc = 0;
  double d_acc = nanx();
  uint64_t nunmapped = 0;
};

class MapSC
{
public:
  struct params
  {
    uint32_t hdist_th = 3;
    uint64_t tau = 1;
    uint64_t bin_shift = 0;
    double chisq = 33.00051; // chi-square(1) survival ~1e-8
    uint64_t sample_size = 200;
    bool enum_only = false;
    vec<double> thresholds_v;
  };

  MapSC(CLI::App& sc);
  void map();
  bool validate_configuration();
  uint64_t get_total_qseq() const { return total_qseq; }

private:
  str target_path;
  std::filesystem::path sketch_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  params params;
  uint64_t total_qseq = 0;
};

#endif
