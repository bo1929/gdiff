#ifndef _MAP_HPP
#define _MAP_HPP

#include "CLI11.hpp"
#include "dim.hpp"
#include "diststat.hpp"
#include "lshf.hpp"
#include "maptils.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

template<typename T>
class QIE
{
  static constexpr size_t WIDTH = std::is_same_v<T, double> ? 1 : RWIDTH;

public:
  QIE(const params_t<T>& params,
      const sketch_sptr_t& sketch,
      const lshf_sptr_t& lshf,
      const vec<str>& seq_batch,
      const vec<str>& qid_batch);
  void map_sequences(std::ostream& sout, const str& rid);

  uint64_t get_n_unmapped() const { return n_unmapped; }                         // emitted intervals with no k-mer hits
  uint64_t get_n_unmapped_samples() const { return diststat->get_n_unmapped(); } // null windows with no k-mer hits

private:
  void search_mers(const char* cseq, uint64_t len, DIM<T>& dim);
  void search_mers(const char* cseq, uint64_t len, DIM<T>& dim_fw, DIM<T>& dim_rc);
  void extract_ordered_intervals(DIM<T>& dim, bool is_rc, uint64_t tau_eff, double d_q_bg);
  void extract_simple_intervals(DIM<T>& dim, bool is_rc, uint64_t tau_eff, double d_q_bg);
  void emit_record(DIM<T>& dim, uint64_t a_bin, uint64_t b_bin, size_t th_ix, bool is_rc, double d_q_bg);
  xy_t get_distance_bin(const record_t& r, const arr<double, WIDTH>& th_sorted) const;
  void report_contiguous(std::ostream& sout, const str& rid) const;

  const params_t<T>& params;
  const sketch_sptr_t sketch;
  const lshf_sptr_t lshf;
  const vec<str>& seq_batch;
  const vec<str>& qid_batch;
  const uint64_t batch_size;
  const uint32_t k;
  const uint32_t h;
  llh_sptr_t<T> llhf;
  bool skip_test;
  bool keep_hist;
  bool enum_only;
  bool coordinates_only;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint64_t onmers; // Number of observed k-mers in current query (e.g., due to Ns)
  uint64_t enmers; // Number of expected k-mers in current query (= len - k + 1)
  uint64_t nbins;  // Number of bins: ceil(enmers / bin_len)
  uint64_t bix;    // Index of the current query in the this batch
  vec<bp_t> bp_v;  // Breakpoint scratch reused across intervals
  diststat_sptr_t<T> diststat;
  vec<record_t> records_v;
  vec<uint64_t> v_scratch;
  vec<uint64_t> v_acc;
  uint64_t u_acc = 0;
  double d_acc = nanx();
  uint64_t n_unmapped = 0; // emitted intervals with no k-mer hits (t == 0)
};

class MapSC
{
public:
  MapSC(CLI::App& sc);
  void map();
  bool validate_configuration();
  uint64_t get_total_qseq() const { return total_qseq; }

private:
  str query_path;
  std::filesystem::path sketch_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  uint32_t hdist_th = 4;
  uint64_t tau = 1;
  uint64_t bin_shift = 0;
  double chisq = 33.00051; // 1e-10
  uint64_t sample_size = 200;
  bool enum_only = false;
  std::vector<double> dist_th;
  uint64_t total_qseq = 0;
};

#endif
