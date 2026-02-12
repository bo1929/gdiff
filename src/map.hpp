#ifndef _COMPARE_H
#define _COMPARE_H

#include "common.hpp"
#include "llh.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include "hm.hpp"
#include "enc.hpp"
#include "types.hpp"
#include "exthash.hpp"

class LLH;

class DIM
{
public:
  DIM(llh_sptr_t llhf, uint64_t en_mers);
  double get_fdt();
  double get_sdt();
  double fdt_at(uint64_t i);
  double sdt_at(uint64_t i);
  void inclusive_scan();
  void optimize_loglikelihood();
  void extract_intervals(uint64_t tau);
  uint64_t expand_intervals(double chisq_th);
  void report_intervals(std::ostream& output_stream, std::string& identifer);
  void aggregate_mer(sketch_sptr_t sketch, uint32_t rix, enc_t enc_lr, uint64_t i);
  void skip_mer(uint64_t i);

private:
  llh_sptr_t llhf;
  const uint64_t en_mers;
  const bool opposite;
  const uint32_t hdist_th;
  uint64_t merhit_count = 0;
  uint64_t merna_count = 0;
  uint64_t mermiss_count = 0;
  std::vector<uint64_t> hdisthist_v;
  double d_llh = std::numeric_limits<double>::quiet_NaN();
  double v_llh = std::numeric_limits<double>::quiet_NaN();
  double fdt = 0;
  double sdt = 0;
  vec<double> fdc_v;
  vec<double> sdc_v;
  vec<double> fdps_v;
  vec<double> sdps_v;
  vec<double> fdpmax_v;
  vec<double> fdsmin_v;
  vec<interval_t> rintervals_v;
  vec<interval_t> eintervals_v;
  vec<double> chisq_v;
};

class SBatch
{
public:
  SBatch(sketch_sptr_t sketch, qseq_sptr_t qs, uint32_t hdist_th, double dist_th, uint64_t min_length, double chisq);
  void map_sequences(std::ostream& output_stream);
  void search_mers(const char* seq, uint64_t len, DIM& or_summary, DIM& rc_summary);

private:
  uint32_t k;
  uint32_t h;
  uint32_t m;
  uint32_t hdist_th;
  lshf_sptr_t lshf;
  sketch_sptr_t sketch;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint64_t en_mers;
  uint64_t onmers;
  uint64_t batch_size;
  uint64_t bix;
  vec<std::string> seq_batch;
  vec<std::string> identifer_batch;
  llh_sptr_t llhf;
  double rho;
  // vec<double> s1;
  // vec<double> c1;
  // vec<double> s2;
  // vec<double> c2;
  vec<double> prefix_sum_c;
  vec<double> prefix_sum_s;
  vec<double> prefix_maxima;
  vec<double> suffix_minima;
  double mp;
  double mp2;
  double csum1 = 0;
  double ssum1 = 0;
  double csum2 = 0;
  double ssum2 = 0;
  double dist_th;
  uint64_t min_length = 0;
  double chisq = 3.841; // 95%
};

#endif
