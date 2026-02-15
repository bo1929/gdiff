#ifndef _MAP_H
#define _MAP_H

#include <simde/x86/avx512.h>
#include "llh.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include "hm.hpp"
#include "enc.hpp"
#include "types.hpp"
#include "exthash.hpp"

template<typename T>
class DIM
{
  static constexpr size_t WIDTH = std::is_same_v<T, double> ? 1 : RWIDTH;

public:
  DIM(llh_sptr_t<T> llhf, uint32_t hdist_th, uint64_t en_mers);
  T get_fdt() const { return fdt; }
  T get_sdt() const { return sdt; }
  static inline double at(T v, size_t idx);
  T fdt_at(uint64_t i) const { return fdc_v[i]; }
  T sdt_at(uint64_t i) const { return sdc_v[i]; }
  // void optimize_loglikelihood();
  void inclusive_scan();
  // void skip_mer(uint64_t i); // TODO: Anything better?
  void aggregate_mer(sketch_sptr_t sketch, uint32_t rix, enc_t enc_lr, uint64_t i);
  void extract_intervals(uint64_t tau, size_t idx = 0);
  uint64_t expand_intervals(double chisq_th, size_t idx = 0);
  interval_t get_interval(uint64_t i, size_t idx = 0);

private:
  const llh_sptr_t<T> llhf;
  const uint64_t en_mers;
  const uint32_t hdist_th;
  uint64_t merhit_count = 0;
  uint64_t merna_count = 0;
  uint64_t mermiss_count = 0;
  vec<uint64_t> hdisthist_v;
  T fdt;
  T sdt;
  vec<T> fdc_v;
  vec<T> sdc_v;
  vec<T> fdps_v;
  vec<T> sdps_v;
  vec<T> fdpmax_v;
  vec<T> fdsmin_v;
  arr<vec<interval_t>, WIDTH> rintervals_v;
  arr<vec<interval_t>, WIDTH> eintervals_v;
  arr<vec<double>, WIDTH> chisq_v;
  double d_llh = std::numeric_limits<double>::quiet_NaN();
  double v_llh = std::numeric_limits<double>::quiet_NaN();
};

template<typename T>
class QIE
{
  static constexpr size_t WIDTH = std::is_same_v<T, double> ? 1 : RWIDTH;

public:
  QIE(sketch_sptr_t sketch, lshf_sptr_t lshf, qseq_sptr_t qs, params_t<T> params);
  void map_sequences(std::ostream& output_stream);

private:
  void search_mers(const char* cseq, uint64_t len, DIM<T>& or_summary, DIM<T>& rc_summary);
  void report_intervals(std::ostream& output_stream, DIM<T>& dim, size_t idx = 0);
  static inline double at(T v, size_t idx);

  const sketch_sptr_t sketch;
  const lshf_sptr_t lshf;
  const uint64_t batch_size;
  const uint32_t k;
  const uint32_t h;
  const uint32_t m;
  const params_t<T> params;
  const uint64_t min_length;
  const double chisq; // 3.841; // 95%
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint64_t onmers;
  uint64_t en_mers;
  uint64_t bix;
  llh_sptr_t<T> llhf;
  vec<std::string> seq_batch;
  vec<std::string> identifier_batch;
};

#define WRITE_CINTERVAL(identifier, a, b, n, dist_th) identifier << '\t' << a << '\t' << b << '\t' << n << '\t' << dist_th

#endif
