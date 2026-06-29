#ifndef _DIM_HPP
#define _DIM_HPP

#include <simde/x86/avx512.h>
#include "llh.hpp"
#include "maptils.hpp"
#include "types.hpp"

template<typename T>
class DIM
{
  static constexpr size_t WIDTH = std::is_same_v<T, double> ? 1 : RWIDTH;

public:
  DIM(const params_t<T>& params, const llh_sptr_t<T>& llhf, uint64_t nbins, uint64_t nmers);
  void release_accumulators() noexcept;
  void inclusive_scan();
  void extrema_scan();
  void compute_prefhistsum();
  // void skip_mer(uint64_t i); // TODO: Anything better than ignoring?
  void aggregate_mer(uint32_t hdist_min, uint64_t i);
  void set_query_distance(double d_q);
  void extract_intervals_mx(uint64_t tau, uint64_t lix, uint64_t rix, size_t ix = 0);
  void extract_intervals_sx(uint64_t tau, uint64_t lix, uint64_t rix, size_t ix = 0);
  void expand_intervals(double chisq_th, size_t ix = 0);
  void total_histogram(vec<uint64_t>& v, uint64_t& u, uint64_t& t) const;
  void extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const;
  uint64_t get_nbins() const { return nbins; }
  uint64_t get_nmers() const { return nmers; }
  const vec<interval_t>& get_intervals(size_t ti) const { return intervals_v[ti]; }
  const vec<size_t>& get_thrank() const { return thrank_v; }
  static inline void add_to(T& dest, const T& source)
  {
    if constexpr (std::is_same_v<T, double>) {
      dest += source;
    } else {
      simde__m512d vd = simde_mm512_loadu_pd(dest.data());
      simde__m512d vs = simde_mm512_loadu_pd(source.data());
      vd = simde_mm512_add_pd(vd, vs);
      simde_mm512_storeu_pd(dest.data(), vd);
    }
  }

private:
  const params_t<T>& params;
  const llh_sptr_t<T> llhf;   // log-likelihood function for all calculations
  const uint64_t nbins;       // number of bins
  const uint64_t nmers;       // number of k-mers in query (for per-k-mer HD tracking)
  const bool keep_hist;       // whether to keep track of the histogram(s) for the query sequence
  uint64_t t_q = 0;           // total number of k-mers hits below hdist_th per query sequence
  uint64_t u_q = 0;           // total number misses per query sequence
  vec<uint64_t> hdisthist_v;  // D[i][j] is the number of hits with HD=j, [(nbins+1) x (hdist_th+1)] row-major; D[0][j]=0
  vec<T> fdc_v;               // The f' contribution c_i of the k-mer (bin) starting at i
  vec<T> sdc_v;               // The f'' contribution s_i of the k-mer (bin) starting at i
  vec<T> fdps_v;              // C[i] = sum(c_0, ..., c_{i}), C[0] = 0 (length n) (shifted by 1 w.r.t. fdc_v)
  vec<T> sdps_v;              // S[i] = sum(s_0, ..., s_{i}), S[0] = 0 (length n) (shifted by 1 w.r.t. sdc_v)
  vec<T> fdpmax_v;            // H[i] = max(C_1, ..., C_{i}), H_0 = -inf, H_{n+1} = inf (length n+1)
  vec<T> fdsmin_v;            // L[i] = min(C_{i}, ..., C_n), L_0 = inf, L_{n+1}= -inf (length n+1)
  vec<size_t> thrank_v;       // extraction order for this strand (depends on d_q)
  arr<bool, WIDTH> thneg_v{}; // which thresholds need a sign flip
  arr<vec<interval_t>, WIDTH> intervals_v; // 1-based inclusive bin coordinates per threshold

  void apply_threshold_signs();
};

#endif
