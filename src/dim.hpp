#ifndef _DIM_HPP
#define _DIM_HPP

#include <simde/x86/avx512.h>
#include "llh.hpp"
#include "stils.hpp"
#include "tpool.hpp"
#include "types.hpp"

class HDHist
{
public:
  HDHist() = default;
  explicit HDHist(uint64_t nbins, uint32_t hdist_th, uint64_t bin_shift);

  template<bool Atomic = false>
  void aggregate_mer(uint32_t hdist_min, uint64_t i);
  void compute_prefhistsum();
  void compute_prefhistsum_parallel(ThreadPool& pool, uint32_t nchunks);
  void extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const;
  void extract_histogram(uint64_t a, uint64_t b, window_counts_t& wc) const;
  void extract_histogram(uint64_t a, uint64_t b, swindow_counts_t& wc, bool is_rc) const;
  [[nodiscard]] uint64_t get_nbins() const { return nbins; }

private:
  uint64_t nbins = 0;
  uint32_t hdist_th = 0;
  uint64_t bin_shift = 0;
  vec<uint64_t> hist_v;
  vec<uint64_t> miss_v;
};

template<typename T>
inline double at(T v, const size_t ix)
{
  if constexpr (std::is_same_v<T, double>) {
    return v;
  } else {
    return v[ix];
  }
}

template<typename T>
class DIM
{
  static constexpr size_t WIDTH = std::is_same_v<T, double> ? 1 : RWIDTH;

public:
  DIM(const params_t<T>& params, const llh_sptr_t<T>& llhf, uint64_t nbins, uint64_t nmers);
  void inclusive_scan();
  void extrema_scan();
  void compute_prefhistsum();
  void skip_mer(uint64_t i);
  void aggregate_mer(uint32_t hdist_min, uint64_t i);
  void set_query_distance(double d_q);
  void extract_intervals_mx(uint64_t tau, uint64_t lix, uint64_t rix, size_t ix = 0);
  void extract_intervals_sx(uint64_t tau, uint64_t lix, uint64_t rix, size_t ix = 0);
  void expand_intervals(double chisq_th, size_t ix = 0);
  void total_histogram(vec<uint64_t>& v, uint64_t& u, uint64_t& t) const;
  void extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const;
  void extract_histogram(uint64_t a, uint64_t b, window_counts_t& wc) const;
  void extract_histogram(uint64_t a, uint64_t b, swindow_counts_t& wc, bool is_rc) const;
  uint64_t get_nbins() const { return nbins; }
  uint64_t get_nmers() const { return nmers; }
  uint64_t get_bin_shift() const { return params.bin_shift; }
  uint64_t get_bin_size() const { return params.bin_size; }
  bool get_has_skips() const { return has_skips; }
  bool is_skip(uint64_t i) const { return has_skips && skip_v[i - 1]; } // 1-based bin
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
  const llh_sptr_t<T> llhf;                // log-likelihood function for all calculations
  const uint64_t nbins;                    // number of bins
  const uint64_t nmers;                    // number of k-mers in query (for per-k-mer HD tracking)
  const bool keep_hist;                    // whether to keep track of the histogram(s) for the query sequence
  uint64_t t_q = 0;                        // total number of k-mers hits below hdist_th per query sequence
  uint64_t u_q = 0;                        // total number misses per query sequence
  HDHist hdhist;                           // per-bin histogram (when keep_hist); delegates aggregate/prefix/extract
  vec<uint64_t> hist_v;                    // flat query-wide histogram (when !keep_hist)
  vec<T> fdc_v;                            // The f' contribution c_i of the k-mer (bin) starting at i
  vec<T> sdc_v;                            // The f'' contribution s_i of the k-mer (bin) starting at i
  vec<T> fdps_v;                           // C[i] = sum(c_0, ..., c_{i}), C[0] = 0 (length n) (shifted by 1 w.r.t. fdc_v)
  vec<T> sdps_v;                           // S[i] = sum(s_0, ..., s_{i}), S[0] = 0 (length n) (shifted by 1 w.r.t. sdc_v)
  vec<T> fdpmax_v;                         // H[i] = max(C_1, ..., C_{i}), H_0 = -inf, H_{n+1} = inf (length n+1)
  vec<T> fdsmin_v;                         // L[i] = min(C_{i}, ..., C_n), L_0 = inf, L_{n+1}= -inf (length n+1)
  vec<size_t> thrank_v;                    // extraction order for this strand (depends on d_q)
  arr<bool, WIDTH> thneg_v{};              // which thresholds need a sign flip
  arr<vec<interval_t>, WIDTH> intervals_v; // 1-based inclusive bin coordinates per threshold
  vec<uint8_t> skip_v;                     // 0-based per-bin N-run break flags (lazy)
  bool has_skips = false;                  // whether any bin was flagged by skip_mer

  void extract_mx(uint64_t tau, uint64_t lix, uint64_t rix, size_t ix);
  void extract_sx(uint64_t tau, uint64_t lix, uint64_t rix, size_t ix);
  void apply_threshold_signs();
};

#endif
