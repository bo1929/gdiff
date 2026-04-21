#ifndef _MAP_HPP
#define _MAP_HPP

#include <simde/x86/avx512.h>
#include "gamma.hpp"
#include "llh.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"

static constexpr double eps = 1e-7;
static constexpr uint32_t hdist_bound = 7;
// Striding parameters for null sample grid
static constexpr double grid_growth = 1.25; // Incremental step size multiplier across the grid
static constexpr double subsample_factor = 0.1; // Baseline step size coefficient for tau in bin sizes
// Constants for weighted distance sampling
static constexpr uint64_t ks = 10;  // Number of distances to sample per interval
static constexpr uint64_t ncand = 30;        // Candidate pool size for weighted resampling

template<typename T>
class QIE;

struct qstride_t
{
  size_t bix;
  uint64_t nbins;
  uint64_t nmers;
  vec<uint64_t> hdisthist_v; // subsampled prefix-sum rows, flattened [(nbins/G + 1) × W]
};

struct snull_t
{
  size_t qidx;    // query index
  uint64_t start; // bin start
  uint64_t end;   // bin end
  double d;       // estimated distance
};

struct segment_t
{
  uint64_t start;
  uint64_t end;
  double d;
  uint8_t mask; // bitmask: bit i set iff threshold i hits this segment
  char sign;    // '<' for positive thresholds, '>' for negative thresholds
};

struct record_t
{
  uint64_t bix;
  uint64_t L, a, b;
  uint64_t nbins_s;
  char strand;
  bool rstrand; // reference strand: true if it the has lower genome-wide distance, false otherwise
  double d, d_q;
  double percentile = std::numeric_limits<double>::quiet_NaN();
  double fold = std::numeric_limits<double>::quiet_NaN();
  uint8_t mask;
  char sign;

  record_t(uint64_t bix, uint64_t L, uint64_t a, uint64_t b, char strand, double d_q, const segment_t& ab, bool rstrand)
    : bix(bix)
    , L(L)
    , a(a)
    , b(b)
    , nbins_s(ab.end - ab.start)
    , strand(strand)
    , rstrand(rstrand)
    , d(ab.d)
    , d_q(d_q)
    , mask(ab.mask)
    , sign(ab.sign)
  {
  }
};

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
  void extract_intervals_mx(uint64_t tau, size_t idx = 0);
  void extract_intervals_sx(uint64_t tau, size_t idx = 0);
  void expand_intervals(double chisq_th, size_t idx = 0);
  [[nodiscard]] interval_t get_interval(uint64_t i, size_t idx = 0) const;
  const vec<uint64_t>& get_hdisthist() const { return hdisthist_v; }
  const vec<segment_t>& get_segments() const { return segments_v; }
  uint64_t get_nbins() const { return nbins; }
  uint64_t get_nmers() const { return nmers; }
  void map_contiguous_segments(uint8_t th_bv, char sign);
  double estimate_interval_distance(uint64_t a, uint64_t b);
  void extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const;
  bool get_rstrand() const { return rstrand; }
  void set_rstrand(bool v) { rstrand = v; }
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
  const llh_sptr_t<T> llhf;  // log-likelihood function for all calculations
  const uint64_t nbins;      // number of bins
  const uint64_t nmers;      // number of k-mers in query (for per-k-mer HD tracking)
  uint64_t t_q = 0;          // total number of k-mers hits below hdist_th per query sequence
  uint64_t u_q = 0;          // total number misses per query sequence
  bool rstrand = false;      // reference strand: set after per-query MLE comparison
  vec<uint64_t> hdisthist_v; // D[i][j] is the number of hits with HD=j, [(nbins+1) × (hdist_th+1)] row-major; D[0][j]=0
  vec<T> fdc_v;              // The f' contribution c_i of the k-mer (bin) starting at i
  vec<T> sdc_v;              // The f'' contribution s_i of the k-mer (bin) starting at i
  vec<T> fdps_v;             // C[i] = sum(c_0, ..., c_{i}), C[0] = 0 (length n) (shifted by 1 w.r.t. fdc_v)
  vec<T> sdps_v;             // S[i] = sum(s_0, ..., s_{i}), S[0] = 0 (length n) (shifted by 1 w.r.t. sdc_v)
  vec<T> fdpmax_v;           // H[i] = max(C_1, ..., C_{i}), H_0 = -inf, H_{n+1} = inf (length n+1)
  vec<T> fdsmin_v;           // L[i] = min(C_{i}, ..., C_n), L_0 = inf, L_{n+1}= -inf (length n+1)
  arr<vec<interval_t>, WIDTH> rintervals_v; // inital-raw intervals without any postprocessing; released after postprocessing
  arr<vec<interval_t>, WIDTH>
    eintervals_v; // final intervals after expansing and merging; used to compute segments or directly reported in enum-only mode
  vec<segment_t>
    segments_v; // segments computed from eintervals_v based on boundaries across different distance thresholds, constitutes the main output
};

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

private:
  double compute_mle_dist(const vec<uint64_t>& v, uint64_t u, uint64_t t);
  void search_mers(const char* cseq, uint64_t len, DIM<T>& dim_fw, DIM<T>& dim_rc);
  void collect_segments(const DIM<T>& dim, double d_q, bool rc);
  void save_qstride(const DIM<T>& dim);
  void fit_gamma_significance();
  void sample_distances(uint64_t L, vec<snull_t>& samples_v) const;
  void emit_segments(std::ostream& sout, const str& rid) const;
  void report_intervals(std::ostream& sout, const str& rid, DIM<T>& dim, bool rc, size_t idx = 0);

  const params_t<T>& params;
  const sketch_sptr_t sketch;
  const lshf_sptr_t lshf;
  const vec<str>& seq_batch;
  const vec<str>& qid_batch;
  const uint64_t batch_size;
  const uint32_t k;
  const uint32_t h;
  const uint32_t m;
  llh_sptr_t<T> llhf;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint64_t onmers; // Number of observed k-mers in current query (e.g., due to Ns)
  uint64_t enmers; // Number of expected k-mers in current query (= len - k + 1)
  uint64_t nbins;  // Number of bins  (= ceil(enmers / bin_len))
  uint64_t bix;    // Index of the current query in the this batch

  double d_acc = std::numeric_limits<double>::quiet_NaN();
  uint64_t u_acc = 0;
  vec<uint64_t> v_acc;

  uint64_t stride_len = 1;

  vec<qstride_t> qstrides_v;
  vec<record_t> records_v;
};

template<typename... Args>
inline std::ostream& write_tsv(std::ostream& os, const Args&... args)
{
  size_t n = 0;
  ((os << (n++ ? "\t" : "") << args), ...);
  return os;
}

#endif
