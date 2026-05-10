#ifndef _MAP_HPP
#define _MAP_HPP

#include <utility>
#include <sstream>
#include <simde/x86/avx512.h>
#include "gamma.hpp"
#include "llh.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"

static constexpr double eps = 1e-7;        // Tolerance for floating-point comparisons
static constexpr uint32_t hdist_bound = 7; // Hamming distance bound for the SIMD alignment

template<typename T>
class QIE;

struct mesh_t
{
  size_t bix;
  uint64_t nbins;
  uint64_t nmers;
  vec<uint64_t> points_v; // sorted bin positions; points_v[0]=0, points_v.back()=nbins
  vec<uint64_t> hists_v;  // flat, size = points_v.size() * (hdist_th + 1);

  mesh_t(size_t bix, uint64_t nbins, uint64_t nmers)
    : bix(bix)
    , nbins(nbins)
    , nmers(nmers)
  {
  }
};

struct record_t
{
  uint64_t bix;
  uint64_t L;
  uint64_t a;     // 1-based reported start coordinate on the source query
  uint64_t b;     // 1-based reported end/boundary coordinate on the source query
  uint64_t nbins; // The number of bins in the interval [a_bin, b_bin)
  uint64_t a_bin; // 1-based inclusive bin coordinate on the source query
  uint64_t b_bin; // 1-based exclusive bin coordinate on the source query
  bool is_rc;     // Is the source query on the reverse complement strand?
  bool is_ref;    // Is this record from the strand with the lower MLE distance?
  double d;       // MLE distance for this interval
  double d_q;     // MLE distance of the source query (e.g., contig)
  double I;       // observed Fisher information I(d)
  uint8_t mask;   // Which thresholds cover the interval [a_bin, b_bin)?
  double d_low;   // Lower bound of the satisfied distance thresholds
  double d_high;  // Upper bound of the satisfied distance threshold
  double percentile = std::numeric_limits<double>::quiet_NaN();
  double fold = std::numeric_limits<double>::quiet_NaN();

  bool is_intact() const { return a == 1 && b == L; }
  // 0-based half-open bin-boundary interval for histogram extraction and null-window overlap.
  interval_t get_interval() const { return {a_bin - 1, b_bin - 1}; }

  record_t(uint64_t bix,
           uint64_t L,
           uint64_t a,
           uint64_t b,
           uint64_t a_bin,
           uint64_t b_bin,
           bool is_rc,
           bool is_ref,
           double d,
           double d_q,
           double I,
           uint8_t mask,
           xy_t d_range)
    : bix(bix)
    , L(L)
    , a(a)
    , b(b)
    , nbins(b_bin - a_bin)
    , a_bin(a_bin)
    , b_bin(b_bin)
    , is_rc(is_rc)
    , is_ref(is_ref)
    , d(d)
    , d_q(d_q)
    , I(I)
    , mask(mask)
    , d_low(d_range.first)
    , d_high(d_range.second)
  {
  }
};

inline interval_t get_rinterval(const interval_t& ab_bin, uint64_t bin_shift, uint64_t enmers, uint32_t k, bool is_last)
{
  const uint64_t a = ((ab_bin.first - 1) << bin_shift) + 1;
  const uint64_t b = is_last ? (enmers + k - 1) : std::min(((ab_bin.second - 1) << bin_shift) + 1, enmers + 1);
  return {a, b};
}

inline interval_t get_einterval(const interval_t& ab, uint64_t bin_shift, uint64_t enmers, uint32_t k)
{
  const uint64_t a = ((ab.first - 1) << bin_shift) + 1;
  const uint64_t b = std::min(ab.second << bin_shift, enmers) + k - 1;
  return {a, b};
}

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
  void extract_intervals_mx(uint64_t tau, size_t ix = 0);
  void extract_intervals_sx(uint64_t tau, size_t ix = 0);
  void expand_intervals(double chisq_th, size_t ix = 0);
  void extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const;
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
  uint64_t get_nbins() const { return nbins; }
  uint64_t get_nmers() const { return nmers; }
  const vec<uint64_t>& get_hdisthist() const { return hdisthist_v; }
  const vec<interval_t>& get_eintervals(size_t ti) const { return eintervals_v[ti]; }
  [[nodiscard]] interval_t
  get_interval(uint64_t i, size_t ix = 0) const; // i-th merged interval in 1-based inclusive bin coordinates

private:
  const params_t<T>& params;
  const llh_sptr_t<T> llhf;  // log-likelihood function for all calculations
  const uint64_t nbins;      // number of bins
  const uint64_t nmers;      // number of k-mers in query (for per-k-mer HD tracking)
  uint64_t t_q = 0;          // total number of k-mers hits below hdist_th per query sequence
  uint64_t u_q = 0;          // total number misses per query sequence
  vec<uint64_t> hdisthist_v; // D[i][j] is the number of hits with HD=j, [(nbins+1) x (hdist_th+1)] row-major; D[0][j]=0
  vec<T> fdc_v;              // The f' contribution c_i of the k-mer (bin) starting at i
  vec<T> sdc_v;              // The f'' contribution s_i of the k-mer (bin) starting at i
  vec<T> fdps_v;             // C[i] = sum(c_0, ..., c_{i}), C[0] = 0 (length n) (shifted by 1 w.r.t. fdc_v)
  vec<T> sdps_v;             // S[i] = sum(s_0, ..., s_{i}), S[0] = 0 (length n) (shifted by 1 w.r.t. sdc_v)
  vec<T> fdpmax_v;           // H[i] = max(C_1, ..., C_{i}), H_0 = -inf, H_{n+1} = inf (length n+1)
  vec<T> fdsmin_v;           // L[i] = min(C_{i}, ..., C_n), L_0 = inf, L_{n+1}= -inf (length n+1)
  arr<vec<interval_t>, WIDTH> rintervals_v; // 1-based inclusive coordinates of initial intervals prior to postprocessing
  arr<vec<interval_t>, WIDTH> eintervals_v; // 1-based inclusive coordinates of final intervals after expanding and merging
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
  void enumerate_intervals(std::ostream& sout, const str& rid, DIM<T>& dim, bool rc, size_t ix = 0);
  void save_mesh(const DIM<T>& dim);
  void collect_contiguous(DIM<T>& dim, uint8_t th_bv, double d_q, bool is_rc, bool is_ref);
  void test_significance(const uint64_t nsamples = 100, const uint64_t ntries = 250);
  bool sample_mesh_distance(const uint64_t mix, const uint64_t bix, const interval_t& ab);
  bool sample_metropolis_hastings(const xy_t& init, uint64_t S = 50, uint64_t B = 100);
  xy_t score_gamma(const record_t& r, const vec<xy_t>& samples_v) const;
  void report_contiguous(std::ostream& sout, const str& rid) const;

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
  vec<mesh_t> meshes_v;
  vec<record_t> records_v;
  vec<xy_t> samples_v;
  vec<uint64_t> v_acc;
  uint64_t u_acc = 0;
  double d_acc = std::numeric_limits<double>::quiet_NaN();
};

template<typename... Args>
inline std::ostream& write_tsv(std::ostream& os, const Args&... args)
{
  size_t n = 0;
  ((os << (n++ ? "\t" : "") << args), ...);
  return os;
}

#endif
