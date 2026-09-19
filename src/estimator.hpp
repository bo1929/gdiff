#ifndef _ESTIMATOR_HPP
#define _ESTIMATOR_HPP

// Estimating a distance by scoring windows against a reference's index: per-window scoring
// of a query sketch's stored windows (run_direction), and the estimator that samples fresh
// windows from query sequences (WindowEstimator). Shared by dist and detect.

#include "distance.hpp"
#include "llh.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include "sym.hpp"
#include "tpool.hpp"

// Scalar LLH for one reference/threshold pair; derivatives are not needed.
inline LLH<double> make_llhf(const Sketch& sketch, uint32_t hdist_th)
{
  return {sketch.get_k(), sketch.get_h(), sketch.get_rho(), hdist_th, 0.0, false};
}

// The nearer of two strand distances, with the strand that produced it.
inline std::pair<double, char> select_strand_distance(double d_fw, double d_rc)
{
  const bool fw_valid = is_valid_distance(d_fw);
  const bool rc_valid = is_valid_distance(d_rc);
  if (!fw_valid && !rc_valid) return {nanx(), '.'};
  if (!rc_valid || (fw_valid && d_fw <= d_rc)) return {d_fw, '+'};
  return {d_rc, '-'};
}

// Scoring one direction of a pair: per-window rows plus its summary. `samples` needs
// --output-samples.
struct direction_t
{
  vec<ds_t> rows_v;
  double d_median = nanx();
  uint64_t n_valid = 0;
  str samples;
};

// Score a query's stored windows against a reference's buckets. Sample rows name the pair in
// canonical (a < b) order, so `is_ba` says the query is the B member; `dist` renders samples.
direction_t run_direction(const Sketch& query, const Sketch& reference, uint32_t hdist_th, bool output_samples, bool is_ba);

// The FASTA-side query path; sketch-side queries use run_direction instead.
class WindowEstimator
{
public:
  WindowEstimator(const Sketch& sketch, const vec<qseq_t>& batch_v, uint64_t tau, uint64_t bin_shift, uint32_t hdist_th);
  void estimate_all(uint64_t sample_size, ThreadPool& pool);
  void estimate_per_sequence(uint64_t sample_size, ThreadPool& pool);

  // Visit every sampled window in scheme order. The callback takes
  // (bix, enmers, start, d, strand, hist, u); hist is null and u zero for an unmapped window,
  // which has no distance and therefore no counts.
  template<typename Fn>
  void for_each_window(Fn&& fn) const
  {
    for (const auto& sch : schemes_v) {
      for (uint64_t s = 0; s < sch.nsamples; ++s) {
        const uint64_t* hist = nullptr;
        uint64_t u = 0;
        if (is_valid_distance(sch.d_v[s])) {
          hist = sch.hist_v.data() + s * (hdist_bound + 1);
          u = sch.u_v[s];
        }
        fn(sch.bix, sch.enmers, sch.starts_v[s], sch.d_v[s], sch.strand_v[s], hist, u);
      }
    }
  }

  // The (d, s) row behind one visited window; s stays NaN for an unmapped window.
  ds_t make_row(double d, const uint64_t* hist, uint64_t u) const;

  uint64_t nwins() const;

private:
  struct scheme_t
  {
    uint64_t bix = 0;
    uint64_t nsamples = 0;
    uint64_t enmers = 0;
    vec<uint64_t> starts_v;
    vec<double> d_v;
    vec<char> strand_v;
    vec<uint64_t> hist_v;
    vec<uint64_t> u_v;

    scheme_t(uint64_t bix, uint64_t nsamples, uint64_t enmers, vec<uint64_t> starts_v)
      : bix(bix)
      , nsamples(nsamples)
      , enmers(enmers)
      , starts_v(std::move(starts_v))
      , d_v(nsamples, nanx())
      , strand_v(nsamples, '+')
      , hist_v(nsamples * (hdist_bound + 1))
      , u_v(nsamples)
    {
    }
  };

  void plan_for_all(uint64_t sample_size);
  void plan_per_sequence(uint64_t sample_size);
  void evaluate(ThreadPool& pool);

  const Sketch& sketch;
  const vec<qseq_t>& batch_v;
  bool canonical;
  uint32_t hdist_th;
  uint64_t tau;
  uint64_t bin_shift;
  uint64_t nwinmers;
  LLH<double> llhf;
  vec<scheme_t> schemes_v;
};

#endif
