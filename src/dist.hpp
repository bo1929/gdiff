#ifndef _DIST_HPP
#define _DIST_HPP

#include "CLI11.hpp"
#include "llh.hpp"
#include "stils.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include "tpool.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

std::pair<double, char> select_strand_distance(double d_fw, double d_rc);

// Linear-interpolation quantile of a sorted, non-empty vector.
double linear_quantile(const vec<double>& v, double p);

xy_t bracket_distance(double d, const vec<double>& th_v);

// Shared binning validation for map/dist/detect w.r.t. its length.
bool validate_binning(uint64_t bin_shift, uint64_t tau);

class DistanceSampler
{
public:
  DistanceSampler(const sketch_sptr_t& sketch,
                  const vec<qseq_t>& batch_v,
                  uint64_t tau,
                  uint64_t bin_shift,
                  uint32_t hdist_th);
  void run_for_all(uint64_t sample_size, bool keep_counts, ThreadPool& pool);
  void run_per_sequence(uint64_t sample_size, bool keep_counts, ThreadPool& pool);
  void collect_distances(vec<double>& d_v) const;
  void collect_distances(vec<vec<double>>& d_vvec) const;

  template<typename Fn>
  void for_each_sample(Fn&& fn) const
  {
    for (const auto& sch : schemes_v) {
      for (uint64_t s = 0; s < sch.nsamples; ++s)
        fn(sch.bix, sch.enmers, sch.starts_v[s], sch.d_v[s], sch.strand_v[s]);
    }
  }

  // When keep_counts was enabled and the sample has a valid distance, hist/u are
  // the selected-strand counts; otherwise hist is nullptr.
  template<typename Fn>
  void for_each_sample_counts(Fn&& fn) const
  {
    for (const auto& sch : schemes_v) {
      for (uint64_t s = 0; s < sch.nsamples; ++s) {
        const uint64_t* hist = nullptr;
        uint64_t u = 0;
        if (sch.keep_counts && is_valid_distance(sch.d_v[s])) {
          hist = sch.hist_v.data() + s * (hdist_bound + 1);
          u = sch.u_v[s];
        }
        fn(sch.bix, sch.enmers, sch.starts_v[s], sch.d_v[s], sch.strand_v[s], hist, u);
      }
    }
  }

  uint64_t get_nsamples() const;
  uint64_t get_nwinmers() const { return nwinmers; }
  const LLH<double>& get_llhf() const { return llhf; }

private:
  struct scheme_t
  {
    uint64_t bix = 0;
    uint64_t nsamples = 0;
    uint64_t enmers = 0;
    uint64_t nbins = 0;
    vec<uint64_t> starts_v;
    vec<double> d_v;
    vec<char> strand_v;
    bool keep_counts = false;
    vec<uint64_t> hist_v;
    vec<uint64_t> u_v;

    scheme_t(uint64_t bix, uint64_t nsamples, uint64_t enmers, uint64_t nbins, vec<uint64_t> starts_v, bool keep_counts)
      : bix(bix)
      , nsamples(nsamples)
      , enmers(enmers)
      , nbins(nbins)
      , starts_v(std::move(starts_v))
      , d_v(nsamples, nanx())
      , strand_v(nsamples, '+')
      , keep_counts(keep_counts)
      , hist_v(keep_counts ? nsamples * (hdist_bound + 1) : 0)
      , u_v(keep_counts ? nsamples : 0)
    {
    }
  };

  void build_for_all(uint64_t sample_size, bool keep_counts);
  void build_per_sequence(uint64_t sample_size, bool keep_counts);
  void evaluate(ThreadPool& pool);

  sketch_sptr_t sketch;
  const vec<qseq_t>& batch_v;
  bool canonical;
  uint32_t hdist_th;
  uint64_t tau;
  uint64_t tau_bin;
  uint64_t bin_shift;
  uint64_t bin_size;
  uint64_t nwinmers;
  LLH<double> llhf;
  uint32_t k;
  vec<scheme_t> schemes_v;
};

class DistSC
{
public:
  explicit DistSC(CLI::App& sc);
  bool validate_configuration();
  void dist();
  void sample_distances(const sketch_sptr_t& sketch, const vec<qseq_t>& batch_v, strstream& sout, ThreadPool& pool);

private:
  str target_path;
  std::filesystem::path sketch_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  bool output_samples = false;
  uint64_t tau = 0;
  uint64_t sample_size = 200;
  uint64_t bin_shift = 0;
  uint32_t hdist_th = 4;
};

#endif
