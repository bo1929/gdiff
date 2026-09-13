#ifndef _DIST_HPP
#define _DIST_HPP

#include "CLI11.hpp"
#include "llh.hpp"
#include "sym.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include "stils.hpp"
#include "tpool.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

std::pair<double, char> select_strand_distance(double d_fw, double d_rc);

xy_t bracket_distance(double d, const vec<double>& th_v);

// Shared binning validation for map/dist/detect w.r.t. its length.
bool validate_binning(uint64_t bin_shift, uint64_t tau);

// One direction: per-window rows plus its summary. `samples` needs --output-samples.
struct dir_result_t
{
  vec<dpoint_t> rows;
  double d_median = nanx();
  uint64_t n_valid = 0;
  uint64_t n_unmapped = 0;
  str samples;
};

// Score a query's stored windows against a reference's buckets.
dir_result_t run_direction(const Sketch& query,
                           const Sketch& reference,
                           const str& dir_label,
                           uint32_t hdist_th,
                           uint64_t nmers_limit,
                           bool output_samples);

// The FASTA-side query path; sketch-side queries use run_direction instead.
class DistanceSampler
{
public:
  DistanceSampler(const Sketch& sketch,
                  const vec<qseq_t>& batch_v,
                  uint64_t tau,
                  uint64_t bin_shift,
                  uint32_t hdist_th);
  void run_for_all(uint64_t sample_size, bool keep_counts, ThreadPool& pool);
  void run_per_sequence(uint64_t sample_size, bool keep_counts, ThreadPool& pool);
  void collect_distances(vec<double>& d_v) const;
  void collect_distances(vec<vec<double>>& d_vvec) const;
  // Per-window (d, lr_ub) rows for reconciliation; requires keep_counts.
  void collect_samples(vec<dpoint_t>& rows) const;
  void collect_samples(vec<vec<dpoint_t>>& rows_per_seq) const;

  template<typename Fn>
  void for_each_sample(Fn&& fn) const
  {
    for (const auto& sch : schemes_v) {
      for (uint64_t s = 0; s < sch.nsamples; ++s)
        fn(sch.bix, sch.enmers, sch.starts_v[s], sch.d_v[s], sch.strand_v[s]);
    }
  }

  // hist/u are the selected-strand counts, or hist is null without keep_counts.
  template<typename Fn>
  void for_each_counts(Fn&& fn) const
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

    scheme_t(uint64_t bix,
             uint64_t nsamples,
             uint64_t enmers,
             uint64_t nbins,
             vec<uint64_t> starts_v,
             bool keep_counts)
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
  dpoint_t sample_row(const scheme_t& sch, uint64_t s) const;

  const Sketch& sketch;
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

private:
  // One record of one container, resolved from a positional or a list file.
  struct member_t
  {
    const SketchFile* file = nullptr;
    uint32_t rec = 0;
    str rname;
  };

  // An unordered pair scheduled as two directional jobs.
  struct pair_t
  {
    uint32_t a = 0;
    uint32_t b = 0;
    dir_result_t ab;
    dir_result_t ba;
  };

  void resolve_source(const std::filesystem::path& path, vec<uint32_t>& out);
  const SketchFile* file_for(const std::filesystem::path& path);
  // Sketch-vs-sketch: all unordered pairs within one set, or the cross product.
  void dist_sketches();
  // FASTA query: forward samples the query, reverse uses its in-memory sketch.
  void dist_fasta();
  void write_pair_row(std::ostream& os, const pair_t& pr, const str& name_a, const str& name_b);
  void emit_header(std::ostream& os) const;

  vec<std::unique_ptr<SketchFile>> files; // stable addresses
  vec<member_t> members;
  vec<uint32_t> set_a;
  vec<uint32_t> set_b;
  std::filesystem::path set_a_path;
  std::filesystem::path set_b_path;
  std::filesystem::path query_list_path;
  std::filesystem::path reference_list_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  bool output_samples = false;
  bool symmetric = true;
  bool no_header = false;
  uint64_t tau = 0;
  uint64_t sample_size = 0;
  uint64_t bin_shift = 0;
  uint32_t hdist_th = 4;
  double lr_th = lr_th_default;
  double min_portion = min_portion_default;
};

#endif
