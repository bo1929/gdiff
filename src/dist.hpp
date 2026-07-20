#ifndef _DIST_HPP
#define _DIST_HPP

#include "CLI11.hpp"
#include "llh.hpp"
#include "maptils.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include "tpool.hpp"
#include "types.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>

class HDHist
{
public:
  // With zero=false the backing store is left uninitialized; the caller must
  // zero it (e.g. zero_range, possibly from multiple threads) before use.
  explicit HDHist(uint64_t nbins, uint64_t nmers, uint32_t hdist_th, uint64_t bin_shift, bool zero = true);

  void aggregate_mer(uint32_t hdist_min, uint64_t i);
  void aggregate_mer_atomic(uint32_t hdist_min, uint64_t i);
  void zero_range(uint64_t r0, uint64_t r1);
  [[nodiscard]] uint64_t storage_rows() const { return nbins + 1; }
  void compute_prefhistsum();
  // Blocked parallel prefix sum over rows; result is identical to
  // compute_prefhistsum() (integer addition is exact and associative here).
  void compute_prefhistsum_parallel(ThreadPool& pool, uint32_t nchunks);
  void extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const;
  [[nodiscard]] uint64_t get_nbins() const { return nbins; }
  [[nodiscard]] uint64_t get_nmers() const { return nmers; }

private:
  uint64_t nbins;
  uint64_t nmers;
  uint32_t hdist_th;
  uint64_t bin_shift;
  std::unique_ptr<uint64_t[]> hdisthist_v;
};

struct dist_summary_t
{
  uint64_t n = 0;
  double mean = nanx();
  double sd = nanx();
  arr<double, 7> quantiles{};
};

dist_summary_t summarize_distances(vec<double> d_v);
std::pair<double, char> select_strand_distance(double d_fw, double d_rc);

// Draws n_samples uniform region starts in [0, npos) (with replacement).
vec<uint64_t> sample_region_starts(uint64_t npos, uint64_t n_samples, std::mt19937& rng);

// Weighted reservoir slot claim: each of the `total` slots is claimed by a
// sequence spanning cumulative-length interval [offset, offset+size) with
// probability size / (offset+size). Returns the indices of claimed slots.
vec<size_t> select_reservoir_slots(uint64_t offset, uint64_t size, uint64_t total, std::mt19937& rng);

class DistSC
{
public:
  explicit DistSC(CLI::App& sc);
  void dist();
  bool validate_configuration();

private:
  void sample_sequences(const sketch_sptr_t& sketch,
                        const vec<str>& seq_batch,
                        const vec<str>& qid_batch,
                        strstream& sout,
                        strstream* samples_sout,
                        ThreadPool& pool);

  str query_path;
  std::filesystem::path sketch_path;
  std::filesystem::path output_path;
  std::filesystem::path samples_output_path;
  std::ofstream output_file;
  std::ofstream samples_output_file;
  std::ostream* output_stream = &std::cout;
  std::ostream* samples_output_stream = nullptr;
  uint64_t tau = 0;
  uint64_t sample_size = 200;
  uint64_t bin_shift = 0;
  uint32_t hdist_th = 4;
};

#endif
