#ifndef _DIST_HPP
#define _DIST_HPP

#include "CLI11.hpp"
#include "llh.hpp"
#include "maptils.hpp"
#include "rqseq.hpp"
#include "sketch.hpp"
#include "types.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>

class HDHist
{
public:
  HDHist(uint64_t nbins, uint64_t nmers, uint32_t hdist_th, uint64_t bin_shift);

  void aggregate_mer(uint32_t hdist_min, uint64_t i);
  void compute_prefhistsum();
  void extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const;
  [[nodiscard]] uint64_t get_nbins() const { return nbins; }
  [[nodiscard]] uint64_t get_nmers() const { return nmers; }

private:
  uint64_t nbins;
  uint64_t nmers;
  uint32_t hdist_th;
  uint64_t bin_shift;
  vec<uint64_t> hdisthist_v;
};

struct dist_summary_t
{
  uint64_t n = 0;
  double mean = nanx();
  double sd = nanx();
  arr<double, 7> quantiles{};
};

struct dist_sample_t
{
  str qid;
  uint64_t L = 0;
  uint64_t a = 0; // 0-based sequence start
  char strand = '.';
  double d = nanx();
};

dist_summary_t summarize_distances(vec<double> d_v);
std::pair<double, char> select_strand_distance(double d_fw, double d_rc);

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
                        strstream* samples_sout);
  void sample_sequence(const sketch_sptr_t& sketch,
                       const str& seq,
                       const str& qid,
                       uint64_t n_samples,
                       vec<double>& d_v,
                       vec<dist_sample_t>* samples_v);
  void search_mers(const sketch_sptr_t& sketch, const char* cseq, uint64_t len, HDHist& hist) const;
  void search_mers(const sketch_sptr_t& sketch, const char* cseq, uint64_t len, HDHist& hist_fw, HDHist& hist_rc) const;
  void search_window(const sketch_sptr_t& sketch,
                     const char* cseq,
                     uint64_t len,
                     uint64_t j0,
                     uint64_t j1,
                     vec<uint64_t>& v_fw,
                     vec<uint64_t>* v_rc) const;

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
