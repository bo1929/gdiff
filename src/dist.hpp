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

class HDHistogram
{
public:
  HDHistogram(uint64_t nmers, uint32_t hdist_th);

  void aggregate_mer(uint32_t hdist, uint64_t pos);
  void inclusive_scan();
  void extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& hist, uint64_t& misses, uint64_t& hits) const;

private:
  uint64_t nmers;
  uint32_t hdist_th;
  uint32_t width;
  vec<uint64_t> hist_v;
};

struct dist_summary_t
{
  uint64_t n = 0;
  double mean = nanx();
  double sd = nanx();
  arr<double, 7> quantiles{};
};

dist_summary_t summarize_distances(vec<double> distances);
vec<uint64_t> sample_region_starts(uint64_t seq_len, uint64_t region_len, uint64_t sample_size, std::mt19937& rng);
std::pair<double, char> select_strand_distance(double d_fw, double d_rc);

class DistSC
{
public:
  explicit DistSC(CLI::App& sc);
  void dist();

private:
  void process_pair(const sketch_sptr_t& sketch, const str& seq, const str& qid);
  void
  search_mers(const sketch_sptr_t& sketch, const char* cseq, uint64_t len, HDHistogram& hist_fw, HDHistogram* hist_rc) const;
  void write_summary(const str& qid, const str& rid, const dist_summary_t& summary);
  void write_sample(const str& qid, uint64_t seq_len, uint64_t start, char strand, const str& rid, double distance);

  str query_path;
  std::filesystem::path sketch_path;
  std::filesystem::path output_path;
  std::filesystem::path samples_output_path;
  std::ofstream output_file;
  std::ofstream samples_output_file;
  std::ostream* output_stream = &std::cout;
  std::ostream* samples_output_stream = nullptr;
  uint64_t region_len = 0;
  uint64_t sample_size = 200;
  uint32_t hdist_th = 4;
};

#endif
