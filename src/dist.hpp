#ifndef _DIST_HPP
#define _DIST_HPP

#include "CLI11.hpp"
#include "container.hpp"
#include "distance.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

struct mle_t
{
  double d = nanx();
  double s = nanx();
};

struct samples_t
{
  vec<mle_t> windows_v;
  double d_median = nanx();
  uint64_t n_valid = 0;
  str samples;
};

samples_t process_samples(const Sketch& query, const Sketch& reference, uint32_t hdist_th, bool output_samples, bool is_ba);

// Reconciled summary of one pair: both sides merged rank-wise.
struct summary_t
{
  double d = nanx();         // reported estimate
  double d_median = nanx();  // median over the distances behind `d`
  double d_mean = nanx();    // plain mean before filtering
  double d_highest = nanx(); // max d over every mapped window
  double d_upper = nanx();   // max d over the kept windows only
  vec<double> d_v;           // distances after symmetrization
  uint64_t n_na = 0;         // windows without a distance estimate
  uint64_t n_ub = 0;         // windows whose likelihood-ratio statistic is zero
  uint64_t n_filtered = 0;   // windows rejected by the likelihood-ratio filter
};

summary_t summarize_symmetric(vec<mle_t> ab_v, vec<mle_t> ba_v, double lr_th, double min_portion);

// Compare sets of sketches: every pair, both sides, reconcile data points after sorting.
class DistSC
{
public:
  struct params
  {
    bool output_samples = false;
    uint32_t hdist_th = 3;
    double lr_th = 6.635;
    double min_portion = 0.66;
  };

  explicit DistSC(CLI::App& sc);
  bool validate_configuration();
  void dist();

private:
  // One sketch of one container, resolved from a positional or a list file.
  struct entry_t
  {
    const Container* file = nullptr;
    uint32_t rix = 0;
    str rname;
  };

  // An unordered pair scheduled as two one-way jobs.
  struct pair_t
  {
    uint32_t a = 0;
    uint32_t b = 0;
    samples_t ab;
    samples_t ba;
  };

  // Every unordered pair within set A, or the cross product of A and B.
  void estimate_distances();
  const Container* container_for(const std::filesystem::path& path);
  void resolve_container(const std::filesystem::path& path, vec<uint32_t>& out);
  void resolve_list(const std::filesystem::path& path, vec<uint32_t>& out);
  void write_pair_line(std::ostream& os, const pair_t& pr, const str& name_a, const str& name_b);
  void emit_header(std::ostream& os) const;

  vec<std::unique_ptr<Container>> containers_v; // stable addresses
  vec<entry_t> entries_v;
  vec<uint32_t> set_a_v;
  vec<uint32_t> set_b_v;
  std::filesystem::path container_a_path;
  std::filesystem::path container_b_path;
  std::filesystem::path list_a_path;
  std::filesystem::path list_b_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  params params;
};

#endif
