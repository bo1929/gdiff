#ifndef _GDIFF_HPP
#define _GDIFF_HPP

#include <cmath>
#include <atomic>
#include <chrono>
#include <ctime>
#include <mutex>
#include <thread>
#include "msg.hpp"
#include "common.hpp"
#include "types.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "map.hpp"
#include "dist.hpp"
#include "sketch.hpp"
#include "hm.hpp"
#include "CLI11.hpp"

#define VERSION "v0.0.0"
#define PRINT_VERSION std::cerr << "??? version: " << VERSION << std::endl;
#define STRSTREAM_PRECISION 4

extern uint32_t num_threads;
extern str invocation;

class BaseLSH
{
public:
  void set_lshf();
  void set_nrows();
  void set_sketch_defaults()
  {
    k = 27;
    w = k + 6;
    h = 11;
    rate = 1.0;
    canonical = true;
    nrows = uint32_t(1) << (2 * h); // recomputed by set_nrows()
  }

protected:
  uint8_t w;
  uint8_t k;
  uint8_t h;
  double rate = 1.0; // k-mer sampling rate on top of minimizers: keep if LSH(x) < rate * 2^(2h)
  bool canonical = true;
  uint32_t nrows;
  lshf_sptr_t lshf = nullptr;
};

class SketchSC : public BaseLSH
{
public:
  SketchSC(CLI::App& sc);
  void create();
  void save();
  bool validate_configuration();
  void write_header(std::ofstream& stream);
  void write_config(std::ofstream& stream);

private:
  str input_path;
  std::filesystem::path sketch_path;
  sfhm_sptr_t sketch_sfhm = nullptr;
  double rho;
};

class MapSC
{
public:
  MapSC(CLI::App& sc);
  void map();
  bool validate_configuration();
  uint64_t get_total_qseq() const { return total_qseq; }

private:
  str query_path;
  std::filesystem::path sketch_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  uint32_t hdist_th = 4;
  uint64_t tau = 1;
  uint64_t bin_shift = 0;
  double chisq = 33.00051; // 1e-10
  uint64_t sample_size = 200;
  bool enum_only = false;
  std::vector<double> dist_th;
  uint64_t total_qseq = 0;
};

class MergeSC
{
public:
  MergeSC(CLI::App& sc);
  void merge();

private:
  std::filesystem::path output_path;
  std::vector<str> sketch_paths;
};

class InfoSC
{
public:
  InfoSC(CLI::App& sc);
  void info();

private:
  std::filesystem::path sketch_path;
};

#endif
