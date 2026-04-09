#ifndef _GDIFF_HPP
#define _GDIFF_HPP

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
#include "sketch.hpp"
#include "hm.hpp"
#include "CLI11.hpp"

#define VERSION "v0.0.0"
#define PRINT_VERSION std::cerr << "??? version: " << VERSION << std::endl;
#define STRSTREAM_PRECISION 4

extern uint32_t num_threads;
extern str invocation;

const auto url_validator = CLI::Validator(
  [](std::string& input) {
    if (match_url(input)) {
      return std::string("");
    } else {
      return "Given URL is not valid: " + input;
    }
  },
  "URL",
  "URL validator");

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
    m = 2;
    r = 1;
    frac = true;
    nrows = pow(2, (2 * h) - 1);
  }

protected:
  uint8_t w;
  uint8_t k;
  uint8_t h;
  bool frac;
  uint32_t m;
  uint32_t r;
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
  double rho;
  str input_path;
  std::filesystem::path sketch_path;
  sfhm_sptr_t sketch_sfhm = nullptr;
};

class MapSC
{
public:
  MapSC(CLI::App& sc);
  void map();
  void write_header();
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
  uint64_t total_qseq = 0;
  double chisq = 33.00051; // 1e-10
  bool enum_only = false;
  std::vector<double> dist_th;
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
