#ifndef _GIDIFF_H
#define _GIDIFF_H

#include "common.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "map.hpp"
#include "sketch.hpp"
#include "table.hpp"
#include <CLI.hpp>

const auto url_validator = CLI::Validator(
  [](std::string& input) {
    const std::regex url_regexp = std::regex(
      R"(^(?:(?:https?|ftp)://)(?:\S+@)?(?:(?!10(?:\.\d{1,3}){3})(?!127(?:\.\d{1,3}){3})(?!169\.254(?:\.\d{1,3}){2})(?!192\.168(?:\.\d{1,3}){2})(?!172\.(?:1[6-9]|2\d|3[0-1])(?:\.\d{1,3}){2})(?:[1-9]\d?|1\d\d|2[01]\d|22[0-3])(?:\.(?:1?\d{1,2}|2[0-4]\d|25[0-5])){2}(?:\.(?:[1-9]\d?|1\d\d|2[0-4]\d|25[0-4]))|(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+(?:\.(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+)*(?:\.(?:[a-z\u00a1-\uffff]{2,})))(?::\d{2,5})?(?:/\S*)?$)");
    if (std::regex_match(input, url_regexp)) {
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
  void save_configuration(std::ofstream& cfg_stream);
  void set_sketch_defaults()
  {
    k = 26;
    w = k + 6;
    h = 10;
    m = 4;
    r = 1;
    frac = true;
    nrows = pow(2, 2 * h - 1);
    sdust_t = 0;
    sdust_w = 0;
  }
  bool validate_configuration()
  {
    bool is_invalid = true;
    if (is_invalid = w < k) {
      std::cerr << "The minimum minimizer window size (-w) is k (-k)." << std::endl;
    }
    if (is_invalid = h < 3) {
      std::cerr << "The minimum number of LSH positions (-h) is 3." << std::endl;
    }
    if (is_invalid = h > 15) {
      std::cerr << "The maximum number of LSH positions (-h) is 15." << std::endl;
    }
    if (is_invalid = k > 31) {
      std::cerr << "The maximum allowed k-mer length (-k) is 31." << std::endl;
    }
    if (is_invalid = k < 19) {
      std::cerr << "The minimum allowed k-mer length (-k) is 19." << std::endl;
    }
    if (is_invalid = (k - h) > 16) {
      std::cerr << "For compact k-mer encodings, h must be >= k-16." << std::endl;
    }
    if (sdust_t == 0 || sdust_w == 0) {
      std::cerr << "Setting --sdust-w or --sdust-t to 0 will disable dustmasker." << std::endl;
    }
    return !is_invalid;
  }

protected:
  uint8_t w;
  uint8_t k;
  uint8_t h;
  bool frac;
  uint32_t m;
  uint32_t r;
  uint32_t nrows;
  uint32_t sdust_t = 0;
  uint32_t sdust_w = 0;
  lshf_sptr_t lshf = nullptr;
};

class SketchTarget : public BaseLSH
{
public:
  SketchTarget(CLI::App& sc);
  void create_sketch();
  void save_sketch();

private:
  double rho;
  std::string input;
  std::filesystem::path sketch_path;
  sflatht_sptr_t sketch_sflatht = nullptr;
};

class MapSketch
{
public:
  MapSketch(CLI::App& sc);
  void load_sketch();
  void map_sequences();
  void header_dreport(strstream& dreport_stream);
  uint32_t get_total_qseq() { return total_qseq; }

private:
  sketch_sptr_t sketch = nullptr;
  std::filesystem::path sketch_path;
  std::string query;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  uint32_t hdist_th = 4;
  uint64_t total_qseq = 0;
  uint64_t min_length = 0;
  double dist_th = 0;
  double chi_sq = 3.841; // 95%
  bool divergent = false;
};

#endif
