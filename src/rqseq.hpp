#ifndef _RQSEQ_HPP
#define _RQSEQ_HPP

#include <zlib.h>
#include <filesystem>
#if defined(_LCURL) && _LCURL == 1
  #include <curl/curl.h>
#endif
#include "common.hpp"
#include "msg.hpp"
#include "types.hpp"
#include "lshf.hpp"
#include "enc.hpp"
#include "buckets.hpp"
#include "exthash.hpp"
#include "hyperloglog.hpp"

class HandlerURL
{
protected:
#if defined(_LCURL) && _LCURL == 1
  static size_t write_data(void* ptr, size_t s, size_t nmb, FILE* fst)
  {
    size_t nitems = fwrite(ptr, s, nmb, fst);
    return nitems;
  }

  str download_url(str url)
  {
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();
    if (!std::filesystem::exists(tmp_dir) || !std::filesystem::is_directory(tmp_dir)) {
      error_exit(str("Failed to get temp directory: ") + tmp_dir.string());
    }
    str hash_str = std::to_string(ghhp(url));
    str tmp_filename = "rseq_" + hash_str + ".tmp";
    std::filesystem::path tmp_path = tmp_dir / tmp_filename;

    FILE* fp = fopen(tmp_path.string().c_str(), "wb");
    if (!fp) {
      error_exit(str("Failed to open temp file for writing: ") + tmp_path.string());
    }
    CURL* curl = curl_easy_init();
    if (!curl) {
      error_exit("Failed to initialize CURL.");
    }
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    CURLcode resb = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    if (resb != CURLE_OK) {
      error_exit(str("CURL download failed: ") + curl_easy_strerror(resb));
    }
    fclose(fp);

    return tmp_path.string();
  }
#endif
};

extern "C"
{
#include "kseq.h"
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-function"
KSEQ_INIT(gzFile, gzread)
#pragma GCC diagnostic pop

struct hmer_t
{
  uint64_t x = 0, y = 0, z = 0;
};

class RSeq : public HandlerURL
{
public:
  RSeq(const str& input, const LSHF& lshf, uint8_t w, uint32_t fracth, bool canonical);
  ~RSeq();
  bool set_curr_seq();
  bool read_next_seq();
  // HLL estimate of the number of distinct canonical k-mers seen so far.
  [[nodiscard]] double get_card() const { return csk.estimate(); }
  [[nodiscard]] const char* get_name() const { return name; }
  [[nodiscard]] uint64_t get_len() const { return len; }
  // Append the selected minimizers of the current sequence as packed keys.
  void extract_mers(vec<uint64_t>& keys_v);

private:
  gzFile gfile;
  kseq_t* kseq;
  bool is_url;
  uint8_t k;
  uint8_t w;
  uint32_t fracth;
  bool canonical;
  char* cseq = nullptr;
  char* name = nullptr;
  uint64_t len = 0;
  const LSHF& lshf;
  uint64_t mask_bp = 0;
  uint64_t mask_lr = 0;
  hll::HyperLogLog csk;
  std::filesystem::path input_path;
};

// One query sequence with its ID (coupled: always same index, same lifetime).
struct qseq_t
{
  str qid;
  str seq;
};

constexpr uint64_t bpmax_batch_default = uint64_t(64) << 20;

class QSeq : public HandlerURL
{
public:
  explicit QSeq(const str& input, uint64_t bpmax_batch = bpmax_batch_default);
  ~QSeq();
  // Append up to one batch of sequences; false once the input is exhausted.
  bool read_next_batch();
  void clear();
  bool is_empty();
  const vec<qseq_t>& get_batch_v() const { return batch_v; }

private:
  gzFile gfile;
  kseq_t* kseq;
  bool is_url;
  vec<qseq_t> batch_v;
  uint64_t rbatch_size = 512;
  uint64_t bpmax_batch;
  std::filesystem::path input_path;
};

#endif
