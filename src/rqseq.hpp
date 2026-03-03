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
#include "hm.hpp"
#include "exthash.hpp"
#include "hyperloglog.hpp"

/* #define CANONICAL */
#define RBATCH_SIZE 512

class HandlerURL
{
protected:
#if defined(_LCURL) && _LCURL == 1
  static size_t write_data(void* ptr, size_t s, size_t nmb, FILE* fst)
  {
    size_t written = fwrite(ptr, s, nmb, fst);
    return written;
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

KSEQ_INIT(gzFile, gzread)

class RSeq : public HandlerURL
{
public:
  RSeq(str input, lshf_sptr_t lshf, uint8_t w, uint32_t r, bool frac);
  ~RSeq();
  bool set_curr_seq();
  bool read_next_seq();
  void compute_rho();
  double get_rho();
  template<typename T>
  void extract_mers(vvec<T>& table);

private:
  gzFile gfile;
  kseq_t* kseq;
  bool is_url;
  uint8_t k;
  uint8_t w;
  uint32_t m;
  uint32_t r;
  bool frac;
  char* cseq;
  char* name;
  uint64_t len;
  lshf_sptr_t lshf;
  uint64_t mask_bp = 0;
  uint64_t mask_lr = 0;
  double n1_est = 0;
  double n2_est = 0;
  double rho = 1.0;
  std::filesystem::path input_path;
};

class QSeq : public HandlerURL
{
  friend class QIE<double>;
  friend class QIE<cm512_t>;

public:
  QSeq(str input);
  ~QSeq();
  bool read_next_batch();
  void clear();
  bool is_empty();
  uint64_t get_cbatch();
  const vec<str>& get_seq_batch() const { return seq_batch; }
  const vec<str>& get_qid_batch() const { return qid_batch; }

private:
  gzFile gfile;
  kseq_t* kseq;
  bool is_url;
  vec<str> seq_batch;
  vec<str> qid_batch;
  uint64_t cbatch_size = 0;
  uint64_t rbatch_size = RBATCH_SIZE;
  std::filesystem::path input_path;
};

#endif
