#include "rqseq.hpp"
// extern "C"
// {
// #include "sdust.h"
// }

RSeq::RSeq(std::string input, lshf_sptr_t lshf, uint8_t w, uint32_t r, bool frac)
  : w(w)
  , r(r)
  , frac(frac)
  , lshf(lshf)
{
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  k = lshf->get_k();
  m = lshf->get_m();
  mask_bp = u64m >> ((32 - k) * 2);
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));

  is_url = std::regex_match(input, url_regexp);
  if (is_url) {
#if defined _WLCURL && _WLCURL == 1
    input_path = download_url(input);
#else
    warn_msg("Failed to download from URL, compiled without libcurl!");
#endif
  } else {
    input_path = input;
  }

  gfile = gzopen(input_path.c_str(), "rb");
  if (gfile == nullptr) {
    error_exit(std::string("Failed to open the file at ") + input_path.string());
  }
  kseq = kseq_init(gfile);
}

RSeq::~RSeq()
{
  kseq_destroy(kseq);
  gzclose(gfile);
  if (is_url) {
    std::filesystem::remove(input_path);
  }
}

/* void RSeq::compute_rho() { rho = static_cast<double>(wcix) / static_cast<double>(wnix); } */
void RSeq::compute_rho() { rho = n2_est / n1_est; }

bool RSeq::read_next_seq() { return kseq_read(kseq) >= 0; }

double RSeq::get_rho() { return rho; }

bool RSeq::set_curr_seq()
{
  name = kseq->name.s;
  cseq = kseq->seq.s;
  len = kseq->seq.l;
  return len >= w;
}

template<typename T>
void RSeq::extract_mers(vvec<T>& table)
{
  uint32_t i, l;
  uint32_t rix, rix_res;
  uint8_t ldiff;
  if (w > k) {
    ldiff = w - k + 1;
  } else {
    ldiff = 1;
    w = k;
  }
  hll::HyperLogLog c1(12);
  hll::HyperLogLog c2(12);
  uint64_t kix = 0, klix = 0;
  uint64_t orenc64_bp, orenc64_lr, rcenc64_bp;
  std::vector<hmer_t> winenc_v(ldiff);
  hmer_t cminimizer, pminimizer;
  uint32_t mrs = 0, mre = len;
  int mn = 0, mi = 0;
  uint64_t* rgs;
  // if (sdust_t > 0 && sdust_w > 0) rgs = sdust(0, (uint8_t*)cseq, -1, sdust_t, sdust_w, &mn);
  // if (mn > 0) {
  //   mre = (uint32_t)(rgs[mi]);
  //   mrs = (uint32_t)(rgs[mi] >> 32);
  // }
  for (i = l = 0; i < len;) {
    if (SEQ_NT4_TABLE[cseq[i]] >= 4) {
      l = 0, i++;
      continue;
    }
    l++, i++;
    if (l < k) {
      continue;
    }
    if (l == k) {
      compute_encoding(cseq + i - k, cseq + i, orenc64_lr, orenc64_bp);
    } else {
      update_encoding(cseq + i - 1, orenc64_lr, orenc64_bp);
    }
    // Add this block for SDUST masking.
    // if ((mi < mn) && ((i + k) > mrs)) {
    //   c1.add(xhur64(orenc64_bp & mask_bp));
    //   if (i < mre) {
    //     continue;
    //   } else {
    //     mi++;
    //     l = 0;
    //     if ((mi < mn)) {
    //       mre = (uint32_t)(rgs[mi]);
    //       mrs = (uint32_t)(rgs[mi] >> 32);
    //     } else {
    //       free(rgs);
    //     }
    //     continue;
    //   }
    // }
    klix = kix % ldiff;
    winenc_v[klix] = {orenc64_bp & mask_bp, orenc64_lr & mask_lr, xhur64(orenc64_bp & mask_bp)};
    c1.add(winenc_v[klix].z);
    kix++;
    if ((l < w) && (i != len)) {
      continue;
    }
    cminimizer = *std::min_element(winenc_v.begin(), winenc_v.end(), [](hmer_t lhs, hmer_t rhs) { return lhs.z < rhs.z; });
    c2.add(cminimizer.z);
#ifdef CANONICAL
    rcenc64_bp = revcomp_bp64(cminimizer.x, k);
    if (cminimizer.x < rcenc64_bp) {
      cminimizer.x = rcenc64_bp;
      cminimizer.y = bp64_to_lr64(rcenc64_bp);
    }
#endif /* CANONICAL */
    rix = lshf->compute_hash(cminimizer.x);
    rix_res = rix % m;
    if (frac ? rix_res <= r : rix_res == r) {
      rix = frac ? rix / m * (r + 1) + rix_res : rix / m;
      table[rix].push_back(lshf->drop_ppos_lr(cminimizer.y));
      wnix++;
      if (cminimizer.x != pminimizer.x) {
        wcix++;
      }
      pminimizer = cminimizer;
    }
  }
  n1_est += c1.estimate();
  n2_est += c2.estimate();
}

QSeq::QSeq(std::string input)
{
  is_url = std::regex_match(input, url_regexp);
  if (is_url) {
#if defined _WLCURL && _WLCURL == 1
    input_path = download_url(input);
#else
    warn_msg("Failed to download from URL, compiled without libcurl!");
#endif
  } else {
    input_path = input;
  }
  gfile = gzopen(input_path.c_str(), "rb");
  if (gfile == nullptr) {
    error_exit(std::string("Failed to open the file at ") + input_path.string());
  }
  kseq = kseq_init(gfile);
}

QSeq::~QSeq()
{
  kseq_destroy(kseq);
  gzclose(gfile);
}

bool QSeq::read_next_batch()
{
  seq_batch.clear();
  identifier_batch.clear();
  seq_batch.reserve(rbatch_size);
  identifier_batch.reserve(rbatch_size);
  bool cont_reading = false;
  uint64_t ix = 0, bpc = 0;
  // Alternatively, use (ix < rbatch_size).
  while ((bpc < bpc_limit) && (cont_reading = kseq_read(kseq) >= 0)) {
    bpc += kseq->seq.l;
    seq_batch.emplace_back(kseq->seq.s);
    identifier_batch.emplace_back(kseq->name.s);
    ix++;
  }
  cbatch_size = ix;
  return cont_reading;
}

bool QSeq::is_batch_finished()
{
  assert(seq_batch.size() == identifier_batch.size());
  return seq_batch.empty() && identifier_batch.empty();
}

void QSeq::clear_curr_batch()
{
  seq_batch.clear();
  identifier_batch.clear();
}

uint64_t QSeq::get_cbatch_size() { return cbatch_size; }

template void RSeq::extract_mers(vvec<enc_t>& table);
