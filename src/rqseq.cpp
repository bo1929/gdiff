#include "rqseq.hpp"

#define HLL_BUCKET_FACTOR 16

RSeq::RSeq(const str& input, const lshf_sptr_t& lshf, uint8_t w, uint32_t frac_th, bool canonical)
  : w(w)
  , frac_th(frac_th)
  , canonical(canonical)
  , lshf(lshf)
  , csk(HLL_BUCKET_FACTOR)
{
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  k = lshf->get_k();
  mask_bp = u64m >> ((32 - k) * 2);
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));

  is_url = std::regex_match(input, urlexp);
  if (is_url) {
#if defined(_LCURL) && _LCURL == 1
    input_path = download_url(input);
#else
    warn_msg("Failed to download from URL, compiled without libcurl!");
#endif
  } else {
    input_path = input;
  }

  gfile = gzopen(input_path.c_str(), "rb");
  if (gfile == nullptr) {
    error_exit(str("Failed to open the file at ") + input_path.string());
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

bool RSeq::read_next_seq() { return kseq_read(kseq) >= 0; }

bool RSeq::set_curr_seq()
{
  name = kseq->name.s;
  cseq = kseq->seq.s;
  len = kseq->seq.l;
  return len >= w;
}

void RSeq::extract_mers(vec<uint64_t>& keys)
{
  uint8_t ldiff;
  if (w > k) {
    ldiff = w - k + 1;
  } else {
    ldiff = 1;
    w = k;
  }
  uint64_t klix = 0;
  uint64_t orenc64_bp = 0, orenc64_lr = 0;
  vec<hmer_t> winenc_v(ldiff);
  hmer_t cminimizer{};
  uint64_t i, l;
  for (i = l = 0; i < len;) {
    if (SEQ_NT4_TABLE[static_cast<uint8_t>(cseq[i])] >= 4) {
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
    const uint64_t enc_bp = orenc64_bp & mask_bp;
    winenc_v[klix] = {enc_bp, orenc64_lr & mask_lr, xhur64(enc_bp)};
    csk.add(canonical ? xhur64(std::max(enc_bp, revcomp_bp64(enc_bp, k))) : winenc_v[klix].z);
    if (++klix == ldiff) klix = 0;
    if (l < w) {
      continue;
    }
    cminimizer = *std::min_element(
      winenc_v.begin(), winenc_v.end(), [](hmer_t lhs, hmer_t rhs) { return lhs.z < rhs.z; });
    if (canonical) {
      uint64_t rcenc64_bp = revcomp_bp64(cminimizer.x, k);
      if (cminimizer.x < rcenc64_bp) {
        cminimizer.x = rcenc64_bp;
        cminimizer.y = bp64_to_lr64(rcenc64_bp);
      }
    }
    const uint32_t rix = lshf->compute_hash_bp(cminimizer.x);
    if (rix < frac_th) {
      keys.push_back(pack_key(rix, lshf->drop_ppos_lr(cminimizer.y)));
    }
  }
}

QSeq::QSeq(const str& input, uint64_t max_batch_bases)
  : max_batch_bases(max_batch_bases)
{
  is_url = std::regex_match(input, urlexp);
  if (is_url) {
#if defined(_LCURL) && _LCURL == 1
    input_path = download_url(input);
#else
    warn_msg("Failed to download from URL, compiled without libcurl!");
#endif
  } else {
    input_path = input;
  }
  gfile = gzopen(input_path.c_str(), "rb");
  if (gfile == nullptr) {
    error_exit(str("Failed to open the file at ") + input_path.string());
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
  bool cont_reading = false;
  uint64_t ix = 0;
  uint64_t nbases = 0;
  while ((ix < rbatch_size) && (nbases < max_batch_bases) && (cont_reading = (kseq_read(kseq) >= 0))) {
    nbases += kseq->seq.l;
    batch_v.push_back({kseq->name.s, kseq->seq.s});
    ix++;
  }
  cbatch_size = ix;
  return cont_reading;
}

bool QSeq::is_empty() { return batch_v.empty(); }

void QSeq::clear() { batch_v.clear(); }
