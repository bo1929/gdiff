#include "rqseq.hpp"

static constexpr uint8_t default_hll_factor = 16;

RSeq::RSeq(const str& input, const LSHF& lshf, uint8_t w, uint32_t fracth, bool canonical)
  : w(w)
  , fracth(fracth)
  , canonical(canonical)
  , lshf(lshf)
  , csk(default_hll_factor)
{
  k = lshf.get_k();
  const lsh_masks_t masks = get_lsh_masks(k);
  mask_bp = masks.bp;
  mask_lr = masks.lr;

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

void RSeq::extract_mers(vec<uint64_t>& keys_v)
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
  uint64_t cx_prev;
  bool is_wmin = false;
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
    const uint64_t enc_lr = orenc64_lr & mask_lr;
    uint64_t x = enc_bp, y = enc_lr;
    if (canonical) {
      const uint64_t rc_bp = revcomp_bp64(enc_bp, k);
      if (x < rc_bp) {
        x = rc_bp;
        y = bp64_to_lr64(rc_bp);
      }
    }
    const uint64_t z = xhur64(x);
    winenc_v[klix] = {x, y, z};
    csk.add(z);
    if (++klix == ldiff) klix = 0;
    if (l < w) {
      continue;
    }
    cminimizer = *std::min_element(winenc_v.begin(), winenc_v.end(), [](hmer_t lhs, hmer_t rhs) { return lhs.z < rhs.z; });
    if (is_wmin && cminimizer.x == cx_prev) continue;
    is_wmin = true;
    cx_prev = cminimizer.x;
    const uint32_t rix = lshf.compute_hash_bp(cminimizer.x);
    if (rix < fracth) {
      keys_v.push_back(pack_key(rix, lshf.drop_ppos_lr(cminimizer.y)));
    }
  }
}

QSeq::QSeq(const str& input, uint64_t bpmax_batch)
  : bpmax_batch(bpmax_batch)
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
  while ((ix < rbatch_size) && (nbases < bpmax_batch) && (cont_reading = (kseq_read(kseq) >= 0))) {
    nbases += kseq->seq.l;
    batch_v.push_back({kseq->name.s, kseq->seq.s});
    ix++;
  }
  return cont_reading;
}

bool QSeq::is_empty() { return batch_v.empty(); }

void QSeq::clear() { batch_v.clear(); }
