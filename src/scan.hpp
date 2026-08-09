#ifndef _SCAN_HPP
#define _SCAN_HPP

#include "dim.hpp"
#include "stils.hpp"
#include "sketch.hpp"
#include <cstdint>
#include <limits>

inline uint32_t bucket_hdist_min(const enc_t* ix1, const enc_t* ix2, const enc_t enc)
{
  uint32_t hdist_min = std::numeric_limits<uint32_t>::max();
  for (; ix1 < ix2; ++ix1) {
    const uint32_t hdist = popcount_lr32((*ix1) ^ enc);
    hdist_min = hdist < hdist_min ? hdist : hdist_min;
  }
  return hdist_min;
}

struct scan_ctx_t
{
  const Sketch* sketch;
  const LSHF* lshf;
  const SFHM* sfhm;
  uint32_t k;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint64_t bin_shift;
  uint32_t hdist_th;
  uint32_t tolerance_len; // N-runs longer than this break intervals
};

inline scan_ctx_t make_scan_ctx(const Sketch& sketch, uint64_t bin_shift, uint32_t hdist_th)
{
  const lshf_sptr_t lshf = sketch.get_lshf_sptr();
  const uint32_t k = lshf->get_k();
  const uint64_t u64max = std::numeric_limits<uint64_t>::max();
  return {&sketch,
          lshf.get(),
          sketch.get_sfhm_sptr().get(),
          k,
          u64max >> ((32 - k) * 2),
          ((u64max >> (64 - k)) << 32) + ((u64max << 32) >> (64 - k)),
          bin_shift,
          hdist_th,
          k};
}

constexpr size_t scan_block = 64;

// C <- canonical or per-strand
// Agg <- function we use for aggregation
template<bool C, typename Agg>
inline void scan_mers_range(const scan_ctx_t& ctx, const char* cseq, const uint64_t j0, const uint64_t j1, Agg&& agg)
{
  const uint32_t k = ctx.k;
  const uint64_t i1 = j1 + k - 1;
  const LSHF* lshf = ctx.lshf;
  const Sketch* sketch = ctx.sketch;
  const SFHM* sfhm = ctx.sfhm;
  uint64_t l = 0;
  uint64_t n_run = 0; // consecutive ambiguous bases so far
  uint64_t enc_lr = 0, enc_bp = 0;
  size_t n = 0;
  uint64_t b_bin[scan_block];
  uint32_t b_off[2][scan_block];
  enc_t b_enc[2][scan_block];

  auto flush = [&]() {
    // Phase A: prefetch the bucket-boundary lines for the whole block.
    for (size_t s = 0; s < n; ++s) {
      if (b_off[0][s] != Sketch::INVALID_BIX) sfhm->prefetch_inc(b_off[0][s]);
      if constexpr (C) {
        if (b_off[1][s] != Sketch::INVALID_BIX) sfhm->prefetch_inc(b_off[1][s]);
      }
    }
    // Phase B: resolve bucket bounds (now cache-warm) and prefetch enc lines.
    const enc_t* beg[2][scan_block];
    const enc_t* end[2][scan_block];
    for (size_t s = 0; s < n; ++s) {
      const uint32_t off0 = b_off[0][s];
      if (off0 == Sketch::INVALID_BIX) {
        beg[0][s] = end[0][s] = nullptr;
      } else {
        beg[0][s] = sfhm->bucket_ptr_start(off0);
        end[0][s] = sfhm->bucket_ptr_next(off0);
        if (beg[0][s] < end[0][s]) __builtin_prefetch(beg[0][s], 0, 0);
      }
      if constexpr (C) {
        const uint32_t off1 = b_off[1][s];
        if (off1 == Sketch::INVALID_BIX) {
          beg[1][s] = end[1][s] = nullptr;
        } else {
          beg[1][s] = sfhm->bucket_ptr_start(off1);
          end[1][s] = sfhm->bucket_ptr_next(off1);
          if (beg[1][s] < end[1][s]) __builtin_prefetch(beg[1][s], 0, 0);
        }
      }
    }
    // Phase C: scan buckets and agg observed k-mers (hits and misses).
    for (size_t s = 0; s < n; ++s) {
      if (beg[0][s] != nullptr) {
        const uint32_t hdist = bucket_hdist_min(beg[0][s], end[0][s], b_enc[0][s]);
        agg(b_bin[s], hdist, false);
      }
      if constexpr (C) {
        if (beg[1][s] != nullptr) {
          const uint32_t hdist = bucket_hdist_min(beg[1][s], end[1][s], b_enc[1][s]);
          agg(b_bin[s], hdist, true);
        }
      }
    }
    n = 0;
  };

  for (uint64_t i = j0; i < i1; ++i) {
    if (__builtin_expect(SEQ_NT4_TABLE[static_cast<uint8_t>(cseq[i])] >= 4, 0)) {
      l = 0;
      if (++n_run > ctx.tolerance_len) agg.skip_mer(i >> ctx.bin_shift);
      continue;
    }
    n_run = 0;
    ++l;
    if (l < k) continue;
    const uint64_t j = i - k + 1;
    if (l == k) {
      compute_encoding(cseq + j, cseq + i + 1, enc_lr, enc_bp);
    } else {
      update_encoding(cseq + i, enc_lr, enc_bp);
    }
    enc_bp &= ctx.mask_bp;
    enc_lr &= ctx.mask_lr;
    if (__builtin_expect(j >= j1, 0)) break;
    const uint64_t rc_bp = revcomp_bp64(enc_bp, k);
    if constexpr (C) {
      b_off[0][n] = sketch->validate_bucket_ix(lshf->compute_hash_bp(enc_bp));
      b_off[1][n] = sketch->validate_bucket_ix(lshf->compute_hash_bp(rc_bp));
      b_enc[0][n] = lshf->drop_ppos_lr(enc_lr);
      b_enc[1][n] = lshf->drop_ppos_lr(bp64_to_lr64(rc_bp));
    } else {
      if (rc_bp < enc_bp) {
        b_off[0][n] = sketch->validate_bucket_ix(lshf->compute_hash_bp(enc_bp));
        b_enc[0][n] = lshf->drop_ppos_lr(enc_lr);
      } else {
        b_off[0][n] = sketch->validate_bucket_ix(lshf->compute_hash_bp(rc_bp));
        b_enc[0][n] = lshf->drop_ppos_lr(bp64_to_lr64(rc_bp));
      }
    }
    b_bin[n] = j >> ctx.bin_shift;
    if (++n == scan_block) flush();
  }
  if (n) flush();
}

// Adapts DIM<T>::aggregate_mer to the scan_mers_range call signature.
template<typename T>
struct dim_agg_t
{
  DIM<T>& fw;
  DIM<T>* rc; // null in canonical mode
  inline void operator()(uint64_t bin, uint32_t hd, bool is_rc) const { (is_rc ? *rc : fw).aggregate_mer(hd, bin); }
  inline void skip_mer(uint64_t bin) const
  {
    fw.skip_mer(bin);
    if (rc) rc->skip_mer(bin);
  }
};

#endif
