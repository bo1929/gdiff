#ifndef _SCAN_HPP
#define _SCAN_HPP

// Shared pipelined k-mer scanner used by dist and detect. Hoisted from
// dist.cpp so both subcommands scan sequences with the same prefetch-blocked,
// branchless implementation.

#include "enc.hpp"
#include "maptils.hpp"
#include "sketch.hpp"
#include "types.hpp"
#include <cstdint>
#include <limits>

// Branchless minimum Hamming distance of enc against a bucket; identical
// result to the early-exit scan but without data-dependent branches.
inline uint32_t bucket_hdist_min(const enc_t* ix1, const enc_t* ix2, const enc_t enc)
{
  uint32_t hmin = std::numeric_limits<uint32_t>::max();
  for (; ix1 < ix2; ++ix1) {
    const uint32_t hd = popcount_lr32((*ix1) ^ enc);
    hmin = hd < hmin ? hd : hmin;
  }
  return hmin;
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
};

// Builds the scanner context for a sketch, including the encoding masks:
// mask_bp keeps the 2k bp bits, mask_lr keeps the k low-resolution bits.
inline scan_ctx_t make_scan_ctx(const Sketch& sketch, uint64_t bin_shift, uint32_t hdist_th)
{
  const lshf_sptr_t lshf = sketch.get_lshf_sptr();
  const uint32_t k = lshf->get_k();
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  return {&sketch,
          lshf.get(),
          sketch.get_sfhm_sptr().get(),
          k,
          u64m >> ((32 - k) * 2),
          ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k)),
          bin_shift,
          hdist_th};
}

// Block size of the software pipeline: hashes/offsets for scan_blk mers are
// resolved first while their prefetches are in flight, hiding memory latency.
constexpr size_t scan_blk = 64;

// Scans mers whose start position j lies in [j0, j1) (j1 <= enmers).
// agg(bin_j, hdist, is_rc) is invoked for every observed k-mer (valid LSH bucket).
// STRAND_AWARE: query both strands; otherwise query only the canonical one.
template<bool STRAND_AWARE, typename Agg>
inline void scan_mers_range(const scan_ctx_t& ctx, const char* cseq, const uint64_t j0, const uint64_t j1, Agg&& agg)
{
  const uint32_t k = ctx.k;
  const uint64_t i1 = j1 + k - 1;
  const LSHF* lshf = ctx.lshf;
  const Sketch* sketch = ctx.sketch;
  const SFHM* sfhm = ctx.sfhm;
  uint64_t l = 0;
  uint64_t enc_lr = 0, enc_bp = 0;
  size_t n = 0;
  uint64_t b_bin[scan_blk];
  uint32_t b_off[2][scan_blk];
  enc_t b_enc[2][scan_blk];

  auto flush = [&]() {
    // Phase A: prefetch the bucket-boundary lines for the whole block.
    for (size_t s = 0; s < n; ++s) {
      if (b_off[0][s] != Sketch::OFF_INVALID) sfhm->prefetch_inc(b_off[0][s]);
      if constexpr (STRAND_AWARE) {
        if (b_off[1][s] != Sketch::OFF_INVALID) sfhm->prefetch_inc(b_off[1][s]);
      }
    }
    // Phase B: resolve bucket bounds (now cache-warm) and prefetch enc lines.
    const enc_t* beg[2][scan_blk];
    const enc_t* end[2][scan_blk];
    for (size_t s = 0; s < n; ++s) {
      const uint32_t off0 = b_off[0][s];
      if (off0 == Sketch::OFF_INVALID) {
        beg[0][s] = end[0][s] = nullptr;
      } else {
        beg[0][s] = sfhm->bucket_ptr_start(off0);
        end[0][s] = sfhm->bucket_ptr_next(off0);
        if (beg[0][s] < end[0][s]) __builtin_prefetch(beg[0][s], 0, 0);
      }
      if constexpr (STRAND_AWARE) {
        const uint32_t off1 = b_off[1][s];
        if (off1 == Sketch::OFF_INVALID) {
          beg[1][s] = end[1][s] = nullptr;
        } else {
          beg[1][s] = sfhm->bucket_ptr_start(off1);
          end[1][s] = sfhm->bucket_ptr_next(off1);
          if (beg[1][s] < end[1][s]) __builtin_prefetch(beg[1][s], 0, 0);
        }
      }
    }
    // Phase C: scan buckets and aggregate observed k-mers (hits and misses).
    for (size_t s = 0; s < n; ++s) {
      if (beg[0][s] != nullptr) {
        const uint32_t hd = bucket_hdist_min(beg[0][s], end[0][s], b_enc[0][s]);
        agg(b_bin[s], hd, false);
      }
      if constexpr (STRAND_AWARE) {
        if (beg[1][s] != nullptr) {
          const uint32_t hd = bucket_hdist_min(beg[1][s], end[1][s], b_enc[1][s]);
          agg(b_bin[s], hd, true);
        }
      }
    }
    n = 0;
  };

  for (uint64_t i = j0; i < i1; ++i) {
    if (__builtin_expect(SEQ_NT4_TABLE[static_cast<uint8_t>(cseq[i])] >= 4, 0)) {
      l = 0;
      continue;
    }
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
    if constexpr (STRAND_AWARE) {
      b_off[0][n] = sketch->partial_offset(lshf->compute_hash(enc_bp));
      b_off[1][n] = sketch->partial_offset(lshf->compute_hash(rc_bp));
      b_enc[0][n] = lshf->drop_ppos_lr(enc_lr);
      b_enc[1][n] = lshf->drop_ppos_lr(bp64_to_lr64(rc_bp));
    } else {
      if (rc_bp < enc_bp) {
        b_off[0][n] = sketch->partial_offset(lshf->compute_hash(enc_bp));
        b_enc[0][n] = lshf->drop_ppos_lr(enc_lr);
      } else {
        b_off[0][n] = sketch->partial_offset(lshf->compute_hash(rc_bp));
        b_enc[0][n] = lshf->drop_ppos_lr(bp64_to_lr64(rc_bp));
      }
    }
    b_bin[n] = j >> ctx.bin_shift;
    if (++n == scan_blk) flush();
  }
  if (n) flush();
}

// Per-window aggregator: counts HD histogram and explicit misses for one
// sampled window, per strand.
struct window_agg_t
{
  uint64_t* fw; // hdist_th + 1 hit counts per strand (rc unused in canonical mode)
  uint64_t* rc;
  uint64_t& u_fw; // miss counters
  uint64_t& u_rc;
  uint32_t hdist_th;
  inline void operator()(uint64_t, uint32_t hd, bool is_rc) const
  {
    if (hd <= hdist_th)
      ++(is_rc ? rc[hd] : fw[hd]);
    else
      ++(is_rc ? u_rc : u_fw);
  }
};

#endif
