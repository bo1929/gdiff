#ifndef _SCAN_HPP
#define _SCAN_HPP

#include "intext.hpp"
#include "sketch.hpp"
#include <algorithm>
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

// Scalar LLH for one reference/threshold pair; derivatives are not needed.
inline LLH<double> make_llhf(const Sketch& sketch, uint32_t hdist_th)
{
  return {sketch.get_k(), sketch.get_h(), sketch.get_rho(), hdist_th, 0.0, false};
}

struct scan_ctx_t
{
  const LSHF* lshf;
  const Buckets* buckets;
  uint32_t k;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint64_t bin_shift;
  uint32_t hdist_th;
  uint32_t tol_len; // N-runs longer than this break intervals
};

inline scan_ctx_t make_scan_ctx(const Sketch& sketch, uint64_t bin_shift, uint32_t hdist_th)
{
  const LSHF* lshf = sketch.get_lshf_sptr().get();
  const uint32_t k = lshf->get_k();
  const lsh_masks_t masks = get_lsh_masks(k);
  return {lshf, &sketch.get_buckets(), k, masks.bp, masks.lr, bin_shift, hdist_th, k};
}

// Canonical: one key per k-mer (true) vs strand-aware fw+rc (false). Agg: hit/miss aggregator.
template<bool Canonical, typename Agg>
inline void scan_mers_range(const scan_ctx_t& ctx, const char* cseq, uint64_t j0, uint64_t j1, Agg&& agg)
{
  const uint32_t k = ctx.k;
  const uint64_t i1 = j1 + k - 1;
  const LSHF* lshf = ctx.lshf;
  const Buckets* buckets = ctx.buckets;
  uint64_t l = 0;
  uint64_t n_skip = 0; // length of the current consecutive ambiguous base run
  uint64_t enc_lr = 0, enc_bp = 0;
  size_t n = 0;

  constexpr size_t bufmer_count = 64;
  uint64_t bufbin[bufmer_count];
  uint32_t bufbix[2][bufmer_count];
  enc_t bufenc[2][bufmer_count];

  auto flush = [&]() {
    // Phase A: prefetch the occupancy lines for the whole block.
    for (size_t s = 0; s < n; ++s) {
      buckets->prefetch(bufbix[0][s]);
      if constexpr (!Canonical) buckets->prefetch(bufbix[1][s]);
    }
    // Phase B: resolve bucket bounds (now cache-warm) and prefetch enc lines.
    const enc_t* ix_start[2][bufmer_count];
    const enc_t* ix_end[2][bufmer_count];
    for (size_t s = 0; s < n; ++s) {
      if (!buckets->range(bufbix[0][s], ix_start[0][s], ix_end[0][s])) {
        ix_start[0][s] = ix_end[0][s] = nullptr;
      } else {
        __builtin_prefetch(ix_start[0][s], 0, 0);
      }
      if constexpr (!Canonical) {
        if (!buckets->range(bufbix[1][s], ix_start[1][s], ix_end[1][s])) {
          ix_start[1][s] = ix_end[1][s] = nullptr;
        } else {
          __builtin_prefetch(ix_start[1][s], 0, 0);
        }
      }
    }
    // Phase C: scan buckets and agg observed k-mers (hits and misses).
    for (size_t s = 0; s < n; ++s) {
      if (ix_start[0][s] != nullptr) {
        const uint32_t hdist = bucket_hdist_min(ix_start[0][s], ix_end[0][s], bufenc[0][s]);
        agg(bufbin[s], hdist, false);
      } else {
        // Out-of-range (--frac) or empty bucket: not in the sketch, a miss.
        agg(bufbin[s], ctx.hdist_th + 1, false);
      }
      if constexpr (!Canonical) {
        if (ix_start[1][s] != nullptr) {
          const uint32_t hdist = bucket_hdist_min(ix_start[1][s], ix_end[1][s], bufenc[1][s]);
          agg(bufbin[s], hdist, true);
        } else {
          agg(bufbin[s], ctx.hdist_th + 1, true);
        }
      }
    }
    n = 0;
  };

  for (uint64_t i = j0; i < i1; ++i) {
    if (__builtin_expect(SEQ_NT4_TABLE[static_cast<uint8_t>(cseq[i])] >= 4, 0)) {
      l = 0;
      if (++n_skip > ctx.tol_len) agg.skip_mer(i >> ctx.bin_shift);
      continue;
    }
    n_skip = 0;
    ++l;
    if (l < k) continue;
    // i + 1 - k: i - k + 1 wraps on the first window of a sequence.
    const uint64_t j = i + 1 - k;
    if (l == k) {
      compute_encoding(cseq + j, cseq + i + 1, enc_lr, enc_bp);
    } else {
      update_encoding(cseq + i, enc_lr, enc_bp);
    }
    enc_bp &= ctx.mask_bp;
    enc_lr &= ctx.mask_lr;
    if (__builtin_expect(j >= j1, 0)) break;
    const uint64_t rcenc_bp = revcomp_bp64(enc_bp, k);
    if constexpr (!Canonical) {
      bufbix[0][n] = lshf->compute_hash_bp(enc_bp);
      bufbix[1][n] = lshf->compute_hash_bp(rcenc_bp);
      bufenc[0][n] = lshf->drop_ppos_lr(enc_lr);
      bufenc[1][n] = lshf->drop_ppos_lr(bp64_to_lr64(rcenc_bp));
    } else if (rcenc_bp < enc_bp) {
      bufbix[0][n] = lshf->compute_hash_bp(enc_bp);
      bufenc[0][n] = lshf->drop_ppos_lr(enc_lr);
    } else {
      bufbix[0][n] = lshf->compute_hash_bp(rcenc_bp);
      bufenc[0][n] = lshf->drop_ppos_lr(bp64_to_lr64(rcenc_bp));
    }
    bufbin[n] = j >> ctx.bin_shift;
    if (++n == bufmer_count) flush();
  }
  if (n) flush();
}

// Total hits accumulated in a per-window histogram (indices 0..hdist_th).
inline uint64_t hist_total(const uint64_t* hist, uint32_t hdist_th) noexcept
{
  uint64_t t = 0;
  for (uint32_t hd = 0; hd <= hdist_th; ++hd)
    t += hist[hd];
  return t;
}

// Per-window HD counts (canonical - strand agnostic).
struct window_counts_t
{
  vec<uint64_t> hist_v;
  uint64_t u = 0;
  uint32_t hdist_th = 0;

  window_counts_t() = default;

  explicit window_counts_t(uint32_t hdist_th)
    : hdist_th(hdist_th)
  {
    hist_v.assign(hdist_bound + 1, 0);
  }

  void clear() noexcept
  {
    std::fill(hist_v.begin(), hist_v.end(), 0);
    u = 0;
  }

  uint64_t* hist() noexcept { return hist_v.data(); }
  const uint64_t* hist() const noexcept { return hist_v.data(); }

  inline void operator()(uint64_t /*bin*/, uint32_t hdist, bool /*is_rc*/) noexcept
  {
    if (hdist <= hdist_th)
      ++hist_v[hdist];
    else
      ++u;
  }

  inline void skip_mer(uint64_t /*bin*/) const noexcept {}
};

// Per-window HD counts with explicit fw/rc halves for strand-based scans.
struct swindow_counts_t
{
  vec<uint64_t> hist_fw_v;
  vec<uint64_t> hist_rc_v;
  uint64_t u_fw = 0;
  uint64_t u_rc = 0;
  uint32_t hdist_th = 0;

  swindow_counts_t() = default;

  explicit swindow_counts_t(uint32_t hdist_th)
    : hdist_th(hdist_th)
  {
    hist_fw_v.assign(hdist_bound + 1, 0);
    hist_rc_v.assign(hdist_bound + 1, 0);
  }

  void clear() noexcept
  {
    std::fill(hist_fw_v.begin(), hist_fw_v.end(), 0);
    std::fill(hist_rc_v.begin(), hist_rc_v.end(), 0);
    u_fw = 0;
    u_rc = 0;
  }

  uint64_t* hist_fw() noexcept { return hist_fw_v.data(); }
  uint64_t* hist_rc() noexcept { return hist_rc_v.data(); }
  const uint64_t* hist_fw() const noexcept { return hist_fw_v.data(); }
  const uint64_t* hist_rc() const noexcept { return hist_rc_v.data(); }

  inline void operator()(uint64_t /*bin*/, uint32_t hdist, bool is_rc) noexcept
  {
    if (hdist <= hdist_th)
      ++(is_rc ? hist_rc_v[hdist] : hist_fw_v[hdist]);
    else
      ++(is_rc ? u_rc : u_fw);
  }

  inline void skip_mer(uint64_t /*bin*/) const noexcept {}
};

template<typename T>
struct intext_agg_t
{
  IntExt<T>& fw;
  IntExt<T>* rc;
  inline void operator()(uint64_t bin, uint32_t hd, bool is_rc) const { (is_rc ? *rc : fw).aggregate_mer(hd, bin); }
  inline void skip_mer(uint64_t bin) const
  {
    fw.skip_mer(bin);
    if (rc) rc->skip_mer(bin);
  }
};

#endif
