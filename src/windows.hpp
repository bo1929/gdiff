#ifndef _WINDOWS_HPP
#define _WINDOWS_HPP

// Windows: the sampled stretches stored in a container (metadata plus the Pool and Seq
// payloads) and the sampling plan that picks which windows to visit.

#include "types.hpp"
#include <algorithm>
#include <cstdint>
#include <random>

// Window payload arrays are owned while a container is built and become mapped views once it
// is loaded; `pick_slice` resolves either representation.
template<typename T>
[[nodiscard]] inline slice_t<T> pick_slice(const slice_t<T>& view, const vec<T>& owned) noexcept
{
  return view.data ? view : slice_t<T>{owned.data(), owned.size()};
}

// Window metadata; the payload arrays alongside share its index.
struct window_t
{
  str qid;
  uint64_t start = 0;     // 0-based k-mer start in the source sequence
  uint64_t end = 0;       // exclusive k-mer end
  uint32_t nvalid_fw = 0; // valid (non-ambiguous) k-mers in the window
  uint32_t nvalid_rc = 0; // same on the reverse strand; 0 when canonical
};

// Window k-mers sorted by packed key; owned while building, viewed once mapped.
struct hash_pool_t
{
  vec<uint64_t> hashes_v;
  vec<uint16_t> win_ix_v;
  slice_t<uint64_t> hashes_view;
  slice_t<uint16_t> win_ix_view;

  void clear()
  {
    hashes_v.clear();
    win_ix_v.clear();
    hashes_view = {};
    win_ix_view = {};
  }

  [[nodiscard]] uint64_t size() const noexcept { return pick_slice(hashes_view, hashes_v).n; }
  [[nodiscard]] const uint64_t* hashes_ptr() const noexcept { return pick_slice(hashes_view, hashes_v).data; }
  [[nodiscard]] const uint16_t* win_ix_ptr() const noexcept { return pick_slice(win_ix_view, win_ix_v).data; }
};

// Concatenated 2-bit bases; base_off_v[i] is window i's bit-pair offset.
struct seq_pack_t
{
  vec<uint64_t> packed_v;
  vec<uint64_t> nmask_v;
  vec<uint64_t> base_off_v; // nwins + 1 entries
  slice_t<uint64_t> packed_view;
  slice_t<uint64_t> nmask_view;
  slice_t<uint64_t> base_off_view;

  [[nodiscard]] const uint64_t* packed_ptr() const noexcept { return pick_slice(packed_view, packed_v).data; }
  // The N mask is optional: a container with no masked base stores none at all.
  [[nodiscard]] const uint64_t* nmask_ptr() const noexcept
  {
    if (nmask_view.data != nullptr) return nmask_view.data;
    return nmask_v.empty() ? nullptr : nmask_v.data();
  }
  [[nodiscard]] const uint64_t* base_off_ptr() const noexcept { return pick_slice(base_off_view, base_off_v).data; }
  // Expand window wi into cseq as ACGT characters plus 'N' where masked.
  void unpack(uint64_t wi, str& cseq) const;
};

// All sampled windows of one sketch plus their payload, in one of the two representations.
struct window_sample_t
{
  vec<window_t> wins_v;
  hash_pool_t pool_fw;
  hash_pool_t pool_rc; // empty when canonical
  seq_pack_t packs;
};

// Which windows to evaluate, sampled from a set of sources: one entry per source that got a
// window, in input order. Starts are 0-based k-mer positions, ascending.
struct window_plan_t
{
  struct source_t
  {
    uint64_t bix = 0; // index into the caller's source list
    uint64_t enmers = 0;
    vec<uint64_t> starts_v;
  };

  vec<source_t> sources_v;
  uint64_t nwinmers = 0; // the window span in k-mers, rounded up to whole bins

  [[nodiscard]] uint64_t nwins() const
  {
    uint64_t n = 0;
    for (const source_t& src : sources_v)
      n += src.starts_v.size();
    return n;
  }
};

// The k-mer span of a sampled window: tau rounded up to whole bins of 2^bin_shift.
inline uint64_t window_nmers(uint64_t tau, uint64_t bin_shift) noexcept
{
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  return std::max<uint64_t>(1, (tau + bin_size - 1) >> bin_shift) << bin_shift;
}

// Sample windows of tau k-mers over `source_lens_v`, binned by 2^bin_shift, drawing at most
// `sample_size` starts in total (one source at a time when source_lens_v holds a single
// length). Sources too short for a full window contribute nothing. The draw is part of the
// reproducible output: neither the count nor the order of rng draws may change.
window_plan_t make_window_plan(const vec<uint64_t>& source_lens_v,
                               uint64_t k,
                               uint64_t tau,
                               uint64_t bin_shift,
                               uint64_t sample_size,
                               std::mt19937& rng);

#endif
