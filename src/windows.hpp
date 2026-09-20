#ifndef _WINDOWS_HPP
#define _WINDOWS_HPP

#include "types.hpp"
#include <algorithm>
#include <cstdint>
#include <random>

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

struct seq_pack_t
{
  vec<uint64_t> packed_v;
  vec<uint64_t> nmask_v;
  vec<uint64_t> boffset_v; // nwins + 1 entries
  slice_t<uint64_t> packed_view;
  slice_t<uint64_t> nmask_view;
  slice_t<uint64_t> boffset_view;

  [[nodiscard]] const uint64_t* packed_ptr() const noexcept { return pick_slice(packed_view, packed_v).data; }
  [[nodiscard]] const uint64_t* nmask_ptr() const noexcept
  {
    if (nmask_view.data != nullptr) return nmask_view.data;
    return nmask_v.empty() ? nullptr : nmask_v.data();
  }
  [[nodiscard]] const uint64_t* boffset_ptr() const noexcept { return pick_slice(boffset_view, boffset_v).data; }
  void unpack(uint64_t wi, str& cseq) const;
};

struct window_sample_t
{
  vec<window_t> wins_v;
  hash_pool_t pool_fw;
  hash_pool_t pool_rc; // empty when canonical
  seq_pack_t packs;
};

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

inline uint64_t get_nwinmers(uint64_t tau, uint64_t bin_shift) noexcept
{
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  return std::max<uint64_t>(1, (tau + bin_size - 1) >> bin_shift) << bin_shift;
}

window_plan_t make_window_plan(const vec<uint64_t>& slen_v,
                               uint64_t k,
                               uint64_t tau,
                               uint64_t bin_shift,
                               uint64_t sample_size,
                               std::mt19937& rng);

#endif
