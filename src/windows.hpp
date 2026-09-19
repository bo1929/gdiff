#ifndef _WINDOWS_HPP
#define _WINDOWS_HPP

// Windows: the sampled stretches stored in a container (their metadata plus the Pool and
// Seq payloads), the per-window count aggregators the scan paths fill, and the sampling
// plan that picks which windows to visit.

#include "distance.hpp"
#include "types.hpp"
#include <algorithm>
#include <cstdint>
#include <random>

// A non-owning run of Ts: the mapped counterpart of vec<T>.
template<typename T>
struct slice_t
{
  const T* data = nullptr;
  uint64_t n = 0;
};

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

  [[nodiscard]] bool is_view() const noexcept { return hashes_view.data != nullptr; }
  [[nodiscard]] uint64_t size() const noexcept { return is_view() ? hashes_view.n : hashes_v.size(); }
  [[nodiscard]] const uint64_t* hashes_ptr() const noexcept { return is_view() ? hashes_view.data : hashes_v.data(); }
  [[nodiscard]] const uint16_t* win_ix_ptr() const noexcept { return is_view() ? win_ix_view.data : win_ix_v.data(); }
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

  [[nodiscard]] bool is_view() const noexcept { return base_off_view.data != nullptr; }
  [[nodiscard]] const uint64_t* packed_ptr() const noexcept { return is_view() ? packed_view.data : packed_v.data(); }
  [[nodiscard]] const uint64_t* nmask_ptr() const noexcept
  {
    return is_view() ? nmask_view.data : (nmask_v.empty() ? nullptr : nmask_v.data());
  }
  [[nodiscard]] const uint64_t* base_off_ptr() const noexcept { return is_view() ? base_off_view.data : base_off_v.data(); }
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

// Reject a bin shift that cannot be represented, or that quantises wider than the window.
bool validate_binning(uint64_t bin_shift, uint64_t tau);

// Total hits accumulated in a per-window histogram (indices 0..hdist_th).
inline uint64_t hist_total(const uint64_t* hist, uint32_t hdist_th) noexcept
{
  uint64_t t = 0;
  for (uint32_t hd = 0; hd <= hdist_th; ++hd)
    t += hist[hd];
  return t;
}

// Per-window HD counts (canonical / one strand). No fw/rc split.
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

  // scan_mers_range aggregator (is_rc ignored; canonical scan never sets it)
  inline void operator()(uint64_t /*bin*/, uint32_t hdist, bool /*is_rc*/) noexcept
  {
    if (hdist <= hdist_th)
      ++hist_v[hdist];
    else
      ++u;
  }

  inline void skip_mer(uint64_t /*bin*/) const noexcept {}
};

// Per-window HD counts with explicit fw/rc halves for strand-aware scans.
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

  // scan_mers_range aggregator
  inline void operator()(uint64_t /*bin*/, uint32_t hdist, bool is_rc) noexcept
  {
    if (hdist <= hdist_th)
      ++(is_rc ? hist_rc_v[hdist] : hist_fw_v[hdist]);
    else
      ++(is_rc ? u_rc : u_fw);
  }

  inline void skip_mer(uint64_t /*bin*/) const noexcept {}
};

#endif
