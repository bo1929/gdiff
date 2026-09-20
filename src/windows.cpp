#include "windows.hpp"

#include "random.hpp"
#include <utility>
window_plan_t make_window_plan(const vec<uint64_t>& source_lens_v,
                               uint64_t k,
                               uint64_t tau,
                               uint64_t bin_shift,
                               uint64_t sample_size,
                               std::mt19937& rng)
{
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  const uint64_t nwinmers = window_nmers(tau, bin_shift);
  const uint64_t xtau = nwinmers + k - 1;

  // Cumulative eligible positions, so one global draw can be split back per source.
  uint64_t total_npos = 0;
  vec<uint64_t> lenc_v(source_lens_v.size());
  for (size_t i = 0; i < source_lens_v.size(); ++i) {
    const uint64_t len = source_lens_v[i];
    if (len >= xtau) total_npos += (len - k + 1 - nwinmers) / bin_size + 1;
    lenc_v[i] = total_npos;
  }

  window_plan_t plan;
  plan.nwinmers = nwinmers;
  if (total_npos == 0) return plan;

  const vec<uint64_t> positions_v = sample_coords(total_npos, sample_size, rng);
  plan.sources_v.reserve(source_lens_v.size());
  size_t pidx = 0;
  for (size_t i = 0; i < source_lens_v.size() && pidx < positions_v.size(); ++i) {
    const uint64_t roff = i == 0 ? 0 : lenc_v[i - 1];
    if (lenc_v[i] == roff) continue; // nothing eligible in this source
    const uint64_t enmers = source_lens_v[i] - k + 1;
    window_plan_t::source_t entry;
    entry.bix = i;
    entry.enmers = enmers;
    while (pidx < positions_v.size() && positions_v[pidx] < lenc_v[i]) {
      entry.starts_v.push_back(positions_v[pidx] - roff);
      ++pidx;
    }
    if (entry.starts_v.empty()) continue;
    plan.sources_v.push_back(std::move(entry));
  }
  return plan;
}

// Expand window wi into cseq as ACGT characters plus 'N' where masked.
void seq_pack_t::unpack(uint64_t wi, str& cseq) const
{
  static const char bases[4] = {'A', 'C', 'G', 'T'};
  const uint64_t* off = base_off_ptr();
  const uint64_t* packed_p = packed_ptr();
  const uint64_t* nmask_p = nmask_ptr();
  const uint64_t b0 = off[wi];
  const uint64_t b1 = off[wi + 1];
  cseq.resize(static_cast<size_t>(b1 - b0));
  for (uint64_t b = b0; b < b1; ++b) {
    if (nmask_p && ((nmask_p[b >> 6] >> (b & 63)) & 1ull)) {
      cseq[static_cast<size_t>(b - b0)] = 'N';
      continue;
    }
    const uint64_t code = (packed_p[b >> 5] >> (2 * (b & 31))) & 3ull;
    cseq[static_cast<size_t>(b - b0)] = bases[code];
  }
}
