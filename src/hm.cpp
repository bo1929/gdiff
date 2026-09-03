#include "hm.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include "rqseq.hpp"

SFHM::SFHM(const sdhm_sptr_t& source)
{
  nkmers = source->nkmers;
  nrows = source->enc_vvec.size();
  if (nkmers > static_cast<uint64_t>(std::numeric_limits<inc_t>::max())) {
    error_exit("There are more k-mers than maximum inc_t. Recompile with an appropriate type or reduce the sketch size.");
  }
  inc_v.resize(nrows);
  enc_v.reserve(nkmers);
  inc_t limit_inc = std::numeric_limits<inc_t>::max();
  inc_t cpinc;
  inc_t lix = 0;
  for (uint32_t rix = 0; rix < nrows; ++rix) {
    cpinc = std::min(limit_inc, static_cast<inc_t>(source->enc_vvec[rix].size()));
    enc_v.insert(enc_v.end(), source->enc_vvec[rix].begin(), source->enc_vvec[rix].begin() + cpinc);
    lix += cpinc;
    inc_v[rix] = lix;
    source->enc_vvec[rix].clear();
  }
}

SFHM::~SFHM()
{
  inc_v.clear();
  enc_v.clear();
}

void SFHM::load(std::ifstream& sketch_stream)
{
  sketch_stream.read(reinterpret_cast<char*>(&nkmers), sizeof(uint64_t));
  if (nkmers > static_cast<uint64_t>(1) << 40) {
    error_exit("Corrupt sketch file, or a compatibility issue!?!");
  }
  enc_v.resize(nkmers);
  sketch_stream.read(reinterpret_cast<char*>(enc_v.data()), nkmers * sizeof(enc_t));
  assert(nkmers == enc_v.size());
  sketch_stream.read(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));
  inc_v.resize(nrows);
  sketch_stream.read(reinterpret_cast<char*>(inc_v.data()), nrows * sizeof(inc_t));
  assert(nrows == inc_v.size());
}

void SFHM::load_mem(const char*& p, const char* end)
{
  auto need = [&](size_t n, const char* what) {
    if (p == nullptr || end == nullptr || static_cast<size_t>(end - p) < n) {
      error_exit(std::string("Truncated SFHM while reading ") + what);
    }
  };
  need(sizeof(uint64_t), "nkmers");
  std::memcpy(&nkmers, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  if (nkmers > static_cast<uint64_t>(1) << 40) {
    error_exit("Corrupt sketch file, or a compatibility issue!?!");
  }
  const size_t enc_bytes = static_cast<size_t>(nkmers) * sizeof(enc_t);
  need(enc_bytes + sizeof(uint32_t), "encodings");
  enc_v.resize(nkmers);
  if (nkmers) std::memcpy(enc_v.data(), p, enc_bytes);
  p += enc_bytes;
  std::memcpy(&nrows, p, sizeof(uint32_t));
  p += sizeof(uint32_t);
  const size_t inc_bytes = static_cast<size_t>(nrows) * sizeof(inc_t);
  need(inc_bytes, "increments");
  inc_v.resize(nrows);
  if (nrows) std::memcpy(inc_v.data(), p, inc_bytes);
  p += inc_bytes;
}

void SFHM::view_mem(const char*& p, const char* end)
{
  auto need = [&](size_t n, const char* what) {
    if (p == nullptr || end == nullptr || static_cast<size_t>(end - p) < n) {
      error_exit(std::string("Truncated SFHM while viewing ") + what);
    }
  };
  need(sizeof(uint64_t), "nkmers");
  std::memcpy(&nkmers, p, sizeof(uint64_t));
  p += sizeof(uint64_t);
  if (nkmers > static_cast<uint64_t>(1) << 40) {
    error_exit("Corrupt sketch file, or a compatibility issue!?!");
  }
  const size_t enc_bytes = static_cast<size_t>(nkmers) * sizeof(enc_t);
  need(enc_bytes + sizeof(uint32_t), "encodings");
  enc_view = reinterpret_cast<const enc_t*>(p);
  p += enc_bytes;
  std::memcpy(&nrows, p, sizeof(uint32_t));
  p += sizeof(uint32_t);
  const size_t inc_bytes = static_cast<size_t>(nrows) * sizeof(inc_t);
  need(inc_bytes, "increments");
  inc_view = reinterpret_cast<const inc_t*>(p);
  p += inc_bytes;
  view_mode_ = true;
  // Drop any pre-existing owned payload; the view is authoritative.
  enc_v.clear();
  inc_v.clear();
}

void SFHM::save(std::ostream& sketch_stream)
{
  const enc_t* enc_ptr = enc_view ? enc_view : enc_v.data();
  const inc_t* inc_ptr = inc_view ? inc_view : inc_v.data();
  sketch_stream.write(reinterpret_cast<const char*>(&nkmers), sizeof(uint64_t));
  sketch_stream.write(reinterpret_cast<const char*>(enc_ptr), sizeof(enc_t) * nkmers);
  sketch_stream.write(reinterpret_cast<const char*>(&nrows), sizeof(uint32_t));
  sketch_stream.write(reinterpret_cast<const char*>(inc_ptr), sizeof(inc_t) * nrows);
}

std::vector<enc_t>::const_iterator SFHM::bucket_iter_start(uint32_t rix)
{
  error_exit("bucket_iter_start on a borrowed-view SFHM");
}

std::vector<enc_t>::const_iterator SFHM::bucket_iter_next(uint32_t rix)
{
  error_exit("bucket_iter_next on a borrowed-view SFHM");
}

const enc_t* SFHM::bucket_ptr_start(uint32_t rix) const noexcept
{
  if (view_mode_) {
    return enc_view + ((rix != 0 && rix <= nrows) ? inc_view[rix - 1] : 0);
  }
  return enc_v.data() + ((rix != 0 && rix <= inc_v.size()) ? inc_v[rix - 1] : 0);
}

const enc_t* SFHM::bucket_ptr_next(uint32_t rix) const noexcept
{
  if (view_mode_) {
    return enc_view + (rix < nrows ? inc_view[rix] : nkmers);
  }
  return enc_v.data() + (rix < inc_v.size() ? inc_v[rix] : nkmers);
}

void SFHM::prefetch_inc(uint32_t rix) const noexcept
{
  if (rix > 0 && rix <= (view_mode_ ? nrows : inc_v.size())) {
    const inc_t* base = view_mode_ ? inc_view : inc_v.data();
    __builtin_prefetch(&base[rix - 1], 0, 1);
  }
}

void SFHM::prefetch_enc(uint32_t rix) const noexcept
{
  const enc_t* start = bucket_ptr_start(rix);
  const enc_t* end = bucket_ptr_next(rix);
  if (start < end) {
    __builtin_prefetch(start, 0, 0);
  }
}

void SFHM::fill_nonempty_bitmap(uint64_t* bits, uint32_t nbits) const noexcept
{
  const uint32_t n = std::min(nbits, nrows);
  const inc_t* inc_ptr = view_mode_ ? inc_view : inc_v.data();
  inc_t prev = 0;
  for (uint32_t i = 0; i < n; ++i) {
    const inc_t cur = inc_ptr[i];
    if (cur > prev) bits[i >> 6] |= (uint64_t(1) << (i & 63));
    prev = cur;
  }
}

void SDHM::sort_columns()
{
  for (uint32_t i = 0; i < enc_vvec.size(); ++i) {
    if (!enc_vvec[i].empty()) {
      std::sort(enc_vvec[i].begin(), enc_vvec[i].end());
    }
  }
}

void SDHM::make_unique()
{
  nkmers = 0;
  for (uint32_t i = 0; i < enc_vvec.size(); ++i) {
    if (!enc_vvec[i].empty()) {
      enc_vvec[i].erase(std::unique(enc_vvec[i].begin(), enc_vvec[i].end()), enc_vvec[i].end());
    }
    nkmers += enc_vvec[i].size();
  }
}

void SDHM::fill_table(uint32_t nrows, const rseq_sptr_t& rs)
{
  enc_vvec.resize(nrows);
  while (rs->read_next_seq()) {
    if (rs->set_curr_seq()) {
      rs->extract_mers(enc_vvec);
    }
  }
  sort_columns();
  make_unique();
}

uint64_t SDHM::get_nmers() const { return nkmers; }
