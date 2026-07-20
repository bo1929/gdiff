#include "sketch.hpp"

#include <utility>

Sketch::Sketch(std::filesystem::path sketch_path)
  : sketch_path(std::move(sketch_path))
{
}

void Sketch::load_from_offset(std::ifstream& stream, uint64_t offset)
{
  if (offset > 0) {
    stream.seekg(offset);
  }

  uint64_t rid_len;
  stream.read(reinterpret_cast<char*>(&rid_len), sizeof(uint64_t));
  rid.resize(rid_len);
  stream.read(&rid[0], rid_len);
  stream.read(reinterpret_cast<char*>(&timestamp), sizeof(uint64_t));

  stream.read(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(&canonical), sizeof(bool));
  stream.read(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));

  vec<uint8_t> ppos_v(h), npos_v(k - h);
  stream.read(reinterpret_cast<char*>(ppos_v.data()), h * sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(npos_v.data()), (k - h) * sizeof(uint8_t));

  lshf = std::make_shared<LSHF>(ppos_v, npos_v);

  stream.read(reinterpret_cast<char*>(&rho), sizeof(double));

  sfhm = std::make_shared<SFHM>();
  sfhm->load(stream);

  check_fstream(stream, "Failed to read the sketch file!", sketch_path);
}

void Sketch::seek_past(std::ifstream& stream)
{
  uint64_t rid_len;
  stream.read(reinterpret_cast<char*>(&rid_len), sizeof(uint64_t));
  // Skip rid (rid_len bytes) + timestamp (8 bytes)
  stream.seekg(static_cast<std::streamoff>(rid_len) + static_cast<std::streamoff>(sizeof(uint64_t)), std::ios::cur);

  uint8_t k, h;
  stream.read(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  // skip: w (1), h is read next, then canonical (1), nrows (4) = 6 bytes
  stream.seekg(sizeof(uint8_t), std::ios::cur);
  stream.read(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  // skip: canonical (1) + nrows (4) = 5 bytes
  stream.seekg(sizeof(bool) + sizeof(uint32_t), std::ios::cur);
  // skip: ppos_v (h bytes) + npos_v ((k-h) bytes) = k bytes total
  stream.seekg(static_cast<std::streamoff>(k), std::ios::cur);
  // skip: rho (8 bytes)
  stream.seekg(static_cast<std::streamoff>(sizeof(double)), std::ios::cur);

  uint64_t nkmers;
  stream.read(reinterpret_cast<char*>(&nkmers), sizeof(uint64_t));
  stream.seekg(static_cast<std::streamoff>(nkmers) * static_cast<std::streamoff>(sizeof(enc_t)), std::ios::cur);
  uint32_t sfhm_nrows;
  stream.read(reinterpret_cast<char*>(&sfhm_nrows), sizeof(uint32_t));
  stream.seekg(static_cast<std::streamoff>(sfhm_nrows) * static_cast<std::streamoff>(sizeof(inc_t)), std::ios::cur);
}

void Sketch::make_rho_partial()
{
  // The sketch keeps the prefix [0, nrows) of the 2^(2h) LSH space.
  rho *= static_cast<double>(nrows) / static_cast<double>(uint64_t(1) << (2 * h));
}

sfhm_sptr_t Sketch::get_sfhm_sptr() { return sfhm; }

lshf_sptr_t Sketch::get_lshf() { return lshf; }

double Sketch::get_rho() const { return rho; }

void Sketch::prefetch_offset_inc(uint32_t offset) const noexcept
{
  if (offset != OFF_INVALID) {
    sfhm->prefetch_inc(offset);
  }
}

void Sketch::prefetch_offset_enc(uint32_t offset) const noexcept
{
  // inc_v must already be in cache for this to be effective
  if (offset != OFF_INVALID) {
    sfhm->prefetch_enc(offset);
  }
}

bool Sketch::scan_bucket(uint32_t offset, enc_t enc_lr, uint32_t& hdist_min) const noexcept
{
  if (offset == OFF_INVALID) return false;
  const enc_t* ix1 = sfhm->bucket_ptr_start(offset);
  const enc_t* ix2 = sfhm->bucket_ptr_next(offset);
  uint32_t hmin = std::numeric_limits<uint32_t>::max();
  for (; ix1 < ix2; ++ix1) {
    const uint32_t hd = popcount_lr32((*ix1) ^ enc_lr);
    hmin = hd < hmin ? hd : hmin;
  }
  hdist_min = hmin;
  return true;
}

void Sketch::canonicalize()
{
  if (canonical) {
    return;
  }
  const uint64_t mask_bp = std::numeric_limits<uint64_t>::max() >> ((32 - k) * 2);
  sdhm_sptr_t sdhm = std::make_shared<SDHM>();
  sdhm->enc_vvec.resize(nrows);
  for (uint32_t off = 0; off < nrows; ++off) {
    const enc_t* ix1 = sfhm->bucket_ptr_start(off);
    const enc_t* ix2 = sfhm->bucket_ptr_next(off);
    for (; ix1 < ix2; ++ix1) {
      const uint64_t bp_ppos = lshf->inv_ppos_bp(off);
      const uint64_t bp_npos = lr64_to_bp64(lshf->inv_ppos_lr(*ix1));
      const uint64_t fw_bp = (bp_ppos | bp_npos) & mask_bp;
      const uint64_t rc_bp = revcomp_bp64(fw_bp, k);
      const uint64_t can_bp = std::max(fw_bp, rc_bp);
      const uint32_t rixn = lshf->compute_hash(can_bp);
      if (rixn >= nrows) continue;
      const enc_t new_enc = lshf->drop_ppos_lr(bp64_to_lr64(can_bp));
      sdhm->enc_vvec[rixn].push_back(new_enc);
    }
  }
  sdhm->sort_columns();
  sdhm->make_unique();
  sfhm = std::make_shared<SFHM>(sdhm);
  canonical = true;
}
