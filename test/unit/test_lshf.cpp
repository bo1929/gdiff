#include "doctest/doctest.h"
#include "lshf.hpp"
#include <set>
#include <random>

TEST_SUITE("LSHF") {

TEST_CASE("constructor sets k, h, m correctly") {
  LSHF lshf(27, 11, 2);
  CHECK(lshf.get_k() == 27);
  CHECK(lshf.get_h() == 11);
  CHECK(lshf.get_m() == 2);
}

TEST_CASE("ppos and npos partition {0, ..., k-1}") {
  LSHF lshf(27, 11, 2);
  auto ppos = lshf.get_ppos();
  auto npos = lshf.get_npos();

  CHECK(ppos.size() == 11);
  CHECK(npos.size() == 16);  // k - h = 27 - 11 = 16

  // All positions should be in range [0, k)
  std::set<uint8_t> all_pos;
  for (auto p : ppos) {
    CHECK(p < 27);
    all_pos.insert(p);
  }
  for (auto p : npos) {
    CHECK(p < 27);
    all_pos.insert(p);
  }
  // Union should be {0, 1, ..., 26}
  CHECK(all_pos.size() == 27);
}

TEST_CASE("compute_hash is deterministic for same input") {
  // Construct with known positions
  vec<uint8_t> ppos = {26, 25, 20, 18, 15, 14, 13, 10, 7, 5, 0};
  vec<uint8_t> npos;
  std::set<uint8_t> pset(ppos.begin(), ppos.end());
  for (uint8_t i = 0; i < 27; ++i) {
    if (pset.find(i) == pset.end()) npos.push_back(i);
  }

  LSHF lshf(2, ppos, npos);

  uint64_t enc_bp = 0x123456789ABCULL;
  uint32_t h1 = lshf.compute_hash(enc_bp);
  uint32_t h2 = lshf.compute_hash(enc_bp);
  CHECK(h1 == h2);
}

TEST_CASE("hash distribution is not degenerate") {
  LSHF lshf(27, 11, 2);
  std::mt19937 rng(42);
  std::uniform_int_distribution<uint64_t> dist(0, (1ULL << 54) - 1);

  std::set<uint32_t> hashes;
  for (int i = 0; i < 1000; ++i) {
    hashes.insert(lshf.compute_hash(dist(rng)));
  }
  // At least 100 distinct hash values out of 1000 random k-mers
  CHECK(hashes.size() > 100);
}

TEST_CASE("constructor from known positions") {
  vec<uint8_t> ppos = {5, 3, 1};
  vec<uint8_t> npos = {0, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                       19, 20, 21, 22, 23, 24, 25, 26};
  LSHF lshf(4, ppos, npos);
  CHECK(lshf.get_k() == 27);
  CHECK(lshf.get_h() == 3);
  CHECK(lshf.get_m() == 4);
}

TEST_CASE("inv_compute_hash + inv_drop_ppos_lr reconstruct bp64") {
  vec<uint8_t> ppos = {26, 25, 20, 18, 15, 14, 13, 10, 7, 5, 0};
  vec<uint8_t> npos;
  std::set<uint8_t> pset(ppos.begin(), ppos.end());
  for (uint8_t i = 0; i < 27; ++i) {
    if (pset.find(i) == pset.end()) npos.push_back(i);
  }
  LSHF lshf(2, ppos, npos);
  const uint64_t mask_bp = (1ULL << (2 * 27)) - 1;

  std::mt19937 rng(7);
  std::uniform_int_distribution<uint64_t> dist(0, mask_bp);
  for (int i = 0; i < 200; ++i) {
    const uint64_t bp = dist(rng) & mask_bp;
    const uint32_t rix = lshf.compute_hash(bp);
    const uint32_t enc_lr = lshf.drop_ppos_lr(bp64_to_lr64(bp));
    const uint64_t recon = (lshf.inv_compute_hash(rix) | lr64_to_bp64(lshf.inv_drop_ppos_lr(enc_lr))) & mask_bp;
    CHECK(recon == bp);
  }
}

} // TEST_SUITE
