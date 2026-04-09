#include "doctest/doctest.h"
#include "enc.hpp"
#include <cstring>

TEST_SUITE("SEQ_NT4_TABLE") {

TEST_CASE("standard bases map to 0-3") {
  CHECK(SEQ_NT4_TABLE['A'] == 0);
  CHECK(SEQ_NT4_TABLE['C'] == 1);
  CHECK(SEQ_NT4_TABLE['G'] == 2);
  CHECK(SEQ_NT4_TABLE['T'] == 3);
}

TEST_CASE("lowercase bases map to 0-3") {
  CHECK(SEQ_NT4_TABLE['a'] == 0);
  CHECK(SEQ_NT4_TABLE['c'] == 1);
  CHECK(SEQ_NT4_TABLE['g'] == 2);
  CHECK(SEQ_NT4_TABLE['t'] == 3);
}

TEST_CASE("N and other characters map to 4") {
  CHECK(SEQ_NT4_TABLE['N'] == 4);
  CHECK(SEQ_NT4_TABLE['n'] == 4);
  CHECK(SEQ_NT4_TABLE['X'] == 4);
  CHECK(SEQ_NT4_TABLE[' '] == 4);
  CHECK(SEQ_NT4_TABLE[0]   == 4);
}

} // TEST_SUITE

TEST_SUITE("compute_encoding") {

TEST_CASE("short k-mer ACGT") {
  const char* seq = "ACGT";
  uint64_t enc_lr, enc_bp;
  compute_encoding(seq, seq + 4, enc_lr, enc_bp);
  // A=00, C=01, G=10, T=11 in bp encoding => 0b00011011 = 0x1B = 27
  CHECK(enc_bp == 0x1B);
  // enc_lr should also be non-zero
  CHECK(enc_lr != 0);
}

TEST_CASE("single base A") {
  const char* seq = "A";
  uint64_t enc_lr, enc_bp;
  compute_encoding(seq, seq + 1, enc_lr, enc_bp);
  CHECK(enc_bp == 0);
  CHECK(enc_lr == 0);
}

TEST_CASE("single base T") {
  const char* seq = "T";
  uint64_t enc_lr, enc_bp;
  compute_encoding(seq, seq + 1, enc_lr, enc_bp);
  CHECK(enc_bp == 3); // T = 11
}

TEST_CASE("all-A sequence") {
  const char* seq = "AAAAAAAAAA"; // 10 A's
  uint64_t enc_lr, enc_bp;
  compute_encoding(seq, seq + 10, enc_lr, enc_bp);
  CHECK(enc_bp == 0);
  CHECK(enc_lr == 0);
}

TEST_CASE("all-T sequence") {
  const char* seq = "TTTTTTTTTT"; // 10 T's
  uint64_t enc_lr, enc_bp;
  compute_encoding(seq, seq + 10, enc_lr, enc_bp);
  // T=11, 10 bases -> 20 bits all 1: 0xFFFFF
  CHECK(enc_bp == 0xFFFFF);
}

} // TEST_SUITE

TEST_SUITE("update_encoding") {

TEST_CASE("update matches fresh computation") {
  // Compute encoding for ACGTG from scratch
  const char* seq = "ACGTG";
  uint64_t enc_lr_full, enc_bp_full;
  compute_encoding(seq, seq + 5, enc_lr_full, enc_bp_full);

  // Compute ACGT first, then update with G
  uint64_t enc_lr, enc_bp;
  compute_encoding(seq, seq + 4, enc_lr, enc_bp);
  update_encoding(seq + 4, enc_lr, enc_bp);

  // After masking to 5 bases (10 bits for bp):
  uint64_t mask_bp = (1ULL << 10) - 1;
  CHECK((enc_bp & mask_bp) == (enc_bp_full & mask_bp));
}

TEST_CASE("sliding window consistency for k=4") {
  const char* seq = "ACGTACGT";
  const uint8_t k = 4;
  const uint64_t mask_bp = (1ULL << (2 * k)) - 1;

  // First window: compute from scratch
  uint64_t enc_lr, enc_bp;
  compute_encoding(seq, seq + k, enc_lr, enc_bp);
  enc_bp &= mask_bp;

  // Slide through remaining positions and verify against fresh computation
  for (size_t i = k; i < 8; ++i) {
    update_encoding(seq + i, enc_lr, enc_bp);
    enc_bp &= mask_bp;

    uint64_t ref_lr, ref_bp;
    compute_encoding(seq + i - k + 1, seq + i + 1, ref_lr, ref_bp);
    ref_bp &= mask_bp;
    CHECK(enc_bp == ref_bp);
  }
}

} // TEST_SUITE

TEST_SUITE("revcomp_bp64") {

TEST_CASE("ACGT reverse complement is ACGT") {
  // ACGT -> rc(ACGT) = ACGT
  // bp: A=00 C=01 G=10 T=11 -> 00011011 = 0x1B
  // rc: A=00 C=01 G=10 T=11 -> 0x1B
  uint64_t enc = 0x1B;
  uint64_t rc = revcomp_bp64(enc, 4);
  CHECK(rc == 0x1B);
}

TEST_CASE("self-complementary palindrome") {
  // AATT -> revcomp = AATT
  uint64_t enc_lr, enc_bp;
  const char* seq = "AATT";
  compute_encoding(seq, seq + 4, enc_lr, enc_bp);
  uint64_t mask = (1ULL << 8) - 1;
  enc_bp &= mask;
  uint64_t rc = revcomp_bp64(enc_bp, 4);
  CHECK(rc == enc_bp);
}

TEST_CASE("double revcomp is identity") {
  const char* seq = "ACGTACGTACGTACGTACGTACGTACG"; // 27-mer
  uint64_t enc_lr, enc_bp;
  compute_encoding(seq, seq + 27, enc_lr, enc_bp);
  uint64_t mask = (1ULL << 54) - 1;
  enc_bp &= mask;

  uint64_t rc1 = revcomp_bp64(enc_bp, 27);
  uint64_t rc2 = revcomp_bp64(rc1, 27);
  CHECK(rc2 == enc_bp);
}

TEST_CASE("A -> T, single base") {
  // A = 00, revcomp should be T = 11
  CHECK(revcomp_bp64(0, 1) == 3);
  // T = 11, revcomp should be A = 00
  CHECK(revcomp_bp64(3, 1) == 0);
  // C = 01, revcomp should be G = 10
  CHECK(revcomp_bp64(1, 1) == 2);
  // G = 10, revcomp should be C = 01
  CHECK(revcomp_bp64(2, 1) == 1);
}

} // TEST_SUITE

TEST_SUITE("bp64_to_lr64") {

TEST_CASE("round-trip consistency with hdist_lr64") {
  // Two identical k-mers should have hdist 0
  const char* seq = "ACGTACGT";
  uint64_t enc_lr, enc_bp;
  compute_encoding(seq, seq + 8, enc_lr, enc_bp);

  // Convert bp -> lr and compare
  uint64_t lr_from_bp = bp64_to_lr64(enc_bp);
  // hdist of identical sequences should be 0
  CHECK(hdist_lr64(lr_from_bp, lr_from_bp) == 0);
}

TEST_CASE("single-base difference gives hdist 1") {
  // ACGT vs ACGA: differ at position 3 (T->A)
  const char* seq1 = "ACGT";
  const char* seq2 = "ACGA";
  uint64_t lr1, bp1, lr2, bp2;
  compute_encoding(seq1, seq1 + 4, lr1, bp1);
  compute_encoding(seq2, seq2 + 4, lr2, bp2);

  uint64_t lr_from_bp1 = bp64_to_lr64(bp1);
  uint64_t lr_from_bp2 = bp64_to_lr64(bp2);
  CHECK(hdist_lr64(lr_from_bp1, lr_from_bp2) == 1);
}

} // TEST_SUITE

TEST_SUITE("hdist_lr64") {

TEST_CASE("identical sequences have hdist 0") {
  uint64_t lr, bp;
  const char* seq = "ACGTACGTACGTACGTACGTACGTACG"; // 27-mer
  compute_encoding(seq, seq + 27, lr, bp);
  CHECK(hdist_lr64(lr, lr) == 0);
}

TEST_CASE("all-A vs all-T has max hdist for length") {
  // A=00 vs T=11: every position differs -> hdist = k
  const char* seqA = "AAAA";
  const char* seqT = "TTTT";
  uint64_t lrA, bpA, lrT, bpT;
  compute_encoding(seqA, seqA + 4, lrA, bpA);
  compute_encoding(seqT, seqT + 4, lrT, bpT);
  CHECK(hdist_lr64(lrA, lrT) == 4);
}

} // TEST_SUITE

TEST_SUITE("extract_bits") {

TEST_CASE("extracts correct bits") {
  // x = 0b11010110, mask = 0b10100100 -> extract bits at positions 2,5,7
  // bit 2 = 1 (-> result bit 0), bit 5 = 0 (-> result bit 1), bit 7 = 1 (-> result bit 2)
  // result = 0b101 = 5
  uint64_t x = 0b11010110;
  uint64_t mask = 0b10100100;
  CHECK(extract_bits(x, mask) == 0b101);
}

TEST_CASE("empty mask gives 0") {
  CHECK(extract_bits(uint64_t(0xFF), uint64_t(0)) == 0);
}

TEST_CASE("full mask gives identity for small values") {
  uint64_t x = 0b1010;
  uint64_t mask = 0b1111;
  CHECK(extract_bits(x, mask) == 0b1010);
}

} // TEST_SUITE
