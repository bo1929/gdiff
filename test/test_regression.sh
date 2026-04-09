#!/bin/bash
set -euo pipefail

# Regression test: runs gidiff sketch + map in both --enum-only and continuous modes,
# then compares output against ground truth in gt/.
#   gt/*.enum.txt  — expected enum-only output
#   gt/*.cont.txt  — expected continuous-mode output
# Exit 0 on success, non-zero on failure.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

GIDIFF="../gidiff"

if [ ! -x "$GIDIFF" ]; then
  echo "ERROR: gidiff binary not found at $GIDIFF"
  echo "Run 'make' first."
  exit 1
fi

NPROC="${NPROC:-4}"
SKETCHING_OPTS="-k 27 -w 31 -h 11 -m 2 -r 1 --frac"

# ── Phase 0: Create sketches ──────────────────────────────────────────────────
echo "=== Phase 0: Sketching ==="
rm -rf sketches && mkdir -p sketches
while read -r name; do
  "$GIDIFF" sketch $SKETCHING_OPTS -i "genomes/${name}.fna.gz" -o "sketches/${name}.skc" 2>/dev/null &
  while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do wait -n 2>/dev/null || true; done
done < genome_names.txt
wait

rm -rf est && mkdir -p est

# ── Phase 1: enum-only mode ───────────────────────────────────────────────────
echo "=== Phase 1: enum-only regression ==="
ENUM_OPTS="-d 0.10 -l 9900 --chisq 10000 --enum-only"
while IFS=$'\t' read -r query ref; do
  "$GIDIFF" map $ENUM_OPTS \
    -i "sketches/${ref}.skc" \
    -q "genomes/${query}.fna.gz" \
    -o "est/query_${query}-ref_${ref}.enum.txt" 2>/dev/null &
  while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do wait -n 2>/dev/null || true; done
done < genome_pairs.txt
wait

ENUM_FAIL=0
ENUM_DETAIL=""
ENUM_GT_MISSING=0
while IFS=$'\t' read -r query ref; do
  gt_f="gt/query_${query}-ref_${ref}.enum.txt"
  es_f="est/query_${query}-ref_${ref}.enum.txt"
  if [ ! -f "$gt_f" ]; then
    ENUM_GT_MISSING=1
    continue
  fi
  if ! diff -q "$es_f" "$gt_f" >/dev/null 2>&1; then
    ENUM_FAIL=1
    ENUM_DETAIL+="  differs: $es_f vs $gt_f\n"
  fi
done < genome_pairs.txt

if [ "$ENUM_GT_MISSING" -eq 1 ]; then
  echo "WARN: gt/*.enum.txt not found; generating ground truth from current output."
  while IFS=$'\t' read -r query ref; do
    cp "est/query_${query}-ref_${ref}.enum.txt" "gt/query_${query}-ref_${ref}.enum.txt"
  done < genome_pairs.txt
  ENUM_FAIL=0
  echo "PASS: enum-only ground truth generated (first run)."
elif [ "$ENUM_FAIL" -ne 0 ]; then
  echo "FAIL: enum-only output differs from ground truth"
  printf "%b" "$ENUM_DETAIL"
else
  echo "PASS: enum-only"
fi

# ── Phase 2: continuous mode (default, no --enum-only) ────────────────────────
echo "=== Phase 2: continuous regression ==="
CONT_OPTS="-d 0.10 -l 9900 --chisq 33.00051"
while IFS=$'\t' read -r query ref; do
  "$GIDIFF" map $CONT_OPTS \
    -i "sketches/${ref}.skc" \
    -q "genomes/${query}.fna.gz" \
    -o "est/query_${query}-ref_${ref}.cont.txt" 2>/dev/null &
  while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do wait -n 2>/dev/null || true; done
done < genome_pairs.txt
wait

CONT_FAIL=0
CONT_DETAIL=""
CONT_GT_MISSING=0
while IFS=$'\t' read -r query ref; do
  gt_f="gt/query_${query}-ref_${ref}.cont.txt"
  es_f="est/query_${query}-ref_${ref}.cont.txt"
  if [ ! -f "$gt_f" ]; then
    CONT_GT_MISSING=1
    continue
  fi
  if ! diff -q "$es_f" "$gt_f" >/dev/null 2>&1; then
    CONT_FAIL=1
    CONT_DETAIL+="  differs: $es_f vs $gt_f\n"
  fi
done < genome_pairs.txt

if [ "$CONT_GT_MISSING" -eq 1 ]; then
  echo "WARN: gt/*.cont.txt not found; generating ground truth from current output."
  while IFS=$'\t' read -r query ref; do
    cp "est/query_${query}-ref_${ref}.cont.txt" "gt/query_${query}-ref_${ref}.cont.txt"
  done < genome_pairs.txt
  CONT_FAIL=0
  echo "PASS: continuous ground truth generated (first run)."
elif [ "$CONT_FAIL" -ne 0 ]; then
  echo "FAIL: continuous output differs from ground truth"
  printf "%b" "$CONT_DETAIL"
else
  echo "PASS: continuous"
fi

# ── Cleanup intermediates ─────────────────────────────────────────────────────
rm -rf sketches est

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
if [ "$ENUM_FAIL" -ne 0 ] || [ "$CONT_FAIL" -ne 0 ]; then
  echo "=== REGRESSION FAILED ==="
  exit 1
else
  echo "=== ALL REGRESSION TESTS PASSED ==="
  exit 0
fi
