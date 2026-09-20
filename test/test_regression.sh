#!/bin/bash
set -euo pipefail

# End-to-end regression for gdiff.
#
# Phases 1-2 pin `map` against the committed ground truth in gt/, comparing the complete row
# (column header included) rather than a stable prefix:
#   gt/*.enum.txt  -- expected enum-only output
#   gt/*.cont.txt  -- expected continuous-mode output
# A legitimate change to any column, including the significance ones, means regenerating gt/:
# delete the files and re-run, and the bootstrap below rewrites them.
#
# The later phases are self-validating: they assert invariants (option equivalence, thread
# invariance, round trips, refused inputs) rather than comparing against new golden files, so
# they widen coverage without adding ground truth to maintain.
#
# Everything is written under one scratch directory that is removed on exit, so a run never
# leaves files in the repository and an interrupted run cleans up after itself.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

GDIFF="../gdiff"
NPROC="${NPROC:-4}"
SKETCHING_OPTS="-k 27 -w 31 -h 11"

if [ ! -x "$GDIFF" ]; then
  echo "ERROR: gdiff binary not found at $GDIFF"
  echo "Run 'make' first."
  exit 1
fi

WORK="$(mktemp -d "${TMPDIR:-/tmp}/gdiff-regression.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

FAIL=0
fail() {
  echo "FAIL: $*"
  FAIL=1
}

wait_slot() {
  while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do wait -n 2>/dev/null || true; done
}

# Every output is headed by two `#` provenance lines: the invocation, then the version with a
# timestamp. Both differ between runs by construction, so comparisons must drop them.
strip_comments() { grep -v '^#' "$1" || true; }

# True when an output holds at least one data line, i.e. more than its provenance header.
has_data() { [ -n "$(strip_comments "$1")" ]; }

first_name="$(head -n 1 genome_names.txt)"
first_query="$(cut -f1 genome_pairs.txt | head -n 1)"
first_ref="$(cut -f2 genome_pairs.txt | head -n 1)"
# Two distinct genomes, for the phases that need a pair of sketches.
names=($(head -n 2 genome_names.txt))

# -- Phase 0: sketch every genome into the scratch directory -------------------
echo "=== Phase 0: Sketching ==="
mkdir -p "$WORK/sketches"
while read -r name; do
  "$GDIFF" sketch $SKETCHING_OPTS -i "genomes/${name}.fna.gz" -o "$WORK/sketches/${name}.gs" 2>/dev/null &
  wait_slot
done < genome_names.txt
wait

# -- Phase 1: enum-only mode ---------------------------------------------------
echo "=== Phase 1: enum-only regression ==="
mkdir -p "$WORK/est"
ENUM_OPTS="-d 0.10 -l 9900 --hdist-th 4 --chisq 10000 --enum-only"
while IFS=$'\t' read -r query ref; do
  "$GDIFF" map $ENUM_OPTS \
    "genomes/${query}.fna.gz" \
    "$WORK/sketches/${ref}.gs" \
    -o "$WORK/est/query_${query}-ref_${ref}.enum.txt" 2>/dev/null &
  wait_slot
done < genome_pairs.txt
wait

ENUM_FAIL=0
ENUM_DETAIL=""
ENUM_GT_MISSING=0
while IFS=$'\t' read -r query ref; do
  gt_f="gt/query_${query}-ref_${ref}.enum.txt"
  es_f="$WORK/est/query_${query}-ref_${ref}.enum.txt"
  if [ ! -f "$gt_f" ]; then
    ENUM_GT_MISSING=1
    continue
  fi
  strip_comments "$es_f" > "$WORK/est_enum"
  if ! diff -q "$WORK/est_enum" "$gt_f" >/dev/null 2>&1; then
    ENUM_FAIL=1
    ENUM_DETAIL+="  differs: $es_f vs $gt_f\n"
  fi
done < genome_pairs.txt

if [ "$ENUM_GT_MISSING" -eq 1 ]; then
  echo "WARN: gt/*.enum.txt not found; generating ground truth from current output."
  while IFS=$'\t' read -r query ref; do
    strip_comments "$WORK/est/query_${query}-ref_${ref}.enum.txt" > "gt/query_${query}-ref_${ref}.enum.txt"
  done < genome_pairs.txt
  ENUM_FAIL=0
  echo "PASS: enum-only ground truth generated (first run)."
elif [ "$ENUM_FAIL" -ne 0 ]; then
  echo "FAIL: enum-only output differs from ground truth"
  printf "%b" "$ENUM_DETAIL"
else
  echo "PASS: enum-only"
fi

# -- Phase 2: continuous mode (default, no --enum-only) ------------------------
echo "=== Phase 2: continuous regression ==="
CONT_OPTS="-d 0.10 -l 9900 --hdist-th 4 --chisq 33.00051"
while IFS=$'\t' read -r query ref; do
  "$GDIFF" map $CONT_OPTS \
    "genomes/${query}.fna.gz" \
    "$WORK/sketches/${ref}.gs" \
    -o "$WORK/est/query_${query}-ref_${ref}.cont.txt" 2>/dev/null &
  wait_slot
done < genome_pairs.txt
wait

CONT_FAIL=0
CONT_DETAIL=""
CONT_GT_MISSING=0
while IFS=$'\t' read -r query ref; do
  gt_f="gt/query_${query}-ref_${ref}.cont.txt"
  es_f="$WORK/est/query_${query}-ref_${ref}.cont.txt"
  if [ ! -f "$gt_f" ]; then
    CONT_GT_MISSING=1
    continue
  fi
  strip_comments "$es_f" > "$WORK/est_cont"
  if ! diff -q "$WORK/est_cont" "$gt_f" >/dev/null 2>&1; then
    CONT_FAIL=1
    CONT_DETAIL+="  differs: $es_f vs $gt_f\n"
  fi
done < genome_pairs.txt

if [ "$CONT_GT_MISSING" -eq 1 ]; then
  echo "WARN: gt/*.cont.txt not found; generating ground truth from current output."
  while IFS=$'\t' read -r query ref; do
    strip_comments "$WORK/est/query_${query}-ref_${ref}.cont.txt" > "gt/query_${query}-ref_${ref}.cont.txt"
  done < genome_pairs.txt
  CONT_FAIL=0
  echo "PASS: continuous ground truth generated (first run)."
elif [ "$CONT_FAIL" -ne 0 ]; then
  echo "FAIL: continuous output differs from ground truth"
  printf "%b" "$CONT_DETAIL"
else
  echo "PASS: continuous"
fi

# -- Phase 3: refused inputs ---------------------------------------------------
echo "=== Phase 3: refused inputs ==="
if "$GDIFF" map --hdist-th 8 -d 0.10 -l 9900 \
  "genomes/${first_name}.fna.gz" "$WORK/sketches/${first_name}.gs" >/dev/null 2>&1; then
  fail "map accepted unsupported --hdist-th 8"
fi
# A distance-threshold width of neither 1 nor 8 lanes has no meaning.
if "$GDIFF" map -d 0.1 -d 0.2 -d 0.3 -l 9900 \
  "genomes/${first_query}.fna.gz" "$WORK/sketches/${first_ref}.gs" >/dev/null 2>&1; then
  fail "map accepted three distance lanes"
fi
# The threshold count is checked before any sequence is read.
if "$GDIFF" map -d 0.1 -d 0.2 -l 9900 \
  "genomes/${first_query}.fna.gz" "$WORK/sketches/${first_ref}.gs" >/dev/null 2>&1; then
  fail "map accepted two distance lanes"
fi
if "$GDIFF" map -d 0.1 -l 9900 \
  "genomes/${first_query}.fna.gz" "$WORK/does_not_exist.gs" >/dev/null 2>&1; then
  fail "map accepted a missing reference container"
fi
echo "PASS: refused inputs"

# -- Phase 4: map option equivalence ------------------------------------------
echo "=== Phase 4: map option equivalence ==="
MAP_BASE="-l 9900 --hdist-th 4 --sample-size 100"
map_once() { # $1 = output file, rest = extra options
  local out="$1"
  shift
  "$GDIFF" map "$@" $MAP_BASE "genomes/${first_query}.fna.gz" "$WORK/sketches/${first_ref}.gs" \
    > "$out" 2>/dev/null
}

map_once "$WORK/map_ref.tsv" -d 0.10
strip_comments "$WORK/map_ref.tsv" > "$WORK/map_ref.data"
has_data "$WORK/map_ref.tsv" || fail "map produced no records at all"

# -o must write exactly what stdout would.
"$GDIFF" map -o "$WORK/map_o.tsv" -d 0.10 $MAP_BASE \
  "genomes/${first_query}.fna.gz" "$WORK/sketches/${first_ref}.gs" >/dev/null 2>&1
strip_comments "$WORK/map_o.tsv" > "$WORK/map_o.data"
cmp -s "$WORK/map_ref.data" "$WORK/map_o.data" || fail "map -o differs from the stdout report"

# The worker count must not change the records.
"$GDIFF" --num-threads 4 map -d 0.10 $MAP_BASE \
  "genomes/${first_query}.fna.gz" "$WORK/sketches/${first_ref}.gs" > "$WORK/map_t4.tsv" 2>/dev/null
strip_comments "$WORK/map_t4.tsv" > "$WORK/map_t4.data"
cmp -s "$WORK/map_ref.data" "$WORK/map_t4.data" || fail "map output depends on --num-threads"

# --levels resolves its own thresholds from the background and must stay well formed.
map_once "$WORK/map_levels.tsv" --levels 0.1 0.05 0.01 0.005
has_data "$WORK/map_levels.tsv" || fail "map --levels produced no output"

# --per-sequence and -b are alternate modes, not failures.
map_once "$WORK/map_perseq.tsv" -d 0.10 --per-sequence
has_data "$WORK/map_perseq.tsv" || fail "map --per-sequence produced no output"
map_once "$WORK/map_b2.tsv" -d 0.10 -b 2
has_data "$WORK/map_b2.tsv" || fail "map -b 2 produced no output"
echo "PASS: map option equivalence"

# -- Phase 5: dist -------------------------------------------------------------
echo "=== Phase 5: dist ==="
sketch_a="$WORK/sketches/${names[0]}.gs"
sketch_b="$WORK/sketches/${names[1]}.gs"

# A self-pair exercises the symmetric reconciliation on real data.
"$GDIFF" dist "$sketch_a" "$sketch_a" > "$WORK/dist_self.tsv" 2>/dev/null
has_data "$WORK/dist_self.tsv" || fail "dist produced no output"

# Two lists select the same cross product as two positionals.
"$GDIFF" dist "$sketch_a" "$sketch_b" > "$WORK/dist_pos.tsv" 2>/dev/null
printf '%s\n' "$sketch_a" > "$WORK/list_a.txt"
printf '%s\n' "$sketch_b" > "$WORK/list_b.txt"
"$GDIFF" dist --list-a "$WORK/list_a.txt" --list-b "$WORK/list_b.txt" > "$WORK/dist_list.tsv" 2>/dev/null
strip_comments "$WORK/dist_pos.tsv" > "$WORK/dist_pos.data"
strip_comments "$WORK/dist_list.tsv" > "$WORK/dist_list.data"
cmp -s "$WORK/dist_pos.data" "$WORK/dist_list.data" ||
  fail "dist --list-a/--list-b differ from the positional form"

# -o must write exactly what stdout would.
"$GDIFF" dist -o "$WORK/dist_o.tsv" "$sketch_a" "$sketch_b" >/dev/null 2>&1
strip_comments "$WORK/dist_o.tsv" > "$WORK/dist_o.data"
cmp -s "$WORK/dist_pos.data" "$WORK/dist_o.data" || fail "dist -o differs from the stdout report"

# --output-samples switches to one line per sampled window.
"$GDIFF" dist --output-samples "$sketch_a" "$sketch_b" > "$WORK/dist_samples.tsv" 2>/dev/null
has_data "$WORK/dist_samples.tsv" || fail "dist --output-samples produced no output"
[ "$(strip_comments "$WORK/dist_samples.tsv" | head -n 1 | cut -f1)" = "config" ] ||
  fail "dist --output-samples header changed"
[ "$(strip_comments "$WORK/dist_samples.tsv" | head -n 1 | cut -f8)" = "direction" ] ||
  fail "dist --output-samples lost the direction column"
echo "PASS: dist"

# -- Phase 6: roll -------------------------------------------------------------
echo "=== Phase 6: roll ==="
"$GDIFF" roll -l 100 -s 500 "genomes/${first_query}.fna.gz" "$sketch_a" > "$WORK/roll.tsv" 2>/dev/null
has_data "$WORK/roll.tsv" || fail "roll produced no output"
[ "$(strip_comments "$WORK/roll.tsv" | head -n 1 | cut -f1)" = "seq" ] || fail "roll header changed"

"$GDIFF" roll -l 100 -s 500 -o "$WORK/roll_o.tsv" \
  "genomes/${first_query}.fna.gz" "$sketch_a" >/dev/null 2>&1
strip_comments "$WORK/roll.tsv" > "$WORK/roll.data"
strip_comments "$WORK/roll_o.tsv" > "$WORK/roll_o.data"
cmp -s "$WORK/roll.data" "$WORK/roll_o.data" || fail "roll -o differs from the stdout report"
echo "PASS: roll"

# -- Phase 7: merge and info ---------------------------------------------------
echo "=== Phase 7: merge and info ==="
"$GDIFF" merge -i "$WORK/sketches/${names[0]}.gs" -i "$WORK/sketches/${names[1]}.gs" \
  -o "$WORK/merged.gs" >/dev/null 2>&1
"$GDIFF" info -i "$WORK/merged.gs" > "$WORK/info.txt" 2>/dev/null

grep -qE '^Sketches:[[:space:]]+2$' "$WORK/info.txt" || fail "merged container does not hold two sketches"
blocks="$(grep -c '^\[Sketch ' "$WORK/info.txt" || true)"
[ "$blocks" = "2" ] || fail "info reported $blocks sketch blocks, expected 2"

# A merged container is an ordinary input, and its within-set comparison is the same
# whether it is given positionally or through --list-a.
"$GDIFF" dist "$WORK/merged.gs" > "$WORK/dist_merged.tsv" 2>/dev/null
has_data "$WORK/dist_merged.tsv" || fail "dist rejected the merged container"
printf '%s\n' "$WORK/merged.gs" > "$WORK/list_merged.txt"
"$GDIFF" dist --list-a "$WORK/list_merged.txt" > "$WORK/dist_merged_list.tsv" 2>/dev/null
strip_comments "$WORK/dist_merged.tsv" > "$WORK/dist_merged.data"
strip_comments "$WORK/dist_merged_list.tsv" > "$WORK/dist_merged_list.data"
cmp -s "$WORK/dist_merged.data" "$WORK/dist_merged_list.data" ||
  fail "dist --list-a differs from the positional form for a within-set pair"

# Containers with different configurations cannot be merged, and a refusal writes nothing.
"$GDIFF" sketch -k 23 -i "genomes/${names[0]}.fna.gz" -o "$WORK/k23.gs" >/dev/null 2>&1
if "$GDIFF" merge -i "$WORK/merged.gs" -i "$WORK/k23.gs" -o "$WORK/refused.gs" >/dev/null 2>&1; then
  fail "merge accepted containers with different configurations"
fi
[ -e "$WORK/refused.gs" ] && fail "a refused merge left an output file behind"
echo "PASS: merge and info"

# -- Phase 8: provenance header ------------------------------------------------
echo "=== Phase 8: provenance header ==="
# Text outputs carry the two lines in the data stream; sketch and merge write a binary container,
# so theirs go to stderr.
for probe in map dist roll info; do
  case "$probe" in
    map)  "$GDIFF" map -d 0.10 $MAP_BASE "genomes/${first_query}.fna.gz" "$WORK/sketches/${first_ref}.gs" \
            > "$WORK/prov_$probe.txt" 2>/dev/null ;;
    dist) "$GDIFF" dist "$sketch_a" "$sketch_b" > "$WORK/prov_$probe.txt" 2>/dev/null ;;
    roll) "$GDIFF" roll -l 100 -s 500 "genomes/${first_query}.fna.gz" "$sketch_a" \
            > "$WORK/prov_$probe.txt" 2>/dev/null ;;
    info) "$GDIFF" info -i "$sketch_a" > "$WORK/prov_$probe.txt" 2>/dev/null ;;
  esac
  [ "$(head -n 1 "$WORK/prov_$probe.txt" | cut -c1-14)" = "# invocation: " ] ||
    fail "$probe output is missing its invocation line"
  [ "$(head -n 2 "$WORK/prov_$probe.txt" | tail -n 1 | cut -c1-11)" = "# version: " ] ||
    fail "$probe output is missing its version line"
done

for probe in sketch merge; do
  if [ "$probe" = "sketch" ]; then
    "$GDIFF" sketch -k 27 -w 31 -h 11 -i "genomes/${first_query}.fna.gz" \
      -o "$WORK/prov_sketch.gs" >/dev/null 2> "$WORK/prov_$probe.txt"
  else
    "$GDIFF" merge -i "$sketch_a" -i "$sketch_b" -o "$WORK/prov_merge.gs" >/dev/null 2> "$WORK/prov_$probe.txt"
  fi
  grep -q '^# invocation: ' "$WORK/prov_$probe.txt" || fail "$probe stderr is missing its invocation line"
  grep -q '^# version: gdiff ' "$WORK/prov_$probe.txt" || fail "$probe stderr is missing its version line"
done
echo "PASS: provenance header"

# -- Phase 9: map output schema ------------------------------------------------
echo "=== Phase 9: map output schema ==="
# `map` names its columns so `plot.py` can find them by name instead of by position; the
# strand-aware variant is not covered by the gt/ files above.
EXPECT_AGNOSTIC="seq seq_len start end reference d mask d_bin d_q d_acc percentile fold qvalue info lr_ub"
EXPECT_AWARE="seq seq_len start end strand is_rc reference d mask d_bin d_q d_diff d_acc percentile fold qvalue info lr_ub"

agnostic_header="$(strip_comments "$WORK/map_ref.tsv" | head -n 1 | tr '\t' ' ')"
[ "$agnostic_header" = "$EXPECT_AGNOSTIC" ] || fail "strand-agnostic map header changed: $agnostic_header"

"$GDIFF" sketch -k 27 -w 31 -h 11 --strand-aware -i "genomes/${first_query}.fna.gz" \
  -o "$WORK/aware.gs" >/dev/null 2>&1
"$GDIFF" map -d 0.10 -l 9900 --hdist-th 4 --sample-size 100 "genomes/${first_query}.fna.gz" \
  "$WORK/aware.gs" > "$WORK/aware.tsv" 2>/dev/null
aware_header="$(strip_comments "$WORK/aware.tsv" | head -n 1 | tr '\t' ' ')"
[ "$aware_header" = "$EXPECT_AWARE" ] || fail "strand-aware map header changed: $aware_header"

# Every data row carries as many fields as its header.
for probe in map_ref aware; do
  fields="$(strip_comments "$WORK/$probe.tsv" | head -n 1 | awk -F'\t' '{print NF}')"
  bad="$(strip_comments "$WORK/$probe.tsv" | awk -F'\t' -v n="$fields" 'NF != n' | wc -l | tr -d ' ')"
  [ "$bad" = "0" ] || fail "$probe has $bad row(s) whose field count differs from the header"
done
echo "PASS: map output schema"

# -- Summary -------------------------------------------------------------------
echo ""
if [ "$FAIL" -ne 0 ] || [ "$ENUM_FAIL" -ne 0 ] || [ "$CONT_FAIL" -ne 0 ]; then
  echo "=== REGRESSION FAILED ==="
  exit 1
fi
echo "=== ALL REGRESSION TESTS PASSED ==="
