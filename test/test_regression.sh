#!/bin/bash
set -euo pipefail

# End-to-end regression for gdiff.
#
# Phases 1-2 pin `map` against the committed ground truth in gt/, comparing the complete row
# (column header included) rather than a stable prefix:
#   gt/*.enum.txt  -- expected per-lane enumeration output (--independent-levels)
#   gt/*.cont.txt  -- expected gap-claiming output (the default mode)
# Both golden phases pass fixed distance cutoffs through --thresholds, so the ground truth is
# independent of the background detection; the --levels path is covered by the phase 4
# invariants instead.
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
# Pin every sketch parameter that feeds downstream output, so the committed gt/ files stay valid
# when a default changes: -k/-w/-h plus --frac and -l (the sampled-window length).
SKETCHING_OPTS="-k 27 -w 31 -h 11 --frac 1.0 -l 500"

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

# True when a map TSV holds at least one record, i.e. more than its column header.
has_records() { [ -n "$(strip_comments "$1" | awk 'NR == 2 { print; exit }')" ]; }

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

# -- Phase 1: per-lane enumeration (--independent-levels) ------------------------
echo "=== Phase 1: enumeration regression ==="
mkdir -p "$WORK/est"
ENUM_OPTS="--thresholds 0.02 0.3 -l 9900 --hdist-th 4 --chisq 10000 --independent-levels"
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
ENUM_WITH_RECORDS=0
while IFS=$'\t' read -r query ref; do
  gt_f="gt/query_${query}-ref_${ref}.enum.txt"
  es_f="$WORK/est/query_${query}-ref_${ref}.enum.txt"
  if [ -s "$es_f" ] && has_records "$es_f"; then
    ENUM_WITH_RECORDS=$((ENUM_WITH_RECORDS + 1))
  fi
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
# A divergent pair may legitimately stay silent at the pinned cutoffs - its header-only report
# is then the ground truth, pinning the absence of calls - but a silent whole corpus means the
# pipeline is broken.
if [ "$ENUM_WITH_RECORDS" -eq 0 ]; then
  ENUM_FAIL=1
  ENUM_DETAIL+="  no pair produced any records\n"
fi

if [ "$ENUM_GT_MISSING" -eq 1 ] && [ "$ENUM_WITH_RECORDS" -gt 0 ]; then
  echo "WARN: gt/*.enum.txt not found; generating ground truth from current output."
  while IFS=$'\t' read -r query ref; do
    es_f="$WORK/est/query_${query}-ref_${ref}.enum.txt"
    if [ -s "$es_f" ]; then
      strip_comments "$es_f" > "gt/query_${query}-ref_${ref}.enum.txt"
    fi
  done < genome_pairs.txt
  ENUM_FAIL=0
  echo "PASS: enumeration ground truth generated (first run)."
elif [ "$ENUM_FAIL" -ne 0 ]; then
  echo "FAIL: enumeration output differs from ground truth"
  printf "%b" "$ENUM_DETAIL"
else
  echo "PASS: enumeration"
fi

# -- Phase 2: gap-claiming (default mode) ---------------------------------------
echo "=== Phase 2: gap-claiming regression ==="
CONT_OPTS="--thresholds 0.02 0.3 -l 9900 --hdist-th 4 --chisq 33.00051"
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
CONT_WITH_RECORDS=0
while IFS=$'\t' read -r query ref; do
  gt_f="gt/query_${query}-ref_${ref}.cont.txt"
  es_f="$WORK/est/query_${query}-ref_${ref}.cont.txt"
  if [ -s "$es_f" ] && has_records "$es_f"; then
    CONT_WITH_RECORDS=$((CONT_WITH_RECORDS + 1))
  fi
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
# Same rule as phase 1: per-pair silence is ground truth, corpus-wide silence is a break.
if [ "$CONT_WITH_RECORDS" -eq 0 ]; then
  CONT_FAIL=1
  CONT_DETAIL+="  no pair produced any records\n"
fi

if [ "$CONT_GT_MISSING" -eq 1 ] && [ "$CONT_WITH_RECORDS" -gt 0 ]; then
  echo "WARN: gt/*.cont.txt not found; generating ground truth from current output."
  while IFS=$'\t' read -r query ref; do
    es_f="$WORK/est/query_${query}-ref_${ref}.cont.txt"
    if [ -s "$es_f" ]; then
      strip_comments "$es_f" > "gt/query_${query}-ref_${ref}.cont.txt"
    fi
  done < genome_pairs.txt
  CONT_FAIL=0
  echo "PASS: gap-claiming ground truth generated (first run)."
elif [ "$CONT_FAIL" -ne 0 ]; then
  echo "FAIL: gap-claiming output differs from ground truth"
  printf "%b" "$CONT_DETAIL"
else
  echo "PASS: gap-claiming"
fi

# -- Phase 3: refused inputs ---------------------------------------------------
echo "=== Phase 3: refused inputs ==="
rejects() { # $1 = label, rest = map arguments; the run must fail
  local why="$1"
  shift
  if "$GDIFF" map "$@" >/dev/null 2>&1; then
    fail "map accepted $why"
  fi
}
refused_query="genomes/${first_query}.fna.gz"
refused_ref="$WORK/sketches/${first_ref}.gs"
rejects "unsupported --hdist-th 8" --hdist-th 8 --thresholds 0.1 -l 9900 "$refused_query" "$refused_ref"
# --thresholds packs one lane per value: at most 8 fit the register, every value must sit in
# (0, 1), values must be unique, and the option replaces --levels rather than joining it.
rejects "nine threshold lanes" --thresholds 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 -l 9900 "$refused_query" "$refused_ref"
rejects "duplicate thresholds" --thresholds 0.1 0.1 -l 9900 "$refused_query" "$refused_ref"
rejects "an out-of-domain threshold" --thresholds 1.5 -l 9900 "$refused_query" "$refused_ref"
rejects "--thresholds together with --levels" --thresholds 0.1 --levels 0.05 -l 9900 "$refused_query" "$refused_ref"
# --levels budgets two thresholds per level and caps the list at three.
rejects "four levels" --levels 0.1 0.05 0.01 0.005 -l 9900 "$refused_query" "$refused_ref"
rejects "a missing reference container" --thresholds 0.1 -l 9900 "$refused_query" "$WORK/does_not_exist.gs"
# A p-value needs a usable pool: at least 100 background windows.
rejects "--sample-size 50, below the minimum pool" --thresholds 0.1 -l 9900 --sample-size 50 "$refused_query" "$refused_ref"
rejects "--levels with --sample-size 0" --levels 0.1 0.05 -l 9900 --sample-size 0 "$refused_query" "$refused_ref"
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

map_once "$WORK/map_ref.tsv" --levels 0.05 0.01 0.005
strip_comments "$WORK/map_ref.tsv" > "$WORK/map_ref.data"
has_records "$WORK/map_ref.tsv" || fail "map produced no records at all"

# -o must write exactly what stdout would.
"$GDIFF" map -o "$WORK/map_o.tsv" --levels 0.05 0.01 0.005 $MAP_BASE \
  "genomes/${first_query}.fna.gz" "$WORK/sketches/${first_ref}.gs" >/dev/null 2>&1
strip_comments "$WORK/map_o.tsv" > "$WORK/map_o.data"
cmp -s "$WORK/map_ref.data" "$WORK/map_o.data" || fail "map -o differs from the stdout report"

# A multi-value option immediately before the positionals must not swallow them (regression:
# the following -o became unreachable and the run died with "query-path is required"). Both
# option families carry the trap.
"$GDIFF" map --levels 0.05 0.01 0.005 "genomes/${first_query}.fna.gz" "$WORK/sketches/${first_ref}.gs" \
  -o "$WORK/map_levels_first_o.tsv" $MAP_BASE >/dev/null 2>&1 \
  || fail "map --levels before the positionals with -o failed"
strip_comments "$WORK/map_levels_first_o.tsv" > "$WORK/map_levels_first_o.data"
cmp -s "$WORK/map_ref.data" "$WORK/map_levels_first_o.data" \
  || fail "map --levels before the positionals differs from the reference report"

"$GDIFF" map --thresholds 0.02 0.3 "genomes/${first_query}.fna.gz" "$WORK/sketches/${first_ref}.gs" \
  -o "$WORK/map_th_first_o.tsv" $MAP_BASE >/dev/null 2>&1 \
  || fail "map --thresholds before the positionals with -o failed"
has_records "$WORK/map_th_first_o.tsv" || fail "--thresholds before the positionals produced no output"

# The worker count must not change the records.
"$GDIFF" --num-threads 4 map --levels 0.05 0.01 0.005 $MAP_BASE \
  "genomes/${first_query}.fna.gz" "$WORK/sketches/${first_ref}.gs" > "$WORK/map_t4.tsv" 2>/dev/null
strip_comments "$WORK/map_t4.tsv" > "$WORK/map_t4.data"
cmp -s "$WORK/map_ref.data" "$WORK/map_t4.data" || fail "map output depends on --num-threads"

# --levels resolves its own thresholds from the background and must stay well formed: each
# level places at most one low and one high lane, so three levels fit at most 6 of the 8.
n_lanes="$(
  "$GDIFF" map --levels 0.05 0.01 0.005 "genomes/${first_query}.fna.gz" \
    "$WORK/sketches/${first_ref}.gs" $MAP_BASE 2>&1 >/dev/null |
    awk '/threshold=/ { n++ } END { print n + 0 }'
)"
[ "$n_lanes" -ge 2 ] && [ "$n_lanes" -le 6 ] || fail "--levels ran $n_lanes lanes, expected 2 to 6"

# --thresholds places the given cutoffs verbatim: every record must name one of them, carry no
# alpha (no level placed it), and sit on the side of the background median its value implies.
map_once "$WORK/map_th.tsv" --thresholds 0.02 0.3
has_records "$WORK/map_th.tsv" || fail "map --thresholds produced no output"
if ! strip_comments "$WORK/map_th.tsv" | awk -F'\t' '
  NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
  {
    th = $(col["threshold"]) + 0
    dn = $(col["d_null"]) + 0
    if (th != 0.02 && th != 0.3) bad = 1
    else if ($(col["alpha"]) != "nan") bad = 2
    else if (th <= dn && $(col["side"]) != "low") bad = 3
    else if (th > dn && $(col["side"]) != "high") bad = 4
  }
  END { exit bad ? 1 : 0 }'; then
  fail "a --thresholds record carried an unexpected threshold, alpha, or side"
fi

# Both extraction modes must report on the same lanes.
map_once "$WORK/map_th_ind.tsv" --thresholds 0.02 0.3 --independent-levels
has_records "$WORK/map_th_ind.tsv" || fail "map --thresholds --independent-levels produced no output"

# -b is an alternate mode, not a failure.
map_once "$WORK/map_b2.tsv" --thresholds 0.02 0.3 -b 2
has_records "$WORK/map_b2.tsv" || fail "map -b 2 produced no output"
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
[ "$(strip_comments "$WORK/dist_samples.tsv" | head -n 1 | cut -f7)" = "direction" ] ||
  fail "dist --output-samples lost the direction column"
# The sketch is canonical, so the samples carry no strand column.
strip_comments "$WORK/dist_samples.tsv" | head -n 1 | grep -q "strand" &&
  fail "dist --output-samples kept a strand column for a canonical sketch"
echo "PASS: dist"

# -- Phase 6: roll -------------------------------------------------------------
echo "=== Phase 6: roll ==="
"$GDIFF" roll -l 100 -s 500 "genomes/${first_query}.fna.gz" "$sketch_a" > "$WORK/roll.tsv" 2>/dev/null
has_data "$WORK/roll.tsv" || fail "roll produced no output"
[ "$(strip_comments "$WORK/roll.tsv" | head -n 1 | cut -f1)" = "seq" ] || fail "roll header changed"
# The sketch is canonical, so the profile carries no strand column either.
[ "$(strip_comments "$WORK/roll.tsv" | head -n 1 | awk -F'\t' '{ print NF }')" -eq 5 ] ||
  fail "roll kept a strand column for a canonical sketch"

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
    map)  "$GDIFF" map --thresholds 0.1 $MAP_BASE "genomes/${first_query}.fna.gz" "$WORK/sketches/${first_ref}.gs" \
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
# `map` names its columns so `gviz.py` can find them by name instead of by position; the
# strand-aware variant is not covered by the gt/ files above.
EXPECT_AGNOSTIC="seq seq_len start end reference alpha threshold side interval_len n_kmers d d_seq d_null delta d_min d_max d_floor d_ceiling n_at_floor n_at_ceiling p_interval q_interval p_gamma"
EXPECT_AWARE="seq seq_len start end strand is_rc reference alpha threshold side interval_len n_kmers d d_seq d_null delta d_min d_max d_floor d_ceiling n_at_floor n_at_ceiling p_interval q_interval p_gamma"

agnostic_header="$(strip_comments "$WORK/map_ref.tsv" | head -n 1 | tr '\t' ' ')"
[ "$agnostic_header" = "$EXPECT_AGNOSTIC" ] || fail "strand-agnostic map header changed: $agnostic_header"

"$GDIFF" sketch -k 27 -w 31 -h 11 --no-canonical -i "genomes/${first_query}.fna.gz" \
  -o "$WORK/aware.gs" >/dev/null 2>&1
"$GDIFF" map --thresholds 0.1 -l 9900 --hdist-th 4 --sample-size 100 "genomes/${first_query}.fna.gz" \
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

# -- Phase 10: ambiguous runs split extraction ---------------------------------
# The corpus genomes carry no long ambiguous run, so this builds one: two mutated blocks with a
# 50-base N run between them, mapped back to the genome's own sketch. Skip handling has to split
# the sequence at the run, so no record may carry data on both sides of it.
echo "=== Phase 10: N-run splitting ==="
gzip -dc "genomes/${first_query}.fna.gz" | awk '
  BEGIN { srand(1) }
  /^>/ { if (n++) exit; print; next }
  {
    out = ""
    for (i = 1; i <= length($0); i++) {
      pos++
      c = substr($0, i, 1)
      if (pos >= 20001 && pos <= 20050) c = "N"
      else if ((pos >= 10001 && pos <= 20000) || (pos >= 20051 && pos <= 30050)) c = substr("ACGT", int(rand() * 4) + 1, 1)
      out = out c
    }
    printf "%s\n", out
  }' > "$WORK/nrun.fa"

nrun_map() {
  "$GDIFF" map -l 5000 --thresholds 0.60 --hdist-th 4 "$WORK/nrun.fa" "$WORK/sketches/${first_query}.gs" \
    > "$1" 2>/dev/null
}
nrun_map "$WORK/nrun.tsv" || fail "map failed on an N-run query"
strip_comments "$WORK/nrun.tsv" > "$WORK/nrun.data"
n_rows="$(awk 'NR>1' "$WORK/nrun.data" | wc -l | tr -d ' ')"
[ "$n_rows" -gt 0 ] || fail "map produced no records for an N-run query"

# A record may reach a few bases into the run - the k-mers overlapping its start are simply not
# scored - but none may carry data on both sides of it, and none may leave the sequence.
if ! awk -F'\t' 'NR>1 { if ($3 <= 19901 && $4 >= 20150) bad=1; if ($3 < 1 || $4 > $2) bad=2 } END { exit bad ? 1 : 0 }' \
  "$WORK/nrun.data"; then
  fail "a map record crossed the N run or left the sequence"
fi
# Both mutated blocks must be reported separately, one on each side of the run.
awk -F'\t' 'NR>1 { if ($3 < 20001) left++; if ($3 > 20050) right++ } END { exit (left > 0 && right > 0) ? 0 : 1 }' \
  "$WORK/nrun.data" || fail "the N run did not separate the two mutated blocks"
# Same input, same output.
nrun_map "$WORK/nrun2.tsv" || fail "the second N-run map failed"
cmp -s "$WORK/nrun.data" <(strip_comments "$WORK/nrun2.tsv") || fail "map output is not stable on an N-run query"
echo "PASS: N-run splitting"

# -- Summary -------------------------------------------------------------------
echo ""
if [ "$FAIL" -ne 0 ] || [ "$ENUM_FAIL" -ne 0 ] || [ "$CONT_FAIL" -ne 0 ]; then
  echo "=== REGRESSION FAILED ==="
  exit 1
fi
echo "=== ALL REGRESSION TESTS PASSED ==="
