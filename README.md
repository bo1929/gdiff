# gdiff

Detect genomic regions where the evolutionary distance between a query and a reference differs from the genome-wide background. gdiff builds a compact sketch of a reference genome, then scans query sequences to find intervals that are significantly closer or more distant than expected.

## Installation

```bash
git clone --recurse-submodules -j8 https://github.com/bo1929/gdiff.git
cd gdiff && make
./gdiff --help
```

Requires a C++17 compiler (clang or gcc), `libz`, and `libcurl`. The binary is self-contained — copy `./gdiff` anywhere on your `$PATH`.

## Quickstart

gdiff works in two steps: **sketch** a reference, then **map** queries against it.

### 1. Sketch a reference genome

```bash
gdiff sketch -i reference.fasta -o reference.gs
```

Multiple references can be sketched into one container in a single step (use `--num-threads` to sketch files in parallel; save order is not fixed):

```bash
gdiff --num-threads 8 sketch -i ref_A.fasta ref_B.fasta ref_C.fasta -o combined.gs
```

Uses sensible defaults (`k=27`, `w=33`, `h=11`). For large genomes, tune the LSH parameters to trade speed for sensitivity (see Options below).

### 2. Map queries to find divergent intervals

```bash
gdiff map queries.fasta reference.skc -d 0.05 -l 500 --sample-size 200 -o intervals.tsv
```

`-d 0.05` looks for regions whose MLE distance departs from the query's own distance toward 0.05: thresholds below the query distance flag closer regions, thresholds above it flag more distant ones. `-l 500` sets the minimum interval to 500 k-mers. `--sample-size 200` controls how many background windows are sampled for significance testing.

You can provide multiple thresholds (up to 8) in one pass:

```bash
gdiff map queries.fasta reference.skc -d 0.02 0.05 0.10 0.20 0.30 0.40 0.50 0.60 -l 500
```

Or derive thresholds from the background null as tail probabilities:

```bash
gdiff map queries.fasta reference.skc --levels 0.05 0.01 -l 500
```

To skip significance testing and just enumerate intervals:

```bash
gdiff map queries.fasta reference.skc -d 0.05 -l 500 --enum-only --sample-size 0
```

### 3. Merge sketches (optional)

Combine sketches from multiple references into one file; gdiff maps against all of them in parallel. Prefer `gdiff sketch -i a.fa b.fa -o combined.gs` when starting from FASTA:

```bash
gdiff merge -i ref_A.gs ref_B.gs ref_C.gs -o combined.gs
```

### 4. Inspect a sketch

```bash
gdiff info -i reference.gs
```

### 5. Sample fixed-length distances

Sample exact-length query regions and summarize their MLE distances:

```bash
gdiff dist -i reference.gs -q queries.fasta -l 500 --sample-size 200 -o distance_summary.tsv
```

Use `--output-samples` to write every valid sampled region and its distance
instead of the per-reference summary (to stdout, or the file given by `-o`). The sample size applies to the entire query file, with
each eligible sequence selected in proportion to its length via weighted
reservoir sampling. Histograms are built only for sequences that claim sample
slots and are discarded afterward. Use `-b/--bin-shift` to bin k-mers (as in
`map`) and `--num-threads` to process references in parallel. Sampling is with
replacement and is controlled by the global `--seed` option.

## Output format

The default (continuous) output is a tab-separated file with these columns:

```
QUERY_ID  SEQ_LEN  INTERVAL_START  INTERVAL_END  STRAND  IS_RC  REF_ID  DIST  MASK  D_INTERVAL  DIST_CONTIG  STRAND_DIFF  DIST_GENOME  PERCENTILE  FOLD  QVALUE  INFO  LR_UB
```

For strand-agnostic references the three strand columns (`STRAND`, `IS_RC`, `STRAND_DIFF`) are omitted, giving 15 columns. An example strand-aware row:

```
contig1  5200  120  3041  +  0  ref_A  0.04  1  (0, 0.05)  0.06  -0.02  0.05  0.003  0.8  0.012  8.3e+06  1200
```

| Column | Meaning |
|--------|--------|
| `QUERY_ID` | Query sequence name |
| `SEQ_LEN` | Query length in base pairs |
| `INTERVAL_START`, `INTERVAL_END` | 1-based coordinates of the detected interval |
| `STRAND` | `+` = closer (lower distance), `-` = farther, `.` = unknown |
| `IS_RC` | `0` = forward strand, `1` = reverse-complement |
| `REF_ID` | Reference genome the interval was found against |
| `DIST` | MLE evolutionary distance for this interval |
| `MASK` | Which distance threshold(s) triggered this interval (bitmask) |
| `D_INTERVAL` | Tightest bounds among satisfied thresholds, e.g. `(0, 0.05)` |
| `DIST_CONTIG` | MLE distance of the entire query contig on this strand |
| `STRAND_DIFF` | Difference between forward and reverse-complement contig distances |
| `DIST_GENOME` | Genome-wide average distance across all queries |
| `PERCENTILE` | Empirical significance of the interval distance within the background null (NaN = not tested) |
| `FOLD` | Fold change: interval distance divided by null median |
| `QVALUE` | Benjamini-Hochberg adjusted p-value (per strand) |
| `INFO` | Observed Fisher information of the interval at `DIST` |
| `LR_UB` | Likelihood-ratio statistic of `DIST` against the maximum estimable distance |

Coordinates are 1-based and inclusive. Adjacent rows may share a boundary; the final row on each query strand spans to `SEQ_LEN`.

In `--enum-only` mode, each row is an independent interval covering the k-mer bins that satisfy the threshold.

## Options

### `gdiff sketch`

| Option | Default | Description |
|--------|--------|-------------|
| `-i, --input-path` | (required) | Input FASTA/FASTQ file(s) or URL (gzip compatible) |
| `-o, --output-path` | (required) | Output container (one or more sketches) |
| `-k, --mer-len` | `27` | k-mer length (19–32) |
| `-w, --win-len` | `k+6` | Minimizer window length (>= k) |
| `-h, --num-positions` | `max(floor(k/2)-2, k-16)` | Number of LSH positions (3–16) |
| `--frac` | `1.0` | Keep k-mer if LSH(x) < frac · 2^2h; subsamples on top of minimizers |
| `--num-threads` | `1` | Parallel input-file sketching threads |

### `gdiff map`

Positional: `<query.fasta> <reference.skc>`.

| Option | Default | Description |
|--------|--------|-------------|
| `-d, --dist-th` | - | Absolute distance threshold(s), up to 8; required unless `--levels` |
| `--levels` | - | Two-sided tail probabilities (up to 4); thresholds are their empirical quantiles of the background; required unless `-d` |
| `-l` | (required) | Minimum interval length in k-mers; also the background window length |
| `-o, --output-path` | stdout | Write output to a file |
| `--hdist-th` | `3` | Max Hamming distance for a k-mer hit (0–7) |
| `--chisq` | `33.00051` | Chi-squared threshold for merging adjacent intervals |
| `-b, --bin-shift` | `0` | Bin size = 2^b; groups consecutive k-mers |
| `--sample-size` | `200` | Background windows sampled per reference (`0` = skip significance) |
| `--enum-only` | off | Simple per-threshold enumeration instead of ordered removal |
| `--per-sequence` | off | Measure significance against a per-query background |
| `--verbosity` | `1` | Background and threshold report detail (0–1) |
| `--num-threads` | `1` | Background-sampling and per-sequence scan threads |

Distance thresholds may be given either directly (`-d`) or as tail probabilities of the background (`--levels`); the two are mutually exclusive. Providing more than one threshold runs them in one SIMD-wide pass.

The background null is the empirical distribution of `-l`-long query windows scored against the reference. `PERCENTILE` and `FOLD` are the empirical rank and median ratio within that null; no parametric fit is used.

### `gdiff dist`

| Option | Default | Description |
|--------|--------|-------------|
| `-q, --query-path` | (required) | Query FASTA/FASTQ file or URL (gzip compatible) |
| `-i, --sketch-path` | (required) | Container to query against |
| `-l, --length` | (required) | Sampled region length in k-mers |
| `-b, --bin-shift` | `0` | Bin size = 2^b; groups consecutive k-mers |
| `--sample-size` | `200` | Regions sampled across the whole query file per reference |
| `--hdist-th` | `4` | Max Hamming distance for a k-mer hit (0-7) |
| `-o, --output-path` | stdout | Write summary output to a file |
| `--output-samples` | off | Write per-sample output instead of the per-reference summary |
| `--num-threads` | `1` | Parallel sketch/reference processing threads |

Summary rows contain:

```
QUERY_FILE  REF_ID  N  D_MED
```

`N` is the number of valid MLE samples. `D_MED` is the median of those
distances. Sampled windows with no matching k-mer (zero hits within
`--hdist-th`) have an undefined distance: they are excluded from the median
and their count is reported on stderr. Sample detail rows contain one row per
sampled window (unmapped windows carry NaN fields):

```
config  genome_a  genome_b  seq  start  end  strand  direction  d  lr_ub
```

`LR_UB` is the likelihood-ratio statistic of the sample's distance against the
sketch's max estimable distance.

Coordinates are 1-based and inclusive. For strand-aware sketches, each sample
uses the lower of the forward and reverse-complement MLE distances and reports
the selected strand.

## Interactive visualization

Use `plot.py` to explore intervals interactively:

```bash
pip install -r plot-reqs.txt
python plot.py --input intervals.tsv --tree tree.nwk --query query_name
```

Supports optional genome annotations (GFF3/GTF/TSV) and legacy enum-only files (`--enum-only`). Run `python plot.py --help` for all options.
