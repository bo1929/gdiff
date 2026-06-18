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
gdiff sketch -i reference.fasta -o reference.skc
```

Uses sensible defaults (`k=27`, `w=33`, `h=11`). For large genomes, tune the LSH parameters to trade speed for sensitivity (see Options below).

### 2. Map queries to find divergent intervals

```bash
gdiff map -i reference.skc -q queries.fasta -d 0.05 -l 500 --sample-size 200 -o intervals.tsv
```

`-d 0.05` looks for regions whose MLE distance is below 0.05 (closer than expected). Use a negative value like `-d -0.05` to find regions that are more distant. `-l 500` sets the minimum interval to 500 bp. `--sample-size 200` controls how many genome-wide null samples are used for significance testing.

You can provide multiple thresholds (up to 8) in one pass:

```bash
gdiff map -i ref.skc -q queries.fasta -d 0.02 0.05 0.10 0.20 0.30 0.40 0.50 0.60 -l 500
```

To skip significance testing and just enumerate intervals:

```bash
gdiff map -i ref.skc -q queries.fasta -d 0.05 -l 500 --enum-only --sample-size 0
```

### 3. Merge sketches (optional)

Combine sketches from multiple references into one file; gdiff maps against all of them in parallel:

```bash
gdiff merge -i ref_A.skc ref_B.skc ref_C.skc -o combined.skc
```

### 4. Inspect a sketch

```bash
gdiff info -i reference.skc
```

## Output format

The default (continuous) output is a tab-separated file with these columns:

```
QUERY_ID  SEQ_LEN  INTERVAL_START  INTERVAL_END  STRAND  IS_RC  REF_ID  DIST  MASK  D_INTERVAL  DIST_CONTIG  STRAND_DIFF  DIST_GENOME  PERCENTILE  FOLD  QVALUE
```

An example row:

```
read1  5200  120  3041  +  ref_A  0.04  1  (0, 0.05)  0.06  -0.02  0.05  0.003  0.8  0.012
```

| Column | Meaning |
|--------|--------|
| `QUERY_ID` | Query sequence name |
| `SEQ_LEN` | Query length in base pairs |
| `INTERVAL_START`, `INTERVAL_END` | 1-based coordinates of the detected interval |
| `STRAND` | `+` = closer (lower distance), `-` = farther, `.` = unknown |
| `IS_RC` | `0` = forward strand, `1` = reverse-complement |
| `REF_ID` | Reference genome the interval was found against |
| `DIST` | MLE evolutionary distance for this interval (NaN in enum lite mode) |
| `MASK` | Which distance threshold(s) triggered this interval (bitmask) |
| `D_INTERVAL` | Tightest bounds among satisfied thresholds, e.g. `(0, 0.05)` |
| `DIST_CONTIG` | MLE distance of the entire query contig on this strand |
| `STRAND_DIFF` | Difference between forward and reverse-complement contig distances |
| `DIST_GENOME` | Genome-wide average distance across all queries |
| `PERCENTILE` | Significance p-value from the null model (NaN = not tested) |
| `FOLD` | Fold change: interval distance divided by null median |
| `QVALUE` | Benjamini-Hochberg adjusted p-value (per strand) |

Coordinates are 1-based and inclusive. Adjacent rows may share a boundary; the final row on each query strand spans to `SEQ_LEN`.

In `--enum-only` mode, each row is an independent interval covering the k-mer bins that satisfy the threshold.

## Options

### `gdiff sketch`

| Option | Default | Description |
|--------|--------|-------------|
| `-i, --input-path` | (required) | Input FASTA/FASTQ file or URL (gzip compatible) |
| `-o, --output-path` | (required) | Output sketch file |
| `-k, --mer-len` | `27` | k-mer length (19–32) |
| `-w, --win-len` | `k+6` | Minimizer window length (>= k) |
| `-h, --num-positions` | `k-16` | Number of LSH positions (3–16) |
| `-m, --modulo-lsh` | `2` | LSH space partition modulus |
| `-r, --residue-lsh` | `1` | Keep k-mer if LSH(x) mod m == r |
| `--frac / --no-frac` | `true` | Keep k-mer if LSH(x) mod m <= r |

### `gdiff map`

| Option | Default | Description |
|--------|--------|-------------|
| `-q, --query-path` | (required) | Query FASTA/FASTQ file or URL (gzip compatible) |
| `-i, --sketch-path` | (required) | Sketch file to query against |
| `-d, --dist-th` | (required) | Distance threshold(s); 1 or 8 values |
| `-l, --min-length` | (required) | Minimum interval length in base pairs |
| `-o, --output-path` | stdout | Write output to a file |
| `--hdist-th` | `4` | Max Hamming distance for a k-mer hit (0–7) |
| `--chisq` | `33.00051` | Chi-squared threshold for merging adjacent intervals |
| `-b, --bin-shift` | `0` | Bin size = 2^b; groups consecutive k-mers |
| `--sample-size` | `200` | Null samples for significance testing (`0` = skip) |
| `--num-threads` | `1` | Parallel sketch processing threads |
| `--enum-only` | off | Simple per-threshold enumeration; with `--sample-size 0` emits coordinates only |

Providing 8 distance thresholds runs all eight in one SIMD-wide pass.

## Interactive visualization

Use `plot.py` to explore intervals interactively:

```bash
pip install -r plot-reqs.txt
python plot.py --input intervals.tsv --tree tree.nwk --query query_name
```

Supports optional genome annotations (GFF3/GTF/TSV) and legacy enum-only files (`--enum-only`). Run `python plot.py --help` for all options.
