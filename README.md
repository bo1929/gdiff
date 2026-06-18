# gdiff
gdiff is a k-mer-based tool for identifying genomic intervals in query sequences where the local k-mer distance to a reference differs significantly from the background.
It builds an LSH-based sketch of a reference genome, then scans each query sequence to detect contiguous intervals whose k-mer profile is closer (or further) to the reference than expected by chance, using a log-likelihood ratio tested with a chi-squared threshold.

## Installation

### Compiling from the source

Clone the repository with its submodules and compile with

```
git clone --recurse-submodules -j8 https://github.com/bo1929/gdiff.git
cd gdiff && make
```

Run `./gdiff --help` to verify. You may then copy the binary to a directory on your `$PATH` (e.g., `cp ./gdiff ~/.local/bin`).

## Documentation

- **[TODO.md](TODO.md)** -- active task list.
- **[docs/project-overview.md](docs/project-overview.md)** -- architecture and pipeline guide for agents.
- **[docs/engineering-review.md](docs/engineering-review.md)** -- algorithm, statistics, known gaps, and suggested next steps.
- **[docs/proof-extract-intervals-sx-mx.md](docs/proof-extract-intervals-sx-mx.md)** -- correctness notes for interval extraction.
- **[docs/plot-notes.md](docs/plot-notes.md)** -- gdiff TSV schema and plot.py integration notes.

## Quickstart

gdiff operates in two stages: build a sketch from a reference sequence, then map query sequences against it to extract intervals.

### 1. Building a sketch

```
gdiff sketch -i reference.fasta -o reference.sketch
```

This uses default parameters (`k=27`, `w=33`, `h=11`, `m=2`, `r=1`).
To tune the index:

```
gdiff sketch -i reference.fasta -o reference.sketch -k 27 -w 33 -h 11 -m 2 -r 1
```

The sketch is saved as a binary file. Multiple sketches can later be merged (see below).

### 2. Mapping query sequences

```
gdiff map -i reference.sketch -q query.fasta -d 0.05 -l 500 --sample-size 200 --num-threads 8 | tee intervals.tsv
```

`-d` sets the distance threshold and `-l` sets the minimum interval length (in base pairs).
`--sample-size 0` skips significance testing (enum lite when combined with `--enum-only`).

The output is a 15-column tab-separated file:

```
QUERY_ID  SEQ_LEN  INTERVAL_START  INTERVAL_END  STRAND  REF_ID  DIST  MASK  D_INTERVAL  DIST_CONTIG  STRAND_DIFF  DIST_GENOME  PERCENTILE  FOLD  QVALUE
read1     5200     120             3041          -       ref_A   0.04  1     (0, 0.05)   0.06         -0.02        0.05         0.003       0.8   0.012
```

| Column | Description |
|--------|-------------|
| `INTERVAL_START`, `INTERVAL_END` | 1-based sequence coordinates (see below) |
| `STRAND` | `'.'` when this strand has no contig MLE context; `'-'` on the lower-distance (reference) strand when both are known; `'+'` on the opposite strand or on the sole informative strand when only one contig MLE exists |
| `DIST` | MLE distance for the segment; NaN in enum lite |
| `MASK` | Bitmask of threshold index (`0` = background gap between segments) |
| `D_INTERVAL` | Parenthesized distance bin, e.g. `(0, 0.05)` |
| `DIST_CONTIG` | Contig-level MLE on this strand |
| `STRAND_DIFF` | `d_fw - d_rc` when both contig MLEs are finite; `-inf` when only fw is finite; `+inf` when only rc is finite; `NaN` when neither is finite |
| `DIST_GENOME` | Genome-wide MLE on the reference strand |
| `PERCENTILE`, `FOLD`, `QVALUE` | Significance vs null model; NaN when not tested (enum lite, full-contig rows, or too few null samples) |

Coordinates are 1-based. In `--enum-only` mode, `INTERVAL_START` and `INTERVAL_END` are inclusive sequence coordinates covering the selected k-mer bins. In continuous mode, rows partition the query by 1-based k-mer-start boundaries; adjacent rows can share a boundary coordinate, and the final row extends to `SEQ_LEN`.

By default, output goes to stdout. Use `-o intervals.tsv` to write to a file instead.

### 3. Merging sketches

Multiple sketches built from different references can be merged into a single file, which gdiff will then process in parallel:

```
gdiff merge -i ref_A.sketch ref_B.sketch ref_C.sketch -o all_refs.sketch
```

### 4. Inspecting a sketch

```
gdiff info -i reference.sketch
```

Prints the name, build date, k-mer length, window length, LSH parameters, subsampling rate (rho), and total k-mer count for each sketch in the file.

## Options

### `gdiff sketch`

| Option | Default | Description |
|---|---|---|
| `-i, --input-path` | — | Input FASTA/FASTQ file or URL (gzip compatible). **Required.** |
| `-o, --output-path` | — | Output binary sketch file. **Required.** |
| `-k, --mer-len` | `27` | k-mer length (19–32). |
| `-w, --win-len` | `k+6` | Minimizer window length (≥ k). |
| `-h, --num-positions` | `k-16` | Number of LSH positions (3–16). |
| `-m, --modulo-lsh` | `2` | Modulo value to partition LSH space. |
| `-r, --residue-lsh` | `1` | Include k-mer x if `LSH(x) mod m == r`. |
| `--frac / --no-frac` | `true` | Include k-mers with `LSH(x) mod m <= r`. |

### `gdiff map`

| Option | Default | Description |
|---|---|---|
| `-q, --query-path` | — | Query FASTA/FASTQ file or URL (gzip compatible). **Required.** |
| `-i, --sketch-path` | — | Sketch file to query. **Required.** |
| `-d, --dist-th` | — | Distance threshold(s); provide exactly 1 or 8 values. **Required.** |
| `-l, --min-length` | — | Minimum interval length (bp). **Required.** |
| `-o, --output-path` | stdout | Write output to a file. |
| `--hdist-th` | `4` | Maximum Hamming distance for a k-mer hit; supported range is 0–7. |
| `--chisq` | `33.00051` | Chi-squared threshold (default ≈ p=1e-10 with 8 d.f.). |
| `-b, --bin-shift` | `0` | Group consecutive k-mers into bins of size 2^b. |
| `--sample-size` | `200` | Null samples for significance test (`0` skips the test). |
| `--num-threads` | `1` | Number of threads for parallel sketch processing. |
| `--enum-only` | off | Simple per-threshold intervals instead of continuous breakpoint merging. With `--sample-size 0`, emits coordinates only (no MLE or significance). |

When exactly 8 distance thresholds are provided via `-d`, gdiff evaluates all eight in one SIMD-wide pass (`cm512_t` path).

Continuous mode (default) runs segment MLE, genome-wide MLE, and significance testing when `--sample-size > 0`. See [docs/engineering-review.md](docs/engineering-review.md) for statistical details.

## Interactive visualization

The `plot.py` script builds an interactive Dash/Plotly app for intervals (and optional tree / annotations). It auto-detects the TSV format (legacy 7-column enum, 15-column enum lite, or continuous).

```bash
python plot.py --input intervals.tsv --tree tree.nwk --query query_name
```

Legacy enum-only files (with `DIST_TH` instead of `DIST`):

```bash
python plot.py --input enum.tsv --tree tree.nwk --enum-only
```

With optional annotations:

```bash
python plot.py --input intervals.tsv --tree tree.nwk --annotation annot.tsv --query query_name
```

See `python plot.py --help` for options (port, host, debug, etc.).
