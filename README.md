# gidiff
gidiff is a k-mer-based tool for identifying genomic intervals in query sequences where the local k-mer distance to a reference differs significantly from the background.
It builds an LSH-based sketch of a reference genome, then scans each query sequence to detect contiguous intervals whose k-mer profile is closer (or further) to the reference than expected by chance, using a log-likelihood ratio tested with a chi-squared threshold.

## Installation

### Compiling from the source

Clone the repository with its submodules and compile with

```
git clone --recurse-submodules -j8 https://github.com/bo1929/gidiff.git
cd gidiff && make
```

Run `./gidiff --help` to verify. You may then copy the binary to a directory on your `$PATH` (e.g., `cp ./gidiff ~/.local/bin`).

## Quickstart

gidiff operates in two stages: build a sketch from a reference sequence, then map query sequences against it to extract intervals.

### 1. Building a sketch

```
gidiff sketch -i reference.fasta -o reference.sketch
```

This uses default parameters (`k=27`, `w=33`, `h=11`, `m=2`, `r=1`).
To tune the index:

```
gidiff sketch -i reference.fasta -o reference.sketch -k 27 -w 33 -h 11 -m 2 -r 1
```

The sketch is saved as a binary file. Multiple sketches can later be merged (see below).

### 2. Mapping query sequences

```
gidiff map -i reference.sketch -q query.fasta -d 0.05 -l 500 --num-threads 8 | tee intervals.tsv
```

`-d` sets the distance threshold and `-l` sets the minimum interval length (in base pairs).
The output is tab-separated:

```
QUERY_ID  SEQ_LEN  INTERVAL_START  INTERVAL_END  STRAND  REF_ID  DIST_TH
read1     5200     120             3041          +       ref_A   0.05
read2     4800     0               2199          -       ref_A   0.05
```

Each row is a genomic interval in the query whose k-mer profile is significantly closer to the reference than the background, as determined by a chi-squared test on the log-likelihood difference.
Coordinates are 0-based; `INTERVAL_END` is inclusive and extends to the last base of the final k-mer.

By default, output goes to stdout. Use `-o intervals.tsv` to write to a file instead.

### 3. Merging sketches

Multiple sketches built from different references can be merged into a single file, which gidiff will then process in parallel:

```
gidiff merge -i ref_A.sketch ref_B.sketch ref_C.sketch -o all_refs.sketch
```

### 4. Inspecting a sketch

```
gidiff info -i reference.sketch
```

Prints the name, build date, k-mer length, window length, LSH parameters, subsampling rate (rho), and total k-mer count for each sketch in the file.

## Options

### `gidiff sketch`

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

### `gidiff map`

| Option | Default | Description |
|---|---|---|
| `-q, --query-path` | — | Query FASTA/FASTQ file or URL (gzip compatible). **Required.** |
| `-i, --sketch-path` | — | Sketch file to query. **Required.** |
| `-d, --dist-th` | — | Distance threshold(s); provide exactly 1 or 8 values. **Required.** |
| `-l, --min-length` | — | Minimum interval length (bp). **Required.** |
| `-o, --output-path` | stdout | Write output to a file. |
| `--hdist-th` | `4` | Maximum Hamming distance for a k-mer hit. |
| `--chisq` | `33.00051` | Chi-squared threshold (default ≈ p=1e-10 with 8 d.f.). |
| `-b, --bin-shift` | `0` | Group consecutive k-mers into bins of size 2^b. |
| `--num-threads` | `1` | Number of threads for parallel sketch processing. |

When exactly 8 distance thresholds are provided via `-d`, gidiff uses AVX-512 SIMD to evaluate all 8 thresholds simultaneously in a single pass over the query.

## Interactive Visualization

The `plot.py` script creates an interactive Dash/Plotly web application for visualizing genomic intervals and phylogenetic trees.

### Basic Usage
```bash
python plotx.py --input <TSV_DATA> --tree <NEWICK_TREE> --query <QUERY_NAME>
```

### With Gene Annotations
Add a gene annotation panel below the interval visualization:
```bash
python plotx.py \
  --input <TSV_DATA> \
  --tree <NEWICK_TREE> \
  --annotation <ANNOTATION_TSV> \
  --query <QUERY_NAME>
```

### Options
- `--input`, `-i`: Path to TSV data file (required)
- `--tree`, `-t`: Path to Newick tree file (required)
- `--query`, `-q`: Query leaf name for distance calculations (optional)
- `--annotation`, `-a`: Path to annotation TSV file (optional)
- `--port`, `-p`: Port number (default: 8080)
- `--host`: Host address (default: 127.0.0.1)
- `--debug`: Enable debug mode
>>>>>>> 1aac9c6 (updated plotting to add gtf/gff annotations)
