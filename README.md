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

Multiple references can be sketched into one container in a single step (use `--num-threads` to sketch files in parallel; save order is not fixed):

```bash
gdiff --num-threads 8 sketch -i ref_A.fasta ref_B.fasta ref_C.fasta -o combined.skc
```

Uses sensible defaults (`k=23`, `w=23`, `h=9`). For large genomes, tune the LSH parameters to trade speed for sensitivity (see Options below).

### 2. Map queries to find divergent intervals

```bash
gdiff map queries.fasta reference.skc -d 0.05 -l 500 --sample-size 200 -o intervals.tsv
```

`-d 0.05` looks for regions whose MLE distance departs from the query's own distance toward 0.05: thresholds below the query distance flag closer regions, thresholds above it flag more distant ones. `-l 500` sets the minimum interval to 500 k-mers. `--sample-size 200` controls how many background windows are sampled for significance testing.

You can provide multiple thresholds (up to 8) in one pass:

```bash
gdiff map queries.fasta reference.skc -d 0.02 0.05 0.10 0.20 0.30 0.40 0.50 0.60 -l 500
```

Without `-d`, thresholds come from the background null as tail probabilities (the default is `0.1 0.05 0.01 0.005`):

```bash
gdiff map queries.fasta reference.skc -l 500
gdiff map queries.fasta reference.skc --levels 0.2 0.1 0.05 0.02 -l 500
```

To skip significance testing and just enumerate intervals:

```bash
gdiff map queries.fasta reference.skc -d 0.05 -l 500 --enum-only --sample-size 0
```

### 3. Merge sketches (optional)

Combine sketches from multiple references into one file; gdiff maps against all of them in parallel. Prefer `gdiff sketch -i a.fa b.fa -o combined.skc` when starting from FASTA:

```bash
gdiff merge -i ref_A.skc ref_B.skc ref_C.skc -o combined.skc
```

### 4. Inspect a sketch

```bash
gdiff info -i reference.skc
```

### 5. Compare sketches (optional)

`dist` compares the sampled windows stored in sketches, so both inputs must have been created with `-l`:

```bash
gdiff dist ref_A.skc ref_B.skc -o distances.tsv
gdiff dist combined.skc                 # every pair within one container
gdiff dist --list-a a.txt --list-b b.txt
```

Use `--output-samples` to write one row per sampled window instead of the per-pair summary. `--num-threads` processes pairs in parallel; the global `--seed` controls sampling at sketch time.

## Options

### `gdiff sketch`

| Option | Default | Description |
|--------|--------|-------------|
| `-i, --input-path` | - | Input FASTA/FASTQ file(s) or URL (gzip compatible) |
| `--input-list` | - | Read input paths from a file, one per line (optional `name<TAB>path`); combines with `-i` |
| `-o, --output-path` | (required) | Output container (one or more sketches) |
| `-k, --mer-len` | `23` | k-mer length (19–31) |
| `-w, --win-len` | `k` | Minimizer window length (>= k) |
| `-h, --num-positions` | `max(floor(k/2)-2, k-16)` | Number of LSH positions |
| `--frac` | `0.5` | Keep a k-mer if LSH(x) < frac · 2^2h; subsampling ratio |
| `--canonical` / `--no-canonical` | `--canonical` | Canonical (strand-agnostic) k-mers, or keep forward/reverse separate |
| `-l` | `333` | Sampled window length in k-mers; **`0` stores buckets only** |
| `--sample-size` | `1000` | Windows sampled across each genome; capped by the eligible starts |
| `--keep-seq` | off | Store sampled windows as 2-bit packed bases instead of pre-resolved keys |

`-l 0` produces a buckets-only sketch. `map` and `roll` work with it (they only need the reference index); `dist` needs sampled windows. By default sketches carry windows.

A window of `-l L` spans `L + k - 1` bases, so a sequence shorter than that contributes no windows and `--sample-size` is an upper bound rather than a promise: the draw cannot exceed the eligible starts across the input. `gdiff sketch` warns when it stores fewer than requested and `gdiff info -i <container>` reports the number actually stored (`Windows:`), which is also the number of rows `dist --output-samples` emits per direction.

### `gdiff map`

Positional: `<query.fasta> <reference.skc>`.

| Option | Default | Description |
|--------|--------|-------------|
| `-d, --dist-th` | - | Absolute distance threshold(s): exactly 1 or 8; overrides `--levels` |
| `--levels` | `0.1 0.05 0.01 0.005` | Two-sided tail probabilities: exactly 4; thresholds are their empirical quantiles |
| `-l` | (required) | Minimum interval length in k-mers; also the background window length |
| `-o, --output-path` | stdout | Write output to a file |
| `--hdist-th` | `3` | Max Hamming distance for a k-mer hit (0–7) |
| `--chisq` | `33.00051` | Chi-squared threshold for merging adjacent intervals |
| `-b, --bin-shift` | `0` | Bin size = 2^b; groups consecutive k-mers |
| `--sample-size` | `500` | Background windows sampled per reference (`0` = skip significance) |
| `--enum-only` | off | Simple per-threshold enumeration instead of ordered removal |
| `--per-sequence` | off | Measure significance against a per-query background |
| `--verbosity` | `1` | Background and threshold report detail (0–1) |

Distance thresholds may be given directly (`-d`, exactly 1 or 8 values) or as tail probabilities of the background (`--levels`, exactly 4 levels; the default). The two are mutually exclusive, and either way the run fills its SIMD lanes exactly: `-d` with 1 value uses the scalar path, everything else uses 8 lanes. Because the lower empirical quantiles often sit at the estimable floor (background windows that match almost perfectly), a floored lower quantile is lifted to the nearest higher distinct background distance and each later level takes the next distinct one; a warning says when that happens. If the background is too coarse to place eight distinct thresholds the run aborts, and `--sample-size` sets how finely the quantiles resolve.

The background null is the empirical distribution of `-l`-long query windows scored against the reference. `PERCENTILE` and `FOLD` are the empirical rank and median ratio within that null; no parametric fit is used. `map` only reads the reference's bucket index, so a buckets-only sketch (`-l 0`) is sufficient.

### `gdiff dist`

Positional: `<sketch-a> [sketch-b]`.

| Option | Default | Description |
|--------|--------|-------------|
| `--list-a`, `--list-b` | - | Read set A / set B as a list file of sketch containers |
| `--hdist-th` | `3` | Maximum Hamming distance for a k-mer to match; capped at 2 for inputs >20 Mbp |
| `--lr-th` | `10.828` | Likelihood-ratio cut for the reconciliation filter |
| `--min-portion` | `0.66` | Apply the filter only if at least this fraction of windows exceeds `--lr-th` |
| `-o, --output-path` | stdout | Write output to a file |
| `--output-samples` | off | Write per-window sample rows instead of per-pair summaries |

Inputs are containers whose sketches were created with `-l` (sampled windows). One positional compares every pair within that container; two compare the cross product of set A and set B. Pairs are reconciled across both directions.

Summary output, one row per pair:

```
genome_a  genome_b  d  d_median  d_mean  d_upper  d_highest  d_ab  d_ba  n_ab  n_ba  n_ub  n_na  n_filtered
```

`d` is the reconciled distance and `d_ab`/`d_ba` the two directional medians. `--output-samples` instead writes one row per sampled window:

```
config  genome_a  genome_b  seq  start  end  strand  direction  d  lr_ub
```

`LR_UB` is the likelihood-ratio statistic of the sample's distance against the sketch's maximum estimable distance. Coordinates are 1-based and inclusive; for strand-aware sketches each sample uses the lower of the forward and reverse-complement distances and reports the selected strand.

### `gdiff roll`

Positional: `<query.fasta> <reference.skc>`.

| Option | Default | Description |
|--------|--------|-------------|
| `-l` | (required) | Window length in k-mers |
| `-s` | `-l` | Step between consecutive window starts |
| `--hdist-th` | `3` | Maximum Hamming distance for a k-mer to match |
| `-o, --output-path` | stdout | Write output to a file |

### `gdiff merge`

| Option | Default | Description |
|--------|--------|-------------|
| `-i, --sketch-paths` | (required) | Input containers to merge |
| `-o, --output-path` | (required) | Path to store the merged container |

### `gdiff info`

| Option | Default | Description |
|--------|--------|-------------|
| `-i, --sketch-path` | (required) | Container to inspect |

### Global options

| Option | Default | Description |
|--------|--------|-------------|
| `--verbose` | off | Report progress even when stderr is not a terminal |
| `--seed` | `0` | Random seed for the LSH and other randomness |
| `--num-threads` | `1` | Worker threads per subcommand |

## Visualization

Two tools, two jobs.

### `gviz.py` — one report for any output

`gviz.py` turns a single `map`, `roll` or `dist` result file into a static HTML
page: two or three panels and a filter bar. The panels never zoom, pan or hover —
changing a filter redraws them, which is what keeps the page fast. Plotly.js comes
from the CDN; `--inline` embeds it.

```bash
python gviz.py map-mus_pwk-v020rc.tsv            # 2 panels, 2.4 MB
python gviz.py roll-mus_pwk-v020rc.tsv           # 3 panels, 3M windows in 4.5 s
python gviz.py roll.tsv --seq chr7 --bins 6000   # one chromosome in fine bins, 0.7 MB
python gviz.py dist-m1-v020rc.tsv --annotations genes.gff3 --open
python gviz.py map.tsv --print-schema
```

Panels, position on x wherever the output has coordinates:

| Input | Panels |
|-------|--------|
| `map` | interval track — one row per reference and direction, colour per reference, area ∝ length · significance along the sequence (-log10 q, q = 0.05 line) |
| `roll` | distance profile per reference, binned and smoothed with a running median · reference x position heat map · unmapped-window map |
| `dist` | pairwise distance matrix — square tiles, rows and columns in one shared order (nearest genomes adjacent), mirrored for all-by-all, values printed in the cells, empty tiles for pairs the file does not contain · estimator spread per pair · window accounting |
| `dist --output-samples` | per-window profile · per-pair ECDF |

Reading the report:

- **Axes are pinned** to the full data extent (log axes to the enclosing powers of
  ten), so filtering never rescales a plot or moves a label. Only the Sequence
  selector changes the x range, because each sequence has its own length; `‹` and
  `›` step through sequences, which are listed longest first.
- **One sequence at a time** — the Sequence filter is a single select, so a 2.5 Gb
  genome never becomes a smear. Sequences below 0.2 % of the assembly are left out
  unless you pass `--keep-all`.
- **Thresholds are named, not numbered in CLI terms**: the Threshold filter lists
  e.g. `3 · closer d ≤ 0.000344`, reconstructed from the bracket each interval fell
  into, and holds any combination.
- **Colour** encodes a value (distance, -log10 q) with viridis; when colour encodes
  a series it is the reference, in muted lichen tones.
- `roll` drops unmapped windows instead of drawing gaps, averages each sequence
  into ~1000 bins and smooths with a running median (`smoothing` control, default
  5 points) so spikes do not drag the line; `log y` is a checkbox.
- `--annotations genes.gff3` adds a gene track; a whole chromosome of genes is
  dense by nature, so at most the 3,000 longest features are drawn and `--seq`
  is the way to zoom in.

Options: `--seq`, `--bins`, `--keep-all`, `--assembly-report`, `--annotations`,
`--title`, `--inline`, `--kind`, `--print-schema`. Adding an output format means
one branch in `kind_of`, a `prepare_*` and a `*_panels` function.

### `plot.py` — phylogeny-aware interval explorer

Use `plot.py` when a tree is part of the question:

```bash
pip install -r plot-reqs.txt
python plot.py --input intervals.tsv --tree tree.nwk --query query_name
```

Supports optional genome annotations (GFF3/GTF/TSV) and legacy enum-only files (`--enum-only`). Run `python plot.py --help` for all options.
