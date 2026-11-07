# gdiff

`gdiff` estimates how sequence distance varies along a genome without alignment.
It compares a query sequence against a reference sketch and uses a pseudo-likelihood model to estimate distances and find unusually close or divergent regions.

The same sketch supports three analyses:

- `dist`: randomly sample windows and estimate distances, summarize the sample into a genome-wide estimate.
- `roll`: estimate local distances in sliding windows -- very fast!
- `map`: map intervals with distances that deviate from the rest of the genome (e.g., horizontal transfers).

## Installation

Requires a C++17 compiler (clang or GCC), Make, and zlib.
libcurl is optional and enables URL input when available at build time.

```bash
git clone --recurse-submodules -j8 https://github.com/bo1929/gdiff.git
cd gdiff
make
./gdiff --help
```

Add the resulting `gdiff` binary to your `PATH` to use the commands below.

## Demo

The [demo walkthrough](demo/) introduces `sketch`, `dist`, `roll`, and `map` with commands you can run one at a time.

## Usage

### Sketch a reference

```bash
gdiff sketch -i reference.fasta -o reference.gds
```

Sketches use canonical k-mers by default (`k=23`, `w=23`, `h=9`, `--frac 0.5`) and store up to 500 sampled windows of 333 k-mers for `dist`.
Use `-l 0` to omit the windows if you only need `map` or `roll`.

Multiple genomes can share one container. Each input file becomes a separate sketch:

```bash
gdiff --num-threads 8 sketch -i ref_A.fasta ref_B.fasta ref_C.fasta -o combined.gds
```

### Pairwise & genome-wide distance comparison

```bash
gdiff dist ref_A.gds ref_B.gds -o distances.tsv
gdiff dist combined.gds -o all_pairs.tsv
gdiff dist ref_A.gds ref_B.gds --output-samples -o samples.tsv
```

Two containers compare every genome in one against every genome in the other; one container compares all pairs within it. Both inputs must contain sampled windows. `--list-a` and `--list-b` accept lists of container paths.

The default output has one row per pair, with a reconciled distance `d`, directional medians `d_ab` and `d_ba`, and window counts.
`--output-samples` reports individual windows in both directions.
Distance is on the `1 - ANI/100` scale.

FASTA/FASTQ and gzip input are supported.
For large sets of inputs, you could also use `--input-list paths.txt`, with one path per line or `name<TAB>path` to assign ref. names.

### Estimate local distances through a sliding window
A standard sliding window view:
```bash
gdiff roll query.fasta reference.gds -l 1000 -s 200 -o windows.tsv
```

Here, `-l` is the window length in k-mers; `-s` is the step between starts (default: `-l`).
Make sure to use multi-threading if you have resources available.

### Find outlier intervals

```bash
gdiff map query.fasta reference.gds -l 500 -o intervals.tsv
```

`-l 500` sets the minimum interval length and background sampling window to 500 k-mers. gdiff samples query windows against each reference, fits a gamma background with a point mass at the distance floor, and derives lower and upper thresholds. If the fit fails, it uses empirical quantiles.

The default two-sided levels are `0.05 0.01 0.005`. For example, `0.05` places thresholds at the 2.5th and 97.5th percentiles. To change the levels or sample more background windows:

```bash
gdiff map query.fasta reference.gds -l 500 --levels 0.01 0.005 --sample-size 5000 -o intervals.tsv
```

On each side, the most extreme level claims intervals first; later levels search the remaining gaps. Use `--independent-levels` to enumerate each level separately, allowing overlapping calls. A container with multiple references produces results for each reference.

### Inspect or merge sketches

```bash
gdiff info -i reference.gds
gdiff merge -i ref_A.gds ref_B.gds -o combined.gds
```

`info` reports sketch parameters and stored window counts. `merge` requires compatible sketch configurations; when starting from sequences, pass multiple inputs to `sketch` instead.

## Output

Results are tab-separated, with a header and provenance comments. Coordinates are 1-based and inclusive. `map` reports query coordinates, not corresponding positions in the reference.

Key `map` columns:

| Column | Meaning |
|--------|---------|
| `seq`, `reference`, `start`, `end` | Query sequence, reference sketch, and interval coordinates |
| `alpha`, `threshold`, `side` | Background level, distance cutoff, and low/high side |
| `d`, `d_seq`, `d_null` | Interval distance, whole-sequence distance, and background median |
| `delta` | `d - d_null` |
| `p_interval`, `q_interval` | Likelihood-ratio p-value against `d_null` and its BH adjustment across reported candidates |
| `p_gamma` | Two-sided tail probability under the fitted background model |

`roll` reports one distance per window, or separate forward and reverse-complement distances for sketches made with `--no-canonical`.

## Common options

Run `gdiff <command> --help` for the full option list.

| Command | Option | Default | Purpose |
|---------|--------|---------|---------|
| `sketch` | `-k`, `-w`, `-h` | `23`, `k`, derived from `k` | k-mer length, minimizer window length, and LSH positions |
| `sketch` | `--frac` | `0.5` | Fraction of LSH buckets retained |
| `sketch` | `-l` | `333` | Sampled window length in k-mers; `0` omits windows |
| `sketch` | `--sample-size` | `1000` | Maximum windows stored per genome |
| `map` | `--levels` | `0.05 0.01 0.005` | One to three two-sided background levels |
| `map` | `--sample-size` | `1000` | Background windows per reference; at least `100` |
| `map` | `-b`, `--bin-shift` | `0` | Group k-mers into bins of size `2^b` to reduce memory and coordinate resolution |
| `map`, `roll` | `--hdist-th` | `3` | Maximum Hamming distance for a k-mer hit |
| `dist` | `--hdist-th` | `3`, or `2` for inputs >20 Mb | Maximum Hamming distance for a k-mer hit |
| Global | `--num-threads` | `1` | Worker threads |
| Global | `--seed` | `0` | Random seed |
| Global | `--verbose` | off | Report progress when stderr is not a terminal |

The number of sampled windows is capped by the available starts. Sequences shorter than `L + k - 1` bases contribute no windows; check `gdiff info` for the count stored in a sketch.

## Visualization

`gviz.py` creates HTML reports with filters for distance profiles (`roll`), pairwise summaries (`dist`), and sampled distances (`dist --output-samples`).

```bash
pip install numpy pandas
python gviz.py windows.tsv --seq chr7 --open
python gviz.py distances.tsv --open
```

Use `--annotations genes.gff3` to add a gene track and `--bins` to change profile resolution. Reports load Plotly.js from a CDN; for offline use, install `plotly` and pass `--inline`.

`plot.py` provides a tree-based explorer for legacy interval tables (`QUERY_ID`, `REF_ID`, `INTERVAL_START`, etc.):

```bash
pip install -r plot-reqs.txt
python plot.py --input legacy_intervals.tsv --tree tree.nwk --query query_name
```

Analysis scripts and datasets from the paper are available in [gdiff-analyses](https://github.com/bo1929/gdiff-analyses).
