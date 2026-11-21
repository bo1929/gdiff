# gdiff

`gdiff` estimates how sequence distance varies along a genome without alignment.
It compares a query sequence against a reference sketch and uses a pseudo-likelihood model to estimate distances and find unusually close or divergent regions.

Analysis scripts and datasets from the paper are available in [gdiff-analyses](https://github.com/bo1929/gdiff-analyses).

The same sketch supports three analyses:

- `dist`: randomly sample windows and estimate distances; summarize the sample into a genome-wide estimate.
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

The default output has one row per pair: the reconciled distance `d`; the unfiltered median, mean and highest window distance per direction (`d_median_ab`, `d_median_ba`, `d_mean_ab`, ...); and window counts.
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

`-l 500` sets the minimum interval length and background sampling window to 500 k-mers. gdiff samples query windows against each reference, fits a gamma background with a point mass at the distance floor, and derives lower and upper thresholds. The fit is fail-closed: if it is rejected, no threshold can be placed and the reference is skipped, so thresholds and `p_gamma` always come from the same model. Every skip is also written to the output as a tab-separated comment line, `# skipped<TAB>reference<TAB>reason`.

Each level `pi` is a one-sided tail probability: it places a lower threshold at the `pi` quantile and an upper threshold at the `1 - pi` quantile of the fitted background (no division by 2). The default levels are `0.01 0.001 0.0001`; smaller levels select more extreme tails. To use more permissive levels or sample more background windows:

```bash
gdiff map query.fasta reference.gds -l 500 --levels 0.05 0.01 0.005 --sample-size 1000 -o intervals.tsv
```

For the stringent default levels, `--sample-size 1000` provides more background support for the fitted tails. Model thresholds outside the valid distance range are dropped, as is a lower threshold that falls inside the point mass at the distance floor. The thresholds that were used, and any that were dropped with the reason, are listed before the records as tab-separated comment lines:

```
# threshold<TAB>reference<TAB>low|high<TAB>distance<TAB>pi
# dropped<TAB>reference<TAB>low|high<TAB>distance<TAB>pi<TAB>reason
```

To skip threshold detection entirely, give distance cutoffs directly:

```bash
gdiff map query.fasta reference.gds -l 500 --thresholds 0.02 0.3 -o intervals.tsv
```

Each value becomes one lane, used as-is (up to 8, mutually exclusive with `--levels`). The background is still sampled for `d_null` and the p-value columns; `pi` reports `nan` because no level placed the threshold. Values must be below 0.5; a lane at or below the distance floor, or above `d_ceiling`, is kept but warned about because it cannot (or can only rarely) be reached.

On each side, the most extreme level claims intervals first; later levels search the remaining gaps. Use `--independent-levels` to enumerate each level separately, allowing overlapping calls. A container with multiple references produces results for each reference.

### Inspect or merge sketches

```bash
gdiff info -i reference.gds
gdiff merge -i ref_A.gds ref_B.gds -o combined.gds
```

`info` reports sketch parameters and stored window counts. `merge` requires compatible sketch configurations; when starting from sequences, pass multiple inputs to `sketch` instead.
