# Demo

Run these commands from the repository root after `make`. Inputs are bundled in `demo/data/`:

| Genome | Role |
|--------|------|
| AG-899-N21 | Query |
| AG-894-M23 | Same-species reference |
| AG-893-E15 | Putative donor of the phnCDEubiE operon |

```bash
mkdir -p demo/output
```

## Sketch

Create reference and query sketches, including sampled windows for `dist`:

```bash
./gdiff sketch \
    -i demo/data/AG-894-M23_contigs.fasta.gz demo/data/AG-893-E15_contigs.fasta.gz \
    -o demo/output/references.gds

./gdiff sketch -i demo/data/AG-899-N21_contigs.fasta.gz \
    -o demo/output/query.gds
```

## Genome-wide distances: `dist`

```bash
./gdiff dist demo/output/query.gds demo/output/references.gds \
    -o demo/output/distances.tsv
```

Output: one row per query–reference pair. `d` is the genome-wide distance; smaller values mean greater similarity.

## Window distances: `roll`

```bash
./gdiff roll demo/data/AG-899-N21_contigs.fasta.gz demo/output/references.gds \
    -l 500 -s 500 -o demo/output/roll.tsv
```

Output: query coordinates and distance for each window against each reference. `-l` sets the window length in k-mers; `-s` sets the step in bases.

## Outlier intervals: `map`

```bash
./gdiff map demo/data/AG-899-N21_contigs.fasta.gz demo/output/references.gds \
    -l 500 --sample-size 20000 -o demo/output/intervals.tsv
```

Output: background-relative outlier intervals at levels `0.01 0.001 0.0001`, with normal merging. `-l` sets the minimum length in k-mers; 20,000 background samples improve tail estimation for this example.

Look near **49,276–52,521 bp** on `AG-899-N21_NODE_6`: distances to the donor drop, while those to the same-species reference rise or have no estimate. All coordinates are 1-based and inclusive.

To view the reported intervals:

```bash
python3 gviz.py demo/output/intervals.tsv --seq AG-899-N21_NODE_6
```

Each reference and outlier direction has its own track; gray intervals have no distance estimate. The donor has eight overlapping calls of roughly 22–24 kb on this 80 kb contig. These maximal intervals cover the operon but do not identify its exact boundaries.

Assemblies are from the GORG dataset used in [gdiff-analyses](https://github.com/bo1929/gdiff-analyses).
