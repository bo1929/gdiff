# TODO

Live backlog. The full code-review worklist (dead code, error paths, robustness, style batch,
prioritized order) is in `docs/REVIEW-20OCT26.md` — items are not duplicated here. Docs index:
`docs/README.md`; retired material: `docs/archive/`.

## Map / detection

* Decide the final `map` output format and either print or drop the per-record statistics that
  are computed but never emitted (Fisher info, `lr_lb`, `lr_ub`, empirical `pvalue`/`fold`;
  REVIEW M1). Two extra Brent fits per record are the single largest wasted cost in `map`.
* Optional direction control: a strand-aware reference always scores both strands and reports
  the selected one; there is no way to restrict the search to a single strand.
* Ns and repeats: ambiguous k-mers break intervals (`IntExt::skip_mer`), but there is no
  low-complexity masking, and repeated sequence still produces duplicate hits.
* `hdist_th` should be clamped to `k - h` after a sketch load, with a warning, instead of being
  taken at face value from the command line (the lr encoding can only express `k - h`
  positions).
* Container loads should reject headers whose counts exceed the remaining file bytes before
  resizing (`nsketches` for the index, `nwins` for the window list).

## Memory

* `Histogram::hist_v` / `miss_v` counters are 64-bit, but they are prefix sums bounded by the
  query k-mer count. Narrowing them to 32-bit would halve the largest per-sequence allocation
  once `extract_histogram` is reworked off its 64-bit SIMD loads.
* `map` has no low-memory or batched mode: the whole query is loaded (`QSeq` with
  `bpmax = uint64_max`) and the per-bin extrema arrays are kept for every bin.
* `skip_v` is one byte per bin, lazily allocated; a bitset would cut it 8x.
* `build_sketch`'s key-buffer reserve is sized from the first sequence only; a short first
  record in a large multi-FASTA causes repeated reallocations of `keys_v`.
* No peak-`IntExt`-bytes estimate and no `-b` recommendation warning.

## Tests / regression

* `test/test_regression.sh` and the `test/gt/` goldens are broken against the current CLI: the
  script invokes removed flags (`-d`, `--iterative-removal`) and the goldens predate the
  22/24-column `map` schema. Update the script's option sets (default vs `--independent-levels`
  modes), then delete `gt/*` and re-run to bootstrap new goldens.
* No unit tests for `rqseq` / hyperloglog: `extract_mers` (minimizer selection + dedup) and the
  cardinality estimate are only exercised indirectly through `sketch`.
* The symmetric merge is pinned by hand-written expectations in `test/unit/test_sym.cpp` rather
  than golden vectors.
* Pool-vs-packs (`--keep-seq`) estimate equivalence was verified by hand but is not automated;
  only the window geometry check is in the suite.
* `dist` and `roll` have no golden files. They are covered by property tests (option
  equivalence, thread invariance, `-o` byte equality) instead of pinned values.

## Tooling / paper

* `plot.py` and `gviz.py` still parse the pre-simple-map `map` schema (`MASK`, `D_INTERVAL`,
  `PERCENTILE`, ...). Resync them to the 22/24-column format, then resurrect the integration
  notes and the smoke test from `docs/archive/` (`plot-notes.md`, `test-gviz.sh`,
  `verify-gviz.js`).
* `docs/algorithm.tex` predates the current null model: it describes a plain moment-fit Gamma
  with percentiles, not the hurdle fit (spike at `d_min` + ML-refined Gamma), the empirical
  fallback, or the `p_interval`/`q_interval`/`p_gamma` columns. Update, then rebuild the pdf
  (`latexmk -pdf algorithm.tex`).

## Behaviour worth knowing

* `dist` peak RSS tracks the total container bytes it touches (about 11 MB per genome at
  defaults), because Darwin ignores `MADV_DONTNEED` for file-backed pages. The pages are clean
  and reclaimable, so this degrades to re-reading rather than OOM, and on Linux the hint should
  actually bound it.
* A container stores a wall-clock timestamp, so two runs over the same input with the same
  `--seed` are not byte-identical. Tests therefore compare content (sampled window starts,
  retained k-mers) rather than bytes. Making the stamp optional would give fully reproducible
  containers.

## Housekeeping

* `docs/` is gitignored — none of the documentation is versioned. Decide whether to track it
  (or at least `subcommands.md` / `architecture.md` / `style-guide.md`).
* Stray build artifact `.d` sits at the repo root (ignored); safe to delete.
