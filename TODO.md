# TODO
make sketch positional as well and merge -i

Live backlog. The open code-review findings (error paths, robustness, dead code, style batch,
prioritized order) are in `docs/REVIEW.md`; items are not duplicated here. Docs index:
`docs/README.md`; retired material: `docs/archive/`.

## Map / detection

* Revisit the default `--levels` and the scan inflation (decided 2026-11-20 to change nothing
  for now): the scan takes maxima over many windows, so the loosest levels can call spurious
  intervals on an unedited genome.
* `map` no longer computes the unprinted per-record statistics (Fisher info, `lr_lb`, `lr_ub`,
  empirical `pvalue`/`fold`). If identifiability diagnostics are wanted, add them as printed
  `p_lb`/`p_ub` columns (see `NOTES.md`, section 3).
* The hurdle fit has no point mass at the ceiling (`UB = 0.5 - 1e-7`): on very divergent
  genomes the null median can sit at the cap.
* Optional direction control: a strand-aware reference always scores both strands and reports
  the selected one; there is no way to restrict the search to a single strand.
* Ns and repeats: ambiguous k-mers break intervals (`IntExt::skip_mer`), but there is no
  low-complexity masking, and repeated sequence still produces duplicate hits.
* `hdist_th` should be clamped to `k - h` after a sketch load, with a warning, instead of being
  taken at face value from the command line (the lr encoding can only express `k - h`
  positions).
* Container loads should reject headers whose counts exceed the remaining file bytes before
  resizing (`nsketches` for the index, `nwins` for the window list).

## dist

* `dist` builds its likelihood without the reference's spurious-hit fractions
  (`make_llhf(reference, hdist_th)`), while `map` and `roll` pass `get_spurious_f(hdist_th)`.
  Decide whether that is intended; applying it would change `dist` output.

## roll

* Add `p`-values, FDR and the gamma fit to `roll` (see `NOTES.md`).

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

* No unit tests for `rqseq` / hyperloglog: `extract_mers` (minimizer selection + dedup) and the
  cardinality estimate are only exercised indirectly through `sketch`.
* The symmetric merge is pinned by hand-written expectations in `test/unit/test_sym.cpp` rather
  than golden vectors.
* Pool-vs-packs (`--keep-seq`) estimate equivalence was verified by hand but is not automated;
  only the window geometry check is in the suite.
* `dist` and `roll` have no golden files. They are covered by property tests (option
  equivalence, thread invariance, `-o` byte equality) instead of pinned values.

## Paper

* `docs/algorithm.tex` predates the current null model: it describes a plain moment-fit Gamma
  with percentiles, not the hurdle fit (point mass at `d_floor` + ML-refined Gamma), the
  `UB = 0.5 - 1e-7` domain, the `d_seq` side rule, or the `p_interval`/`q_interval`/`p_gamma`
  columns. Update, then rebuild the pdf (`latexmk -pdf algorithm.tex`).

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

* `docs/` is gitignored; none of the documentation is versioned. Decide whether to track it
  (or at least `subcommands.md` / `architecture.md` / `style-guide.md`).
