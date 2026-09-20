# TODO

- Come up with a better output format for map subcommand -- for all paths and modes.
- Better tests and perhaps change the test data.
- Refactor and review:
  * modified:   src/map.cpp
  * modified:   src/map.hpp
  * modified:   src/records.hpp

## Map / detection

* Optional direction control for `map`: a strand-aware reference always scores both strands and
  reports the selected one; there is no way to restrict the search to a single strand.
* Ns and repeats: ambiguous k-mers break intervals (`IntExt::skip_mer`), but there is no
  low-complexity masking, and repeated sequence still produces duplicate hits.
* `hdist_th` should be clamped to `k - h` after a sketch load, with a warning, instead of being
  taken at face value from the command line.
* Record and table loads should reject headers whose counts exceed the remaining file bytes.

## Memory

* `Histogram::hist_v` / `miss_v` counters are 64-bit, but they are prefix sums bounded by the
  query k-mer count. Narrowing them to 32-bit would halve the largest per-sequence allocation
  once `extract_histogram` is reworked off its 64-bit SIMD loads.
* `map` has no low-memory or batched mode: a large query is loaded whole and the per-bin extrema
  arrays are kept for every bin.
* `skip_v` is one byte per bin, lazily allocated; a bitset would cut it 8x.
* `RSeq::extract_mers` does not reserve its packed key buffer, so the largest queries reallocate
  it repeatedly.
* No peak-`IntExt`-bytes estimate and no `-b` recommendation warning.

## Tests

* No unit test for `rqseq` / hyperloglog: `extract_mers` and the cardinality estimate are only
  exercised indirectly through `sketch`.
* The symmetric merge is pinned by hand-written expectations in `test/unit/test_sym.cpp` rather
  than golden vectors.
* The Pool-vs-Seq estimate equivalence was verified by hand but is not automated; only the
  window geometry check is in the suite.
* `dist` and `roll` have no golden files. They are covered by property tests (option
  equivalence, thread invariance, `-o` byte equality) instead of pinned values.

## Behaviour worth knowing

* `dist` peak RSS tracks the total container bytes it touches (about 11 MB per genome at
  defaults), because Darwin ignores `MADV_DONTNEED` for file-backed pages. The pages are clean
  and reclaimable, so this degrades to re-reading rather than OOM, and on Linux the hint should
  actually bound it.
* A container stores a wall-clock timestamp, so two runs over the same input with the same
  `--seed` are not byte-identical. Tests therefore compare content (sampled window starts,
  retained k-mers) rather than bytes. Making the stamp optional would give fully reproducible
  containers.

## Docs

* `docs/plan-*.md`, `docs/source-review.md` and `docs/easy-improvements-plan.md` describe designs
  that have since been retired (`detect`, `gamma`, `estimator`, `sym`). They are kept as history
  and are not maintained.
