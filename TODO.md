* Refactor `gamma.cpp/hpp` at some point.
* Revisit output formats and CLI.
* Have a map mode where you can specify the direction, too.
* What about Ns and repeats?
* `sdust` is not wired up; when it is, fix the exclusion zone (`rqseq.cpp`: `(i > mrs && i - k < mre)`)
  and reset the minimizer window on mask exit, otherwise stale pre-mask k-mers leak through.
## Outstanding from the GDSK/symmetric rewrite

Memory items that were planned but not done:
* `HDHist` counters are still 64-bit; they are prefix sums bounded by the query k-mer count, so
  32-bit would halve that allocation once `extract_histogram` is reworked off its u64 SIMD loads.
* `map` has neither `--low-memory` nor `--batch-bases`; only `detect` got them, so a large query
  through `map` still allocates the extrema arrays and loads the whole FASTA.
* `dist` in FASTA mode still loads the whole query rather than doing two gzip passes and keeping
  only the sampled windows.
* No peak-DIM-bytes estimate and no `-b` recommendation warning.
* `DIM::sample_random_intervals` still does a full `iota` + `shuffle` (dim.cpp:574) instead of
  sparse rejection sampling; `skip_v` is a byte per bin rather than a bitset; `Detector` still
  double-buffers its output; `extract_mers` does not reserve its packed key buffer.

Test and doc gaps:
* No `test_hll.cpp`, and the symmetric merge is pinned by hand-written expectations in
  `test/unit/test_sym.cpp` rather than golden vectors from the compiled C reference.
* The mx-vs-sx equivalence and the Pool-vs-Seq estimate equivalence were verified by hand but are
  not automated; only the Pool/Seq window geometry check is in the suite.
* `README.md` and `docs/subcommands.md` do not document `--window-repr`, `--symmetric`, `--lr-th`,
  `--min-portion`, `--low-memory`, or `--batch-bases`.

* `sym_est_t::max_distance` is the max over rows that survived the `lr_ub` filter, so when the
  portion rule falls back to the unfiltered mean it describes a set that had no bearing on the
  reported `distance` -- and can print smaller than it. Worth renaming or documenting in the output.

Behaviour worth knowing:
* `dist` peak RSS tracks the total container bytes it touches (about 11 MB per genome at defaults),
  because Darwin ignores `MADV_DONTNEED` for file-backed pages. The pages are clean and reclaimable,
  so this degrades to re-reading rather than OOM, and on Linux the hint should actually bound it.
* `hdist_th` should be clamped to k-h after a sketch load, with a warning.
* Record/table loads should reject headers whose counts exceed the remaining file bytes.

Reconsider:
- `detect.cpp`
- `detect.hpp `
- `dim.cpp`
- `dim.hpp`
- `dist.cpp`
- `dist.hpp`
- `gamma.cpp`
- `gamma.hpp`
- `gdiff.cpp` (?)
- `gdiff.hpp` (?)
- `map.cpp` (?)
- `map.hpp` (?)
- `sketch.cpp` (?)
- `sketch.hpp` (?)
- `types.hpp` (?)
- `scan.hpp` (?)
- `stils.hpp` (?)
- `sym.cpp`
- `sym.hpp`
