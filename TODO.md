# TODO

Actionable items from the current review of `map.cpp`/`map.hpp`, `gamma.hpp`, `gdiff.cpp`/`gdiff.hpp`, and tests.

## Critical Bug Fixes

- **Make sampled/skipped query k-mers consistent between interval detection and distance estimation.** `search_mers()` skips `partial_offset() == max` k-mers, while `extract_histogram()` later counts all non-hit k-mers in the bin span as misses. Decide whether unsampled query k-mers are unobserved or misses, then make the likelihood contributions and histogram counts agree.

## P2 / Next Fixes

- **Make `MapSC::map()` truly batch queries or intentionally load all queries.** `QSeq::read_next_batch()` appends and never clears, so `map()` currently loads the whole query file into memory before threading. Either process sketches per batch, or rename/document this as whole-query loading and simplify the counter.
- **Seed per-thread RNGs deterministically.** `gen` is `thread_local`; `--seed` only seeds the main thread, so worker threads use default-identical streams. Derive worker seeds from `seed` and worker index/sketch index.
- **Audit strand/reference output.** Current output reports both `+` and `-` segments while genome-wide mesh/null sampling uses only the lower-distance strand. Clarify whether downstream users should compare strands, collapse strands, or treat the non-reference strand as inversion evidence.
- **Add regression coverage for multi-contig edge cases.** Include at least two query contigs with different bin counts, one intact/full-sequence segment, one no-hit segment, and one reverse-complement-selected reference strand.
- **Add end-to-end boundary tests for interval coordinates.** Unit helpers now cover conversion math; add CLI/QIE tests for `len == k`, `nbins == 2`, `tau_bin > nbins`, full-span `[1, nbins]`, and last-segment extension.

## Recently Addressed

- **Reference-strand significance:** `record_t` now stores whether a segment is on the null-sampled strand; that strand gets the two-sided p-value.
- **Intact segment significance skip:** records now use record-local site coordinates for full-sequence detection.
- **Unsupported `--hdist-th`:** CLI validation now rejects values above the fixed 8-lane histogram capacity.
- **Null-window overlap coordinate conversion:** null sampling converts record bin coordinates to 0-based half-open bin-boundary coordinates before overlap rejection.
- **Coordinate conventions:** coordinate conversion helpers now define enum-only interval coverage and continuous-mode boundary rows; README documents the public output convention.
- **Direct unit coverage:** added direct tests for reference-strand p-value adjustment, half-open null-window overlap boundaries, and bin-to-output coordinate conversion.

## Optimizations

- **Avoid full-file query buffering for large inputs.** The current map path scales memory with all query sequences plus all per-sketch result buffers. A batch-oriented design would also let continuous significance state be released earlier.
- **Reduce per-record null sampling cost.** `test_significance()` fits a gamma model for every segment and runs MH draws for every accepted mesh sample. Cache reusable null samples by query/reference context where statistically valid, and benchmark fewer MH draws per mesh with calibration tests.
- **Use a faster histogram path for common `hdist_th <= 4`.** The AVX-512/SIMDe path always handles 8 lanes. A small scalar loop or specialized 5-lane path may be faster on Apple Silicon where SIMDe emulates AVX-512.
- **Release interval and prefix buffers after use.** `release_accumulators()` exists but is disabled; revisit lifetime now that continuous mode saves only mesh histograms and records.
- **Benchmark `extract_intervals_mx()` vs `extract_intervals_sx()`.** Tests show agreement; keep one implementation unless both have a measured use.

## Design Cleanup

- **Move significance configuration into CLI/config.** `test_significance(params.nsamples, 250)` hard-codes `ntries`; expose or derive it from `nsamples` and expected rejection rate.
- **Replace magic constants with named config.** Examples: `hdist_bound = 7`, `npoints_min = 8`, `npoints_max = 128`, MH `S = 50`, MH `B = 100`, and distance upper bound `0.5`.
- **Clean dead/stale comments.** Remove commented-out code in `map.cpp` (`extract_intervals_sx`, `chisq_val`, `sample_box_muller` second draw, warnings) once decisions are made.
- **Document output schema.** README and tests should state coordinate indexing, interval end inclusivity, `MASK`, `D_INTERVAL`, `PERCENTILE`, and `FOLD`.
- **Clarify gamma model choice.** `GammaModel::fit()` supports noise-aware MLE/QMATCH, but `score_gamma()` uses `fit_from_samples()`. Keep both only if both are used or planned; otherwise remove the unused path or add tests for the intended public method.

## Validation

- **Run sanitizers on map/significance paths.** Use ASan/UBSan with tests that set `--hdist-th` near the supported maximum and include continuous-mode significance.
- **Run a deterministic multi-thread regression.** Same seed and inputs should produce identical output for `--num-threads 1`, `2`, and higher counts, modulo row ordering if explicitly allowed.
- **Refresh golden regression outputs only after coordinate/significance semantics are fixed.**
- **Compare segment distances against an independent implementation.** Use krepp or a small scalar reference for per-segment MLE and Fisher information.

## Lower Priority / Research

- **Apply FDR correction to segment p-values.**
- **Evaluate whether sdust filtering is still needed.**
- **Evaluate removing minimizer/remainder splitting if benchmarks show it is slower than direct sketching.**
- **Explore inversion detection from strand-specific high/low distance patterns.**
- **Investigate distance-vector models for tree-distance estimation.**

## Plotting

- **Make `plot.py` header-aware for the current output schema.**
- **Use contig length when scaling interval plots.**
- **Improve navigation/selection without zooming away from context.**
- **Add color handling for NaN/no-match and bounded distance intervals.**
