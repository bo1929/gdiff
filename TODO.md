# TODO

## Next TODOs

- Apply FDR for p-value correction.
- Make sure that sampled intervals do not overlap with the tested one(s).
- Make sure that examples are fine and correct.
- Make sure that the test is two-sided.

## Critical / Bugs

- **Parameter validation**: Add checks for edge cases, especially minimum tau.
- **Boundary handling in interval detection**: Audit and clarify boundary logic — currently confusing and potentially incorrect.
- **Reverse complement coordinates**: Fix/improve how reverse complement coordinates are reported.
- **Strand separation**: Strands are kept separate (fw/rc), but the test uses the lower distance — verify this is correct.
- **`extract_intervals` correctness**: Review the algorithm for correctness.

## High Priority

- **Output format**: Design a clean output format with a header for both modes.
- **Interval merging**: Add option to merge overlapping intervals per distance; consider aggressive merging.
- **Code review**: Review and polish `gdiff.cpp`/`gdiff.hpp` — improve messages, naming, and readability.
- **Write tests**: Add automated tests for core functionality.

## Open Questions

- **Remove sdust?** Evaluate whether sdust filtering is still needed.
- **Remove m/r splitting?** Sketching without minimizer/remainder splitting may be significantly faster.
- **How to handle Ns and skip positions?**
- **How to report interval positions?** Use positions of k-mers?
- **Inversion detection**: Can strand comparison (ups/downs of fw vs rc) be used to detect inversions?
- **Postprocessing for fw/rc**: What downstream handling is needed for +/- strand results?

## Improvements

- **Performance optimization**: Further optimize hot paths.
- **Second derivative analysis**: Investigate whether the minimum second-derivative interval is meaningful.
- **Overlapping metric**: Define a metric to measure interval overlap before/after merging.
- **One-pass merge strategy**: Do one pass, take min p-value, merge, check neighbors, repeat.
- **Memory usage for intervals**: Intervals use a lot of memory; keep temp data and clear across thresholds.
- **Early skip for intervals**: Find conditions to skip interval checking (e.g., save prefmax/suffmin).

## I/O and Sketching

- **Binary I/O**: Improve write/read of binary blocks — use buffers, avoid large single-chunk reads.
- **Save k-mer positions**: Store positions of k-mers and query sketch; requires metadata (lengths). Use `rho1*rho2`?
- **Sketch-to-sketch mode?**

## Static Analysis and Sanitizers

- Run Clang's Static Analyzer (`scan-build`).
- Run AddressSanitizer and ThreadSanitizer.
- Verify memory accesses in krepp (related project) are bug-free.
- Validate against krepp?

## plot.py

- **Header-aware input**: Make the plotting input compatible with the new output format.
- **Contig size in plots**: Take contig size into consideration.
- **Color scheme**: Experiment with non-linear color scales.
  - No match: 0.75+
  - One color range: 0.5–0.75
  - Add toggle for switching to a continuous color scale.
- **Selection/navigation**: Show selection without zooming into intervals; enable scroll-back and mouse interaction.
- **Annotations**: Improve annotation styling and colors. Maybe show the name when it's close.

## Future / Research Ideas

- **Distance vector model**: Given a per-segment distance vector and an ANI estimate for a genome, build a model to estimate a tree distance vector. Handle missing data. Constrain vectors. Possibly incorporate genome-wide distance embeddings.
- **krepp distance per segment**: Use krepp to compute/validate per-segment distances.
