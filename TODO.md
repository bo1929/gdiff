
// TODO: Further optimize, and polish.
// TODO: A better output format with a header? For both modes...
// TODO: If we will compute the distances, perhaps have an option to merge all overlapping intervals per distance.
// TODO: Merge aggressively as per above reasons.
// TODO: Something interesting with the second derivatives? Is the minimum second-derivative interval interesting?


Make Windusrf write tests.
Update the makefile for gamma and review.

Output format!!!!! Not very clean. No header.
The entire boundary issue is so confusing. Check it and make it make sense.
## TODOs



Define a metric for overlapping start measuring how they overlapping:
before the test merging
after the testing, different lower alpha

- subsampling?
- backtracking
- multiple-pass
- change alpha
Reverse complement

Krepp distance per segment

Do one pass take min p val, merge, check neighbors, repeat

Take contig size into consideration in plots

no sdust?

no r/m splitting, sketch single already?

Implement a model such that given a distance vector with in ANI estimate for a genome, estimate tree distance vector. A simple model perhaps, with missing data (this is important). Get local distance vectors, a compare against its ANI vector. Make sure that initialization etc. doesn't matter. So the vectors are constrained. Maybe something smarter. Incorporating the estimated embedding for the genome-wide distances.


**Importante**
! Parameter checks, for example especially min tau!
! Strands kept separate for fw/rc, but the test is against the lower distance.
! Inversion detection based on strand comparison ups and downs?


! Review gidiff.cpp/hpp and polish. messages and names.

## names
digg
digit
dimes
dint
distil
drip
gdiff
gdip
gind
glint
grid
grind
mint
paint
pind
ping
pint
plint
point
rind

Save positions of k-mers and query sketch to skecht for this youu need metadata (len). just do rho1*rho2?
Sketch to sketch?
Better write/read binary blocks, perhaps with buffers, large chunks at once is bad for IO

https://claude.ai/share/c6995cef-f98d-4d4b-aeb1-352587ac1543

Find conditions to not even check intervals (save some of the prefmax/suffmin ?)
lots of intervals use a lot of memory, keep it temp, at least across threshold clear stuff

? extract_intervals? correctness of the algorithm
? How to report intervals -- positions of k-mers?
? Remove m/r and splitting. Might make it much faster.
? How to treat Ns and skip positions
? Be careful about reverse complement coordinates!
? Postprocessing or be careful for +/- rc/fw in downstream analysis


grep
src/map.cpp:      l = 0; // TODO: What to do for missing ones?
src/map.cpp:      continue; // TODO: How to propagate the missing ones?
src/map.cpp:{ // TODO: Revisit?
src/map.hpp:  // void optimize_loglikelihood(); // TODO: Is ever needed? If so, correct?
src/map.hpp:  // void skip_mer(uint64_t i); // TODO: Anything better than ignoring?
src/sketch.cpp:  vec_enc_it ix1 = sfhm->bucket_iter_start(offset); // TODO: Use pointers?


Also for krepp, see if memory accesses are bug-free
Clang-Tidy
Clang's Static Analyzer (scan-build)
AddressSanitizer and ThreadSanitizer

plot.py
=====
learn more about annotations? Colors?
Maybe don't ever zoom in on intevarls but show seletion also enable scrolling back mouse motions?
make the input compatbile and header aware.
====


=====
### Correctness of the algorithm:
For an interval (a, b) to be maximal (not properly contained in any larger negative-sum interval):

Left endpoint must be a running max: If ∃ a' < a with fdps_v[a'] >= fdps_v[a], then (a', b) has sum fdps_v[b] - fdps_v[a'] <= fdps_v[b] - fdps_v[a] < 0. So (a', b) also has negative sum and properly contains (a, b) → not maximal.
Right endpoint must be the farthest: If ∃ b' > b with (a, b') negative, then (a, b') properly contains (a, b) → not maximal.
Algorithm Correctness
The algorithm does exactly this:

Line 96-98: Only considers a where fdps_v[a] > fdpmax_v[a-1] (running maxima)
Lines 100-102: Finds rightmost b_star where (a, b_star) could have negative sum (using suffix minima)
Line 109-111: Verifies (a, b_star) actually has negative sum
Line 113-115: b_prev tracking ensures each b_star is claimed by the leftmost valid a
Line 117: Left-maximality check confirms no (a', b_star) with a' < a has negative sum
Linear Time
Outer loop: O(n)
b_curr only increases (never resets), total increments ≤ n−1
Total: O(n)
Verification with Counterexample
My earlier example fdc_v = [1, -3, 2, -1]:

(3, 4) has sum = −1, but (1, 4) also has negative sum and contains (3, 4)
So (3, 4) is NOT maximal under your definition
The algorithm correctly outputs only (1, 4) ✓
Conclusion
The algorithm is correct for finding all maximal (not contained in larger) negative-sum intervals in O(n) time. The comments are somewhat misleading (e.g., "a is a prefix maxima" should say "skip if a is NOT a new record high"), but the logic is sound.
#### Completeness
Yes, it finds **all** maximal intervals. Here's the proof:

### Completeness Proof

Let `(a*, b*)` be any maximal negative-sum interval. We show the algorithm outputs it.

**Step 1**: `a*` must be a running maximum.
- If not, `∃ a' < a*` with `fdps_v[a'] ≥ fdps_v[a*]`
- Then `(a', b*)` has sum ≤ `(a*, b*)` sum < 0
- So `(a', b*)` contains `(a*, b*)` with negative sum → contradicts maximality

**Step 2**: `b*` must equal `b_star` (the rightmost found by the algorithm for `a*`).
- `b_star` = rightmost `b` where `(a*, b)` has negative sum
- If `b* < b_star`, then `(a*, b_star)` contains `(a*, b*)` with negative sum → contradicts maximality
- So `b* = b_star`

**Step 3**: The `b_prev` check passes.
- If `b_prev == b*`, some earlier `a < a*` claimed `b*`
- That means `(a, b*)` was output, which contains `(a*, b*)` with negative sum → contradicts maximality
- So `b_prev ≠ b*`

**Step 4**: The left-maximality check passes.
- If `fdps_v[b*] < fdpmax_v[a*-1]`, then `∃ a' < a*` with `fdps_v[a'] > fdps_v[b*]`
- So `(a', b*)` has negative sum and contains `(a*, b*)` → contradicts maximality
- So the check passes

**Conclusion**: All checks pass → `(a*, b*)` is output. ∎

---
The algorithm is **complete**: every maximal interval is found, none are missed.
======
