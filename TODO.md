distil
dimes
grind
grid
drip
dint
gdip
plint
glint
hind
gind
ping
pind
joint
mint
pint
hint
lint
paint
point
dint
rind
digit
digg

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
prune per sequence? Maybe? Optional?
Make sure that arrows look nice.
learn more about annotations? Colors?
Maybe don't ever zoom in on intevarls but show seletion also enable scrolling back mouse motions?
make compatible with this:  daler.github.io/gffutils/ffutils
make the input compatbile and header aware.
====
