Distil
Remove thread safety overhead.
Remove m/r and splitting. Might make it much faster.
be careful about reverse complement coordinates?
How to treat Ns? skip_mer
How to report intervals -- positions of k-mers?
extract_intervals? correctness of the algorithm
make lsh lookup cache efficient, perhaps keep a buffer
Find conditions to not even check intervals (save some of the prefmax/suffmin ?)
Save positions of k-mers and query sketch to skecht for this youu need metadata (len). just do rho1*rho2?
