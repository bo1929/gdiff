#!/bin/bash

rm -rf sketches && mkdir -p sketches
rm -rf est && mkdir -p est

SECONDS=0
export sketching_options='-k 27 -w 31 -h 11 -m 2 -r 1 --frac'
cat genome_names.txt | \
  xargs -I{} -P 16 bash -c \
  '../gdiff sketch `echo ${sketching_options}` -i genomes/{}.fna.gz -o sketches/{}.skc > /dev/null 2>&1'

export mapping_options="-d 0.10 -l 9900 --chisq 10000 --enum-only"
cat genome_pairs.txt | \
   xargs -n 2 -L 1 -P 16 bash -c \
   '../gdiff map `echo ${mapping_options}` -i sketches/"$1".skc -q genomes/"$0".fna.gz -o est/query_"$0"-ref_"$1".enum.txt > /dev/null 2>&1'

printf '%dh:%dm:%ds\n' $((SECONDS/3600)) $((SECONDS%3600/60)) $((SECONDS%60))

diff <(cd est && ls *.enum.txt | sort) <(cd gt && ls *.enum.txt | sort)
for f in est/*.enum.txt; do diff -q "$f" "gt/$(basename "$f")"; done
