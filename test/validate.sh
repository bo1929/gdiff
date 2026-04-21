#!/bin/bash

rm -rf sketches && mkdir -p sketches
rm -rf est && mkdir -p est

export sketching_options='-k 27 -w 31 -h 11 -m 2 -r 1 --frac'
cat genome_names.txt | \
  xargs -I{} -P 16 bash -c \
  '../gdiff sketch `echo ${sketching_options}` -i genomes/{}.fna.gz -o sketches/{}.skc > /dev/null 2>&1'

for arg in "$@"; do
	# The segmentation mode with --enum-only
	if [[ "$arg" == "enum" ]]; then
		SECONDS=0
		export mapping_options="-d 0.10 -l 9900 --chisq 10000 --enum-only"
		cat genome_pairs.txt | \
			xargs -n 2 -L 1 -P 16 bash -c \
			'../gdiff map `echo ${mapping_options}` -i sketches/"$1".skc -q genomes/"$0".fna.gz -o est/query_"$0"-ref_"$1".enum.txt > /dev/null 2>&1'

		printf '%dh:%dm:%ds\n' $((SECONDS/3600)) $((SECONDS%3600/60)) $((SECONDS%60))
		diff <(cd est && ls *.enum.txt | sort) <(cd gt && ls *.enum.txt | sort)
		for f in est/*.enum.txt; do diff -q "$f" "gt/$(basename "$f")"; done
	fi

	# The default mode with MLE distances
	if [[ "$arg" == "cont" ]]; then
		SECONDS=0
		export mapping_options="-d 0.10 -l 9900 --chisq 33.00051"
		cat genome_pairs.txt | \
			xargs -n 2 -L 1 -P 16 bash -c \
			'../gdiff map `echo ${mapping_options}` -i sketches/"$1".skc -q genomes/"$0".fna.gz -o est/query_"$0"-ref_"$1".cont.txt > /dev/null 2>&1'

		printf '%dh:%dm:%ds\n' $((SECONDS/3600)) $((SECONDS%3600/60)) $((SECONDS%60))

		diff <(cd est && ls *.cont.txt | sort) <(cd gt && ls *.cont.txt | sort)
		# for f in est/*.cont.txt; do diff -q "$f" "gt/$(basename "$f")"; done
		for f in est/*.cont.txt; do 
      # cut -f1,2,3,4,5,7,8 "$f" | tail -n+2 > /tmp/est
      cut -f1,2,3,4,5,6,7 "$f" | tail -n+2 > /tmp/est
      cut -f1,2,3,4,5,6,7 "gt/$(basename "$f")" | tail -n+2 > /tmp/gt
      diff -q /tmp/est /tmp/gt
    done
	fi
done
