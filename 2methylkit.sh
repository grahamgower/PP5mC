#!/bin/sh

# https://github.com/al2na/methylKit

echo -e "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT"

awk 'BEGIN {fr["+"]="F"; fr["-"]="R"}
	NR>1 {printf "%s.%s\t%s\t%d\t%s\t%d\t%.2f\t%.2f\n",
		$1, $3, $1, $3, fr[$4], $6+$7, 100*$7/($6+$7), 100*$6/($6+$7)}'
