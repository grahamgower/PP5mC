#!/bin/sh

# https://github.com/dpryan79/PileOMeth
# 1) The chromosome/contig/scaffold name
# 2) The start coordinate
# 3) The end coordinate
# 4) The methylation percentage rounded to an integer
# 5) The number of alignments/pairs reporting methylated bases
# 6) The number of alignments/pairs reporting unmethylated bases

echo 'track type="bedGraph" description="proportion of methylated Cytosines"'

awk 'NR>1 {printf "%s\t%d\t%d\t%d\t%d\t%d\n",
		$1, $2, $3, 100*$7/($6+$7), $7, $6}'
