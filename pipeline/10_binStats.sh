#!/bin/sh

set -eu -o pipefail -o verbose

module load BEDTools/2.26.0-foss-2016b

bins=/localscratch/grg/fold-bovids/Cow_UMD3_1.1kb.bins.bed

bedgraph() {
	f=$1
	unpigz -c $f \
	| awk -v OFS="\t" '$1~/^chr/ {print $0,$5+$6}' \
	| sort -S 50G -k1,1 -k2,2n
}

for i in *.Cow_UMD3_1.pileOmeth.*.txt.gz; do
	echo $i
	bedgraph $i \
	| bedtools map \
		-a $bins \
		-b - \
		-c 4,5,5,5,5,5,5,7,7,7,7,7,7 \
		-o mean,min,max,mean,median,count,sum,min,max,mean,median,count,sum \
	| pigz \
	> ${i/\.txt\.gz/.binStats.gz}
done

