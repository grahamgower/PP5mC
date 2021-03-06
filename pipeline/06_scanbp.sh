#!/bin/sh

source ./fold.settings

scan() {
	bam=$1
	pfx=$2

	scanbp \
		${bam} \
		> ${pfx}.pairs.txt \
		|| die "${pfx}: scanbp"

	plot_nt_pairing.py \
		--title ${pfx} \
		${pfx}.pairs.txt \
		${pfx}.pairs.pdf \
		|| die "${pfx}: plot_nt_pairing.py"
}

for bam in *.realigned.calmd.bam; do
	pfx=${bam%.dedup.realigned.calmd.bam}
	scan $bam $pfx &
done
wait
