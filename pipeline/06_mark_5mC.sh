#!/bin/sh

source ./fold.settings
module load grg-utils

methlist() {
	bam=$1
	pfx=$2
	ref=$3

	mark_5mC \
		-p $hairpin \
		-5 10 \
		-3 10 \
		${bam} \
		${ref} \
		| pigz -p 2 -c \
		> ${pfx}.methlist.txt.gz \
	|| die "${pfx}: mark_5mC"

	frobmethlist.py \
		--all \
		--gzip \
		${pfx}.methlist.txt.gz \
		${pfx} \
	|| die "${pfx}: frobmethlist.py"
}

for sample in $samples; do
	for ref in $refs; do
		refname=$(basename $ref)
                refname=${refname%.fasta}
		pfx=${sample}.${refname}
		bam=${pfx}.dedup.realigned.calmd.bam
		methlist $bam $pfx $ref &
	done
done
wait
