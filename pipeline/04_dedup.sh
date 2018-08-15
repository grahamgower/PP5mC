#!/bin/sh

source ./fold.settings

dedup() {
	pfx=$1
	bam_in=${pfx}.bam
	bam_out=${pfx}.dedup.bam

	rmdup_collapsed.py \
		--remove-duplicates \
		< $bam_in > $bam_out \
	|| die "${pfx}: paleomix rmdup_collapsed"

	samtools index ${bam_out} || die "${pfx}: samtools index"
}

for sample in $samples; do
	for ref in $refs; do
		refname=$(basename $ref)
		refname=${refname%.fasta}
		dedup $sample.$refname &
	done
done
wait
