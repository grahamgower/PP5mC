#!/bin/sh

source ./fold.settings
module load paleomix-meta
module load SAMtools/1.4.1-foss-2016b

dedup() {
	pfx=$1
	bam_in=${pfx}.bam
	bam_out=${pfx}.dedup.bam

	paleomix rmdup_collapsed \
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
