#!/bin/sh
# Align folded (single end) reads

source ./fold.settings
module load BWA/0.7.15-foss-2016b SAMtools/1.4.1-foss-2016b

map() {
	pfx=$1
	ref=$2
	refname=$3

	fq=${pfx}.folded.fq.gz
	sm=$(echo $pfx | cut -d_ -f1)
	lb=$(echo $pfx | cut -d_ -f2)
	id=${pfx}
	bam_out=${pfx}.${refname}.bam

	if [ ! -f ${fq} ]; then
		die "${pfx}: cannot find ${fq}"
	fi

	bwa mem \
		-t $bwa_threads_per_sample \
		-R "@RG\tID:${id}\tSM:${sm}\tLB:${lb}\tPL:ILLUMINA" \
		-C \
		$ref \
		$fq \
	 | samtools view \
		-q 25 \
		-Sbu \
		- \
	 | samtools sort \
	 	-O bam \
	 	-m 1G \
		-@ $sort_threads_per_sample \
		-T tmp.sort.${pfx}.$refname \
		- \
		> ${bam_out} \
	|| die "${pfx}: bwa mem|samtools view|samtools sort"

	samtools index ${bam_out} || die "${pfx}: samtools index"
}

for s in $sample_lib_pairs; do
	for ref in $refs; do
		refname=$(basename $ref)
		refname=${refname%.fasta}
		map $s $ref $refname &
	done
done
wait
