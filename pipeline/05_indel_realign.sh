#!/bin/sh

source ./fold.settings
module load SAMtools/1.4.1-foss-2016b
module load Java/1.8.0_101 GATK/3.6-Java-1.8.0_101

# NOTE: -Xmx parameter is per thread
java="java -Xmx2g"
gdir=$EBROOTGATK
gatk="$java -jar $gdir/GenomeAnalysisTK.jar"

indel_realign() {
	pfx=$1
	pfx_out=${pfx}.realigned.calmd
	bam_in=${pfx}.bam

	if [ ! -f "$bam_in" ]; then
		die "${pfx}: cannot file ${bam_in}"
	fi

	# train realigner
	$gatk \
		-T RealignerTargetCreator \
		-nt $realigner_threads_per_sample \
		-R $ref \
		-I $bam_in \
		-o ${pfx}.intervals \
	|| die "${pfx}: RealignerTargetCreator"

	# do indel realignment
	$gatk \
		-T IndelRealigner \
		-R $ref \
		--bam_compression 0 \
		--disable_bam_indexing \
		-targetIntervals ${pfx}.intervals \
		-I $bam_in \
		-o ${pfx}.realigned.bam \
	|| die "${pfx}: IndelRealigner"

	# regenerate MD tag to match the indel realignment
	samtools calmd \
		-b \
		${pfx}.realigned.bam \
		$ref \
		> ${pfx}.realigned.calmd.bam \
	|| die "${pfx}: samtools calmd"

	samtools index \
		${pfx}.realigned.calmd.bam \
	|| die "${pfx}: samtools index"

	rm ${pfx}.intervals ${pfx}.realigned.bam
}

for sample in $samples; do
	for ref in $refs; do
		refname=$(basename $ref)
		refname=${refname%.fasta}
		indel_realign ${sample}.${refname}.dedup #&
	done
done
#wait
