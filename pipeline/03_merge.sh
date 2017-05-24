#!/bin/sh
# merge different libraries for the same sample

source ./fold.settings
module load SAMtools/1.4.1-foss-2016b

merge()
{
	sample=$1
	refname=$2
	bam_out=$3
	samtools merge \
		${bam_out} \
		${sample}_*.${refname}.bam \
	|| die "${sample}: samtools merge"
	samtools index ${bam_out} || die "${samtools}: samtools index"
}

for sample in $samples; do
	for ref in $refs; do
		refname=$(basename $ref)
		refname=${refname%.fasta}
		bam_out=${sample}.${refname}.bam
		n_bams=`ls ${sample}_*.${refname}.bam 2>/dev/null | wc -l`

		case $n_bams in
			0)
				die "${sample}: cannot find bam for ${refname}"
				;;
			1)
				# one library for this sample, so symlink files
				ln -s ${sample}_*.${refname}.bam ${bam_out}
				ln -s ${sample}_*.${refname}.bam.bai ${bam_out}.bai
				if [ ! -f ${bam_out}.bai ]; then
					rm -f ${bam_out}.bai
					die "${sample}: cannot find index for ${sample}_*.${refname}.bam"
				fi
				;;
			*)
				# multiple libraries
				merge $sample $refname $bam_out &
				;;
		esac
	done
done
wait
