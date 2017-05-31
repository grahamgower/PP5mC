#!/bin/sh

source ./fold.settings
module load grg-utils

fold() {
	pfx=$1
	fq_out=${pfx}.folded.fq.gz

	f1=`mktemp -u /tmp/foldreads.$pfx.fifo.XXXXXX`
	f2=`mktemp -u /tmp/foldreads.$pfx.fifo.XXXXXX`
	mkfifo $f1 $f2 || die "${pfx}: mkfifo $f1 $f2"

	if [ `ls $data/${pfx}*R1.fastq.gz 2>/dev/null | wc -l` = "0" ]; then
		die "${pfx}: cannot find $data/${pfx}*R1.fastq.gz"
	fi
	if [ `ls $data/${pfx}*R2.fastq.gz 2>/dev/null | wc -l` = "0" ]; then
		die "${pfx}: cannot find $data/${pfx}*R2.fastq.gz"
	fi

	unpigz -c $data/${pfx}*R1.fastq.gz > $f1 &
	unpigz -c $data/${pfx}*R2.fastq.gz > $f2 &

	foldreads \
		-p $hairpin \
		-m ${pfx}.metrics \
		-1 $f1 \
		-2 $f2 \
		| pigz -p 2 -c - > ${fq_out} \
	|| die "${pfx}: foldreads|pigz"

	rm $f1 $f2

	plot_fraglen_from_fasta.py \
		--title ${pfx} \
		${fq_out} \
		${pfx}.folded.length-histogram.pdf \
	|| die "${pfx}: plot_fraglen_from_fasta.py"
}

for s in $sample_lib_pairs; do
	fold $s &
done

wait
