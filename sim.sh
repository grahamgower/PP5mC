#!/bin/sh
# Simulate datasets to compare PP5mC with HBS-tools.

export odir=simout
mkdir -p $odir

export PATH=$HOME/src/pp5mC:$HOME/src/pp5mC/HBS-tools/hbs_tools.v1.1:$HOME/src/bowtie-1.2.2-linux-x86_64:$PATH
export ref=/localscratch/Refs/H_sapiens/chr1-100mb/chr1-100mb.fa
nseqs=100000
readlen=150
profile1=profile.INS.r1.txt
profile2=profile.INS.r2.txt

fraglens="4.0|0.25   5.0|0.4"
moltypes="hp bs pal"

for fraglen in $fraglens; do

	# ancient DNA fraglen dist (median ~55)
	# R> x=20:100; y=dlnorm(x,4.0,0.25); plot(x,y)
	# modern DNA fraglen dist (median ~150)
	# R> x=20:350; y=dlnorm(x,5.0,0.4); plot(x,y)
	mu=$(echo $fraglen | cut -f1 -d\|)
	sigma=$(echo $fraglen | cut -f2 -d\|)

	for moltype in $moltypes; do
		pfx="sim-r${readlen}-u${mu}-s${sigma}-${moltype}"

		echo simhbs \
			-m $moltype \
			-u $mu \
			-s $sigma \
			-E $profile1,$profile2 \
			-n $nseqs \
			-o $odir/$pfx \
			$ref

	done
done | parallel

# Concatenate $pfx's into $sims. Can't do this above, as the for loop runs
# in a separate process when piping, so $sims loses its value.
sims=""
for fraglen in $fraglens; do
	mu=$(echo $fraglen | cut -f1 -d\|)
	sigma=$(echo $fraglen | cut -f2 -d\|)
	for moltype in $moltypes; do
		pfx="sim-r${readlen}-u${mu}-s${sigma}-${moltype}"
		sims="$sims $pfx"
	done
done

do_hbs() {
	sim=$1

	/usr/bin/time \
	hbs_process \
		--hpmeth \
		${odir}/${sim}.r1.fq \
		${odir}/${sim}.r2.fq \
		adapters.fa \
		hairpins.fa \
		> ${odir}/${sim}.hbs_process.txt 2> ${odir}/${sim}.hbs_process-stderr.txt

	/usr/bin/time \
	hbs_mapper \
		--keep-temp 1 \
		-q ${ref%.fa} \
		${odir}/${sim}.r1.fq.processed \
		${odir}/${sim}.r2.fq.processed \
		> ${odir}/${sim}.hbs_mapper.txt 2> ${odir}/${sim}.hbs_mapper-stderr.txt

	# Don't really need this; can get stats from hbs_mapper's *.report.txt
	#/usr/bin/time \
	#hbs_methylation_extractor \
	#	--merge_non_CpG \
	#	${odir}/${sim}.r1.fq.processed.recovered.sam \
	#	> ${odir}/${sim}.hbs_methylation_extractor.txt 2> ${odir}/${sim}.hbs_methylation_extractor-stderr.txt
	#mv CpG_context_${sim}.r1.fq.processed.recovered.sam.txt ${odir}/
	#mv Non_CpG_context_${sim}.r1.fq.processed.recovered.sam.txt ${odir}/

	# check for empty file
	[ -s ${odir}/${sim}.r1.fq.processed.ori ] || return

	if /bin/false; then
		fq2html.py \
			--latex \
			--ori ${odir}/${sim}.r1.fq.processed.ori \
			${odir}/${sim}.r1.fq \
			${odir}/${sim}.r2.fq \
			> ${odir}/${sim}.HBS.tex
		pdflatex \
			--output-directory ${odir} \
			${odir}/${sim}.HBS.tex
	fi

	# map with bwa-mem, to match pp5mC
	bwa mem \
		-t 2 \
		$ref \
		${odir}/${sim}.r1.fq.processed.ori \
		> ${odir}/${sim}.hbs.bwa.sam
}
export -f do_hbs

do_pp5mC() {
	sim=$1

	/usr/bin/time \
	foldreads \
		-1 ${odir}/${sim}.r1.fq \
		-2 ${odir}/${sim}.r2.fq \
		-m ${odir}/${sim}.metrics \
		-u ${odir}/${sim}.unfolded \
		> ${odir}/${sim}.folded.fq 2> ${odir}/${sim}.foldreads-stderr.txt

	# check for empty file
	[ -s ${odir}/${sim}.folded.fq ] || return

	if /bin/false; then
		fq2html.py \
			--latex \
			${odir}/${sim}.folded.fq \
			> ${odir}/${sim}.folded.tex
		pdflatex \
			--output-directory ${odir} \
			${odir}/${sim}.folded.tex

		fq2html.py \
			--latex \
			${odir}/${sim}.unfolded_r1.fq \
			${odir}/${sim}.unfolded_r2.fq \
			> ${odir}/${sim}.unfolded.tex
		pdflatex \
			--output-directory ${odir} \
			${odir}/${sim}.unfolded.tex
	fi

	bwa mem \
		-C \
		-t 2 \
		$ref \
		${odir}/${sim}.folded.fq \
	| samtools view \
		-Sbu - \
	| samtools sort \
		-O bam \
		-@ 2 \
		-m 1G \
		-T tmp.sort.$sim \
		- \
	> ${odir}/${sim}.folded.bwa.bam

	samtools index ${odir}/${sim}.folded.bwa.bam
	samtools view ${odir}/${sim}.folded.bwa.bam > ${odir}/${sim}.folded.bwa.sam

	# Map using bowtie1, with the same parameters used by HBStools.
	# Actually, hbs_mapper uses `-k 2', and then excludes reads for
	# which two equally good alignments exist.  In bwa, MAPQ filtering
	# achieves this goal, but bowtie does not set MAPQ appropriately.
	# So we use `-k 1' here (single best alignment reported only), and
	# avoid MAPQ filtering in bwa for comparison purposes.
	bowtie \
		-q -n 2 -l 28 -k 1 --best --chunkmbs 512 \
		--sam --no-unal \
		-p 2 \
		${ref%.fa} \
		${odir}/${sim}.folded.fq \
		> ${odir}/${sim}.folded.bt1.sam

	# check for empty file
	[ -s ${odir}/${sim}.folded.bwa.sam ] || return

	/usr/bin/time \
	mark5mC \
		${odir}/${sim}.folded.bwa.bam \
		$ref \
		> ${odir}/${sim}.folded.bwa.methlist.txt 2> ${odir}/${sim}.mark5mC-stderr.txt

	scanbp \
		${odir}/${sim}.folded.bwa.bam \
		> ${odir}/${sim}.folded.pairs.txt

	plot_nt_pairing.py \
		--title $sim \
		${odir}/${sim}.folded.pairs.txt \
		${odir}/${sim}.pairs.pdf
}
export -f do_pp5mC


for sim in $sims; do
	echo do_hbs $sim
	echo do_pp5mC $sim
done | parallel -j 8


# Print stats for comparison.
for sim in $sims; do
	# molecules recovered
	mrhbs=$(awk '/Recovered .* original sequences/ {print 100*$2/'$nseqs'}' ${odir}/${sim}.hbs_mapper-stderr.txt)
	mrfold=$(awk '/^NF/ {print 100*$2/'$nseqs'}' ${odir}/${sim}.metrics)

	# correctly mapped (BWA)
	if [ -s ${odir}/${sim}.hbs.bwa.sam ]; then
		cm1hbs=$(checkseq2.py ${odir}/${sim}.r1.fq ${odir}/${sim}.hbs.bwa.sam | awk '{print 100*$3/'$nseqs'}')
	else
		cm1hbs=0
	fi
	if [ -s ${odir}/${sim}.folded.bwa.sam ]; then
		cm1fold=$(checkseq2.py ${odir}/${sim}.r1.fq ${odir}/${sim}.folded.bwa.sam | awk '{print 100*$3/'$nseqs'}')
	else
		cm1fold=0
	fi

	# correctly mapped (Bowtie1)
	if [ -s ${odir}/${sim}.r1.fq.processed.recovered.sam ]; then
		cm2hbs=$(checkseq2.py ${odir}/${sim}.r1.fq ${odir}/${sim}.r1.fq.processed.recovered.sam | awk '{print 100*$3/'$nseqs'}')
	else
		cm2hbs=0
	fi
	if [ -s ${odir}/${sim}.folded.bt1.sam ]; then
		cm2fold=$(checkseq2.py ${odir}/${sim}.r1.fq ${odir}/${sim}.folded.bt1.sam | awk '{print 100*$3/'$nseqs'}')
	else
		cm2fold=0
	fi

	printf "%s\t& %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\\\\ \n" \
		$sim $mrhbs $mrfold $cm1hbs $cm1fold $cm2hbs $cm2fold
done
echo
for sim in $sims; do
	# correct bases
	if [ -s ${odir}/${sim}.r1.fq.processed.ori ]; then
		cbhbs=$(checkseq.py ${odir}/${sim}.r1.fq ${odir}/${sim}.r1.fq.processed.ori | awk '{print 100*$3}')
	else
		cbhbs=0
	fi
	if [ -s ${odir}/${sim}.folded.fq ]; then
		cbfold=$(checkseq.py ${odir}/${sim}.r1.fq ${odir}/${sim}.folded.fq | awk '{print 100*$3}')
	else
		cbfold=0
	fi

	# error rates for mC vs C
	# (all CpG's were methylated to start with, all other C's were unmethylated)
	if [ -s ${odir}/${sim}.folded.bwa.methlist.txt ]; then
		fpfold=$(awk '	$8~/^CG/     {cg_C+=$6; cg_mC+=$7}
				$8~/^C[CAT]/ {ch_C+=$6; ch_mC+=$7}
				END {
					if (cg_mC) a=cg_C/cg_mC ; else a=0;
					if (ch_C) b=ch_mC/ch_C ; else b=0;
					print a, b
				}' ${odir}/${sim}.folded.bwa.methlist.txt)
		fpfold1=$(echo $fpfold | awk '{print 100*$1}')
		fpfold2=$(echo $fpfold | awk '{print 100*$2}')
	else
		fpfold1=0
		fpfold2=0
	fi
	if [ -s ${odir}/${sim}.r1.fq.processed.recovered.sam.report.txt ]; then
		fphbs=$(awk '   / methylated CpG/ {cg_mC=$6}
				/ methylated CH/ {ch_mC+=$6}
				/ unmethylated CpG/ {cg_C=$6}
				/ unmethylated CH/ {ch_C+=$6}
				END {
					if (cg_mC) a=cg_C/cg_mC ; else a=0;
					if (ch_C) b=ch_mC/ch_C ; else b=0;
					print a, b
				}' ${odir}/${sim}.r1.fq.processed.recovered.sam.report.txt)
		fphbs1=$(echo $fphbs | awk '{print 100*$1}')
		fphbs2=$(echo $fphbs | awk '{print 100*$2}')
	else
		fphbs1=0
		fphbs2=0
	fi

	printf "%s\t& %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\\\\ \n" \
		$sim $cbhbs $cbfold $fphbs1 $fpfold1 $fphbs2 $fpfold2
done

