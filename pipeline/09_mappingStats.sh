#!/bin/sh

module load SAMtools/1.4.1-foss-2016b


stats() {
	f=$1
	samtools depth $f > ${f}.depth
	sample=`basename $f | cut -d . -f1`
	ref=`basename $f | cut -d . -f2`
	tot_pairs=`awk '/^Total read pairs/ {a+=$NF} END {print a}' ${sample}_*.metrics`
	folded=`awk '/folded read pairs/ {a+=$(NF-1)} END {print a}' ${sample}_*.metrics`
	complete_hp=`awk '/read pairs with complete/ {a+=$NF} END {print a}' ${sample}_*.metrics`
	mapped=`samtools view -c ${f/\.dedup\.realigned\.calmd/}`
	dedup=`samtools view -c $f`
	covered=`wc -l ${f}.depth | awk '{print $1}'`
	genome_size=`samtools idxstats $f | awk -v FS="\t" '{sum += $2} END {print sum}'`
	perc_covered=`echo "$covered/$genome_size*100" | bc -l | awk '{printf "%.2f\n",$1}'`
	depth=(`awk -v var="$genome_size" \
			'{count+=$3; countsq+=$3*$3} END \
			{printf("%.4f %.4f %.4f %.4f\n", \
			count/NR, \
			sqrt(countsq/NR - (count/NR)**2), \
			count/var, \
			sqrt(countsq/var - (count/var)**2))}' \
			${f}.depth`)
	length=(`samtools view $f \
			| awk 'BEGIN \
			{count=0; avg=0; std=0} \
			{count=count+1; lgth=lgth+length($10); std=std+(length($10)-lgth/count)*(length($10)-lgth/count)} END \
			{printf("%.2f %.2f\n", \
			lgth/count, \
			sqrt((std)/(count-1)))}'`)
	printf "$sample\t\
$ref\t\
$tot_pairs\t\
$folded\t\
$complete_hp\t\
$mapped\t\
$dedup\t\
$perc_covered\t\
${depth[0]}\t\
${depth[1]}\t\
${depth[2]}\t\
${depth[3]}\t\
${length[0]}\t\
${length[1]}\t\
\n" > mappingStats.$sample.$ref.txt
	rm ${f}.depth
}

for i in *.dedup.realigned.calmd.bam; do
	stats $i &
done
wait

printf "SampleID\t\
Reference\t\
Pairs\t\
Folded\t\
Hairpin\t\
Mapped\t\
Unique\t\
Coverage\t\
Depth_covered\t\
Depth_covered_Std\t\
Depth_genome\t\
Depth_genome_Std\t\
Length\t\
Length_Std\t\
\n" > mappingStats.txt

cat mappingStats.*.*.txt >> mappingStats.txt
rm -f mappingStats.*.*.txt
