#!/bin/sh


conversion() {
	f=$1
	sample=`basename $f | cut -d . -f1`
	ref=`basename $f | cut -d . -f2`
	(
		printf "$sample\t$ref"
		unpigz -c $f | \
		awk 'NR>1 {
				if ($8~/CG[ACGTN]/) {
					cpg_C += $6
					cpg_mC += $7
				} else if ($8~/C[ACT][ACT]/) {
					chh_C += $6
					chh_mC += $7
				} else if ($8~/C[ACT]G/) {
					chg_C += $6
					chg_mC += $7
				}
			} END {
				C = cpg_C + chh_C + chg_C
				mC = cpg_mC + chh_mC + chg_mC
				conv = 100*C/(C+mC)
				printf "\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n",
					cpg_C, chh_C, chg_C, cpg_mC, chh_mC, chg_mC, conv
			}'
	) >> conversionStats.$sample.$ref.txt
}

for i in *.REFbisonWMG.methlist.txt.gz; do
	conversion $i &
done
wait

printf "SampleID\t\
Reference\t\
CpG\t\
CHH\t\
CHG\t\
mCpG\t\
mCHH\t\
mCHG\t\
ConversionRate\t\
\n" > conversionStats.txt
cat conversionStats.*.*.txt >> conversionStats.txt
rm -f conversionStats.*.*.txt
