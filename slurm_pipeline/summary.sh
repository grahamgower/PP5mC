#!/bin/sh

json=$1

if [ -z "$1" ]; then
	echo "usage: $0 config.json"
	exit 1
fi

count_reads()
{
	awk '
		NR==1 {n=$1}
		$4~/supplementary/ {sup=$1}
       		END {print n-sup}
	' $1
}

summarise_metrics()
{
	awk '
		/^(NR)|(NF)/ { a[$1]+=$2 }
		/^LN/ { ln+=$2; nln+=1 }
		END { print a["NR"], a["NF"], ln, nln }
	'
}

cov_depth()
{
	stats=$1
	len=$2
	awk '
		/COV/ { dp+=$3*$4; cov+=$4 }
		END { print cov/'$len', dp/'$len', dp/cov }
	' $stats
}

jf=$(mktemp /tmp/jdata.XXXXXX)
jsonhelper.py $json > $jf
samples=$(cut -f1 $jf | sort -Vu)
refs=$(jsonhelper.py $json refs | cut -f1)

echo -e "sample\tlib\tnreads\tnfolded\tfoldratio\tavelen\tref\tmapped_raw\tmapped_dedup\tendog\tclon\tmt_cov\tmt_dp0\tmt_dp1\tchrX_cov\tchrX_dp0\tchrX_dp1\taut_cov\taut_dp0\taut_dp1"
for ref in $refs; do

	fai=$(jsonhelper.py $json refs | awk '/'$ref'/ {print $2}').fai
	mt_len=$(awk '$1~/^(MT|Mt|M)$/ {print $2; exit}' $fai)
	aut_len=$(awk '$1~/^([Cc]hr)?[1-9][0-9]?$/ {a+=$2} END {print a}' $fai)
	chrX_len=$(awk '$1~/^([Cc]hr)?X$/ {print $2; exit}' $fai)

for sample in $samples; do
	libs=$(awk '$1=="'$sample'" {print $2}' $jf | sort -Vu)

	s_nr=0
	s_nf=0
	s_l_nf=0
	s_ln=0
	s_nln=0
	s_mapped_raw=0
	s_mapped_dd=0

	nr=0; nf=0; l_nf=0; ln=0; nln=0;
	mapped_raw=0
	mapped_dd=0

	for lib in $libs; do
		pfx=$sample/$lib/${sample}_${lib}
		runs_nolo=$(awk '$1=="'$sample'" && $2=="'$lib'" && $3!~/^LO/ {print $3}' $jf | sort -Vu)
		runs_lo=$(awk '$1=="'$sample'" && $2=="'$lib'" && $3~/^LO/ {print $3}' $jf | sort -Vu)

		# exclude `LO' runs from foldratio/avelen stats
		set -- $(for r in $runs_nolo; do cat ${pfx}_${r}.metrics; done 2>/dev/null | summarise_metrics)
		nr=$1; nf=$2; ln=$3; nln=$4;

		if [ "$nr" == "0" ]; then
			continue
		fi

		l_nf=0; l_nr=0
		if [ ! -z "$runs_lo" ]; then
			# count the `LO' runs
			set -- $(for r in $runs_lo; do cat ${pfx}_${r}.metrics; done 2>/dev/null | summarise_metrics)
			l_nr=$1; l_nf=$2; l_ln=$3; l_nln=$4;
		fi

		#A3020-AncRAD8k/MCS20-3/A3020-AncRAD8k_MCS20-3.UMD311.bam.flagstat.txt
		mapped_raw=$(count_reads ${pfx}.${ref}.bam.flagstat.txt)
		mapped_dd=$(count_reads ${pfx}.${ref}.dedup.bam.flagstat.txt)

		# sum over all libs
		s_nr=$((s_nr+nr))
		s_nf=$((s_nf+nf))
		s_l_nf=$((s_l_nf+l_nf))
		s_ln=$(perl -e "print($s_ln+$ln)")
		s_nln=$((s_nln+nln))
		s_mapped_raw=$((s_mapped_raw+mapped_raw))
		s_mapped_dd=$((s_mapped_dd+mapped_dd))

		foldratio=$(perl -e "print ($nr ? $nf/$nr : 0)")
		avelen=$(perl -e "print ($nln ? exp($ln/$nln) : 0)")

		endog=$(perl -e "print ($nf+$l_nf ? $mapped_raw/($nf+$l_nf) : 0)")
		clon=$(perl -e "print ($mapped_raw ? ($mapped_raw-$mapped_dd)/$mapped_raw : 0)")

		echo -e "$sample\t$lib\t$nr\t$nf\t$foldratio\t$avelen\t$ref\t$mapped_raw\t$mapped_dd\t$endog\t$clon\t*\t*\t*\t*\t*\t*\t*\t*\t*"

		# stats for the ``LO'' runs of libraries.
		if [ "$l_nr" == "0" ]; then
			continue
		fi

		foldratio=$(perl -e "print ($l_nr ? $l_nf/$l_nr : 0)")
		avelen=$(perl -e "print ($l_nln ? exp($l_ln/$l_nln) : 0)")

		echo -e "$sample\t${lib}-LO\t$l_nr\t$l_nf\t$foldratio\t$avelen\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*"
	done

	# recalculate for sample summary
	foldratio=$(perl -e "print ($s_nr ? $s_nf/$s_nr : 0)")
	avelen=$(perl -e "print ($s_nln ? exp($s_ln/$s_nln) : 0)")
	endog=$(perl -e "print ($s_nf+$s_l_nf ? $s_mapped_raw/($s_nf+$s_l_nf) : 0)")
	clon=$(perl -e "print ($s_mapped_raw ? ($s_mapped_raw-$s_mapped_dd)/$s_mapped_raw : 0)")

	# coverage and depth
	set -- $(cov_depth ${sample}.${ref}.bam.stats.MT.txt $mt_len)
	mt_cov=$1; mt_dp0=$2; mt_dp1=$3
	set -- $(cov_depth ${sample}.${ref}.bam.stats.Aut.txt $aut_len)
	aut_cov=$1; aut_dp0=$2; aut_dp1=$3
	set -- $(cov_depth ${sample}.${ref}.bam.stats.chrX.txt $chrX_len)
	chrX_cov=$1; chrX_dp0=$2; chrX_dp1=$3

	echo -e "$sample\t*\t$s_nr\t$s_nf\t$foldratio\t$avelen\t$ref\t$s_mapped_raw\t$s_mapped_dd\t$endog\t$clon\t$mt_cov\t$mt_dp0\t$mt_dp1\t$chrX_cov\t$chrX_dp0\t$chrX_dp1\t$aut_cov\t$aut_dp0\t$aut_dp1"
done
done

rm $jf

