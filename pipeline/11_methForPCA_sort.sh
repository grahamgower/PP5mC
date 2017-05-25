#!/bin/sh

#Output: autosomes only, sorted with sort -k, input already 0 base start

set -eu -o pipefail -o verbose

methCov1F() {
        f=$1
	unpigz -c $f \
        | awk -v OFS='\t' 'NR>1 && $2~/^chr[0-9]/ && $4~/F/ {print $1,$6}' \
        | sort -k1,1 \
        | pigz \
        > ${f/\.txt\.gz/.txt.1covF.gz}
}

methCov3F() {
        f=$1
        unpigz -c $f \
        | awk -v OFS='\t' 'NR>1 && $2~/^chr[0-9]/ && $4~/F/ && $5>=3 {print $1,$6}' \
        | sort -k1,1 \
        | pigz \
        > ${f/\.txt\.gz/.txt.3covF.gz}
}

methCov5F() {
        f=$1
	unpigz -c $f \
        | awk -v OFS='\t' 'NR>1 && $2~/^chr[0-9]/ && $4~/F/ && $5>=5 {print $1,$6}' \
        | sort -k1,1 \
        | pigz \
        > ${f/\.txt\.gz/.txt.5covF.gz}
}

methCov1R() {
	f=$1
	unpigz -c $f \
	| awk -v OFS='\t' 'NR>1 && $2~/^chr[0-9]/ && $4~/R/ {print $1,$6}' \
	| sort -k1,1 \
	| pigz \
	> ${f/\.txt\.gz/.txt.1covR.gz}
}

methCov3R() {
	f=$1
	unpigz -c $f \
	| awk -v OFS='\t' 'NR>1 && $2~/^chr[0-9]/ && $4~/R/ && $5>=3 {print $1,$6}' \
	| sort -k1,1 \
	| pigz \
	> ${f/\.txt\.gz/.txt.3covR.gz}
}

methCov5R() {
	f=$1
	unpigz -c $f \
	| awk -v OFS='\t' 'NR>1 && $2~/^chr[0-9]/ && $4~/R/ && $5>=5 {print $1,$6}' \
	| sort -k1,1 \
	| pigz \
	> ${f/\.txt\.gz/.txt.5covR.gz}
}


for s in *.Cow_UMD3_1.methylkit.CpG.txt.gz; do
	echo $s
	methCov1F $s &
	methCov3F $s &
	methCov5F $s &
	methCov1R $s &
	methCov3R $s &
	methCov5R $s
	wait
done


