#!/bin/sh

set -eu -o pipefail -o verbose

module load R/3.3.1-foss-2016b

for s in *.binStats.gz; do
	echo $s
	unpigz $s
	Rscript GlobalStatsPlotter.r ${s/\.gz/}
	pigz ${s/\.gz/}
done

