#!/bin/sh


trap 'kill 0' EXIT

./01_fold.sh
./02_map.sh
./03_merge.sh
./04_dedup.sh
./05_indel_realign.sh
./06_mark_5mC.sh
./07_scan_pairs.sh
./08_conversionStats.sh
./09_mappingStats.sh
./10_binStats.sh
./11_methForPCA.sh
./12_GlobalStatsPlotter.sh
