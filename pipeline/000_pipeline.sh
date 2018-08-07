#!/bin/sh

trap 'kill 0' EXIT

./01_fold.sh
./02_map.sh
./03_merge.sh
./04_dedup.sh
./05_indel_realign.sh
./06_scanbp.sh
./07_mark5mC.sh
