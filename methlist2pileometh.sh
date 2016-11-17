#!/bin/sh

# input:
# chrom	pos-0	pos-1	context	+C	+mC	-C	-mC

# output:
# https://github.com/dpryan79/PileOMeth
# 1) The chromosome/contig/scaffold name
# 2) The start coordinate
# 3) The end coordinate
# 4) The methylation percentage rounded to an integer
# 5) The number of alignments/pairs reporting methylated bases
# 6) The number of alignments/pairs reporting unmethylated bases

awk '
BEGIN {
  pfx = ARGV[2]
  if (pfx == "") {
    pfx = "pileOmeth"
  }
  ARGV[2] = ""
  print "track type=\"bedGraph\" description=\"CpG methylation levels\"" > pfx".CpG.txt"
  print "track type=\"bedGraph\" description=\"CHG methylation levels\"" > pfx".CHG.txt"
  print "track type=\"bedGraph\" description=\"CHH methylation levels\"" > pfx".CHH.txt"
}

$4=="CpG" {
  dp = $5+$6+$7+$8;
  printf "%s\t%d\t%d\t%d\t%d\t%d\n",
  	 $1, $2, $3, 100*($6+$8)/dp, $6+$8, $5+$7 >> pfx".CpG.txt"
}
$4=="CHG" {
  dp = $5+$6+$7+$8;
  printf "%s\t%d\t%d\t%d\t%d\t%d\n",
  	 $1, $2, $3, 100*($6+$8)/dp, $6+$8, $5+$7 >> pfx".CHG.txt"
}
$4=="CHH" {
  f_dp = $5+$6;
  r_dp = $7+$8
  if (f_dp) {
    printf "%s\t%d\t%d\t%d\t%d\t%d\n",
  	 $1, $2, $3, 100*$6/f_dp, $6, $5 >> pfx".CHH.txt"
  } else if (r_dp) {
    printf "%s\t%d\t%d\t%d\t%d\t%d\n",
  	 $1, $2, $3, 100*$8/r_dp, $8, $7 >> pfx".CHH.txt"
  }
}
' - $@

