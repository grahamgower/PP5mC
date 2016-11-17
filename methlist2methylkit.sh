#!/bin/sh

# input:
# chrom	pos-0	pos-1	context	+C	+mC	-C	-mC

# output:
# https://github.com/al2na/methylKit

awk '
BEGIN {
  pfx = ARGV[2]
  if (pfx == "") {
    pfx = "methykit"
  }
  ARGV[2] = ""
  print "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT" > pfx".CpG.txt"
  print "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT" > pfx".CHG.txt"
  print "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT" > pfx".CHH.txt"
}

$4=="CpG" {
  dp = $5+$6+$7+$8;
  printf "%s.%s\t%s\t%d\t+\t%d\t%.2f\t%.2f\n",
  	 $1, $2, $1, $2+1, dp, 100*($6+$8)/dp, 100*($5+$7)/dp >> pfx".CpG.txt"
}
$4=="CHG" {
  dp = $5+$6+$7+$8;
  printf "%s.%s\t%s\t%d\t+\t%d\t%.2f\t%.2f\n",
  	 $1, $2, $1, $2+1, dp, 100*($6+$8)/dp, 100*($5+$7)/dp >> pfx".CHG.txt"
}
$4=="CHH" {
  f_dp = $5+$6;
  r_dp = $7+$8
  if (f_dp) {
    printf "%s.%s\t%s\t%d\t+\t%d\t%.2f\t%.2f\n",
  	 $1, $2, $1, $2+1, f_dp, 100*$6/f_dp, 100*$5/f_dp >> pfx".CHH.txt"
  } else if (r_dp) {
    printf "%s.%s\t%s\t%d\t-\t%d\t%.2f\t%.2f\n",
  	 $1, $2, $1, $2+1, r_dp, 100*$8/r_dp, 100*$7/r_dp >> pfx".CHH.txt"
  }
}
' - $@
