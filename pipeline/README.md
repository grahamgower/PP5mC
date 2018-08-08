This folder contains a Bourne shell script pipeline for preprocessing of
HBS-seq data.  It is not intended to be a 100% robust pipeline, but instead
provides a reasonable starting point for local customisations.  It has been
tested only on a CentOS Linux system, but should work with few if any changes
in other environments (e.g. MacOSX).

Start by editing `fold.settings`, and **PLEASE** read the scripts before
running them.  In particular, fastq naming conventions may not match your data,
and paths to programs such as GATK likely need to be changed.
