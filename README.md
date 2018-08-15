# PP5mC

PP5mC is a pre-processing toolkit for reconstructing target molecules
and calling 5-methylcytosine (5mC) state, for hairpin-ligated
bisulfite-treated sequencing (HBS-seq) experiments.

# Prerequisites
The core parts of PP5mC are written in C, and thus requires a C compiler
(tested with **GCC** v5.5.0).  The Slurm pipeline, and several helper scripts
were written in **Python** (tested with **v2.7.15**).  While PP5mC was
developed and tested only on Linux, its components are expected to work
on \*BSD, MacOSX, etc. without major modifications (feedback welcome!).
In addition, the following libraries are required:

* [zlib](https://zlib.net/), used by `foldreads` for decompressing gzipped
  fastq files;
* [htslib](https://github.com/samtools/htslib), used by `scanbp` and `mark5mC`
  to parse bam files; and
* [matplotlib](https://matplotlib.org/), used by plotting scripts.

In addition, the full processing pipelines require
[BWA](https://github.com/lh3/bwa),
[samtools](http://www.htslib.org/), and
[GATK](https://software.broadinstitute.org/gatk/).


# Installation
Clone the git repository, then build with `make`.  You may wish to edit
the `Makefile` to statically link with htslib, and/or to specify its path.

# Usage

The [shell script pipeline](pipeline/) for unix workstations implements a
minimal working pipeline that:
(1) reconstructs reads with `foldreads`;
(2) aligns reconstructed reads with BWA-mem;
(3) merges different libraries for the same sample;
(4) removes PCR duplicates;
(5) realigns reads around indels using `GATK`;
(6) records and plots nucleotide pairing frequencies with `scanbp`; and
(7) calls methylation status at all CpG, CHG, and CHH contexts with `mark5mC`.

The Slurm pipeline follows roughly the same stages, and is further documented
[here](slurm_pipeline/README.md).

Most of the tools included with PP5mC provide some usage information if
they are invoked without any parameters.


## Reconstructing original target molecules
The `foldreads` command attempts to "fold" input sequences back together
at the hairpin, thereby reconstructing the original sequences.  It takes two
HBS-seq FASTQ files as input, and also outputs the reconstructed sequences
in FASTQ format.
Multiple hairpin sequences may be specified (`-p`/`-P` parameters can be
specified multiple times), and `foldreads` searches for hairpins in the order
in which they are specified.

```
$ ./foldreads 
Error: must specify input files.
foldreads v12
usage: ./foldreads [...] -1 IN1.FQ -2 IN2.FQ
 -o OUT.FQ         Fastq output file [stdout]
 -m FILE           Metrics output file [stderr]
 -u PREFIX         Filename prefix for unfolded reads []
 -p SEQ            The hairpin SEQuence [ACGCCGGCGGCAAGTGAAGCCGCCGGCGT]
 -P SEQ            The hairpin SEQ (unmethylated and will be C->T converted)
 -1 IN1.FQ[.GZ]    R1 fastq input file
 -2 IN2.FQ[.GZ]    R2 fastq input file
 -T SEQ            Adapter SEQuence trailing R1 (p7) []
 -B SEQ            Adapter SEQuence trailing R2 (p5) []
```

`foldreads` also collects some metrics from the reads that it processes,
which are printed to ``stderr`` by default, or to a file if specified with
the `-m` parameter.  In particular, `foldreads` records the number of reads
processed (NR), the number of reads successfully reconstructed (NF), the
number in which the full hairpin sequence was observed (NH), the number
of occasions the hairpin was found in different positions in read 1 and
read 2 (ND), and the number of hairpins identified for each hairpin sequence
specified with `-p` or `-P`.

A read length histogram is also recorded, and is fitted to a truncated
lognormal distribution.  Observed fragment lengths are no longer than the
read length (and may be shorter if a hairpin is found, or if bases are trimmed
due to coincidental matches with a hairpin sequence), so the distribution
is truncated near the read length.  Fitting to a lognormal works well when
the fragment length distribution is not truncated by the read length too
heavily  (e.g. see ``test_fit.c``).

## Read alignment and optional SAM fields.
Reads are aligned using `bwa mem -C`.
The [BWA manual](http://bio-bwa.sourceforge.net/bwa.shtml) states the
following about using the `-C` flag:
```
-C 	Append FASTA/Q comment to SAM output. This option can be used to
	transfer read meta information (e.g. barcode) to the SAM output. Note
	that the FASTA/Q comment (the string after a space in the header line)
	must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead
	to incorrect SAM output. 
```

We note that `bwa aln` contains no equivalent flag, however we have submitted
a [pull request](https://github.com/lh3/bwa/pull/203) to add this feature.
The FASTQ file output by `foldreads` contains additional fields following the
@\<read name\>.  After being transferred into the SAM/BAM file by `bwa`, the
additional fields are used by `scanbp` and `mark5mC`.

```
r1:Z:<read 1 sequence>
q1:Z:<PHRED+33 quality scores for read 1>
r2:Z:<read 2 sequence>
q2:Z:<PHRED+33 quality scores for read 2>
hp:Z:<hairpin sequence specified with -p or -P that matched this read>
```

## M-bias and methylation calling
Because HBS-seq produces observations for both the top and bottom strands of
target molecules, we take advantage of the additional information when creating
M-bias plots.  These plots are used to look at the frequency of methylation
calls as a function of read position, and can indicate sequencing artefacts
such as ligation bias.  The `scanbp` command uses the read pileup to record
the frequency of different nucleotide pairs at each position within reads,
as well as positions upstream and downstream of where the read mapped.

```
$ ./scanbp file.bam > ntpairs.txt
$ ./plot_nt_pairing.py ntpairs.txt ntpairs.pdf
```

Regular BS-seq tools are not directly usable for methylation calling from
HBS-seq data, so we provide the `mark5mC` command to do this.  From a visual
inspection of the nucleotide pairing plots, one can identify how many bases
at each end of the reads where methylation calls may not be trustworthy, and
specify their exclusion with `-5` and `-3` parameters to `mark5mC`.

```
$ ./mark5mC -5 10 -3 10 file.bam ref.fasta | gzip -c - > methlist.txt.gz
```

The output file is a (bed compatible) tab separated file, indicating the
``depth``, and number of ``C``/``mC`` observations, and the
``strand``/``context`` for the observations.  The ``C`` and ``mC`` columns
may not sum to the ``depth``, due to heterozygosity or differences between
the reference and the sequenced individual.
```
$ zcat methlist.txt.gz | head
chrom   pos-0   pos-1   strand  depth   C       mC      context
chr1    3000241 3000242 +       1       1       0       CAT
chr1    3000244 3000245 +       1       1       0       CCA
chr1    3000245 3000246 +       1       1       0       CAG
chr1    3000247 3000248 -       1       1       0       CTG
chr1    3000248 3000249 -       1       1       0       CCT
chr1    3000253 3000254 +       2       2       0       CCT
chr1    3000254 3000255 +       2       2       0       CTG
chr1    3000256 3000257 -       2       2       0       CAG
chr1    3000257 3000258 -       2       2       0       CCA
```

This list of methylation calls may be converted to the
[pileOmeth](https://github.com/dpryan79/MethylDackel) format, or the
[methylkit](https://github.com/al2na/methylKit) format,
using the `frobmethlist.py` script provided.

## Simulating HBS-seq reads

### Quality score profiles

### ``simhbs``
