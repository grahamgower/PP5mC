Pre-processing reads is a multi-stage endeavour that encompasses read merging,
mapping, and calling methylation state.  This may be computationally
intensive, and thus it is advantageous to do the processing on a cluster.
This folder contains scripts implementing a flexible pipeline for doing such
processing using the [Slurm](https://slurm.schedmd.com/) workload manager.

Pipeline stages:

1.  Recovery of original sequence from Hairpin-ligated Bisulfite-treated reads
    (`PP5mC/foldreads`).
2.  Map reads (`BWA-mem`).
3.  Merge sequencing runs for the same library (`samtools merge`).
4.  Duplicate removal (`PP5mC/rmdup_collapsed.py` imported from
    [Paleomix](https://github.com/MikkelSchubert/paleomix)).
5.  Merge libraries for the same sample (`samtools merge`).
6.  Realign around indels (1. `$gatk -T RealignerTargetCreator`,
    2. `$gatk -T IndelRealigner`, 3. `samtools calmd`).
7.  Mark C/5mC states for each locus, from the pileup (`PP5mC/mark_5mC`).
8.  Summary statistics and plots.

The pipeline does not do variant calling.

Installation
============

1.  Install the dependencies:

	* [Samtools](http://samtools.github.io/)
	* [BWA](https://github.com/lh3/bwa)
	* [GATK](https://software.broadinstitute.org/gatk/)

2.  Clone the PP5mC git repository, and build the C source files (type `make`).

3.  Export the `gatk` environment variable, specifying how to run GATK (see note below).

Configuring the pipeline
========================
Configuration is separated into two logical components, data to be processed,
and cluster resource usage for each pipeline stage.  Both components are
specified in a [JSON](https://tools.ietf.org/rfc/rfc7159.txt) format
configuration file.

Data specification
------------------
A "samples" object contains stanzas for each sample to be processed.
Each sample stanza may contain one or more library, and each library contains
one or more sequencing run.  Each sequencing run specifies two input files,
the first must be the read 1 fastq file (R1) and the second must be the
read 2 fastq file (R2).  It is assumed that these files contain the strings
'R1' and 'R2'.

If all fastq files are contained under a common subdirectory, this may be
specified with the "prefix" object.  Sample ids, library ids, and run ids are
all free form fields wherin punctuation is allowed.  However, these ids will
be used in filenames, so there may be limitations imposed by the filesystem.

```
    "prefix": "/path/to/fastq/files",
    "samples": {
        "Sample1": {
            "Lib1": {
                "MiSeq-20170801" : [
                     "miseq-20170801/sample1-lib1-R1.fastq.gz",
                     "miseq-20170801/sample1-lib1-R2.fastq.gz"
                ],
                "HiSeq-20171004" : [
                     "hiseq/sample1-lib1-20171004-R1.fastq.gz",
                     "hiseq/sample1-lib1-20171004-R2.fastq.gz"
                ]
            },
            "Lib2": {
                "MiSeq-20170801" : [
                     "miseq-20170801/sample1-lib2-R1.fastq.gz",
                     "miseq-20170801/sample1-lib2-R2.fastq.gz"
                ]
            }
        },
        "Sample2": {
            "Lib1": {
                "MiSeq-20170801" : [
                     "miseq-20170801/sample2-lib1-R1.fastq.gz",
                     "miseq-20170801/sample2-lib1-R2.fastq.gz"
                ]
            }
        }
    },
```

A "refs" object specifies the set of genome references to which the data should
be aligned.  Each reference is given a "key" which will be used in filenames
to indicate the data corresponds to a particular reference.

```
    "refs": {
        "UMD311": "/path/to/refs/Bos_taurus/UMD_3.1.1/UMD311.fasta"
        "Lambda": "/path/to/refs/Lambda/Lambda.fasta"
    },
```


Cluster resources
-----------------

A basic familiarity with [Slurm](https://slurm.schedmd.com/) is assumed.

Slurm provides a convenient framework for cluster administrators to limit
resource utilisation such that each user is given fair access to resources.
Unfortunately, this also makes the job of the user more complex.  In
particular, it is not possible *a priori* to make accurate guesses about the
resources that are required for a given job.  The best one can do is run a
small job on a subset of one's data, requesting more resources than are
thought to be required, then look at the resources that were actually used and
scale that up for a job with the full amount of data.  We first describe how
to specify the resource utilisation in the PP5mC pipeline, then discuss how
a test job can be created for fine tuning resource requests on your cluster.

### Specifying resources ###
In the config file, a "resources" object contains a stanza for each stage of
the pipeline.  Within each stanza, resources may be specified as name-value
pairs.  The resources that may be requested are memory usage ("mem"), cpu
count ("cpus"), and run time ("secs", "mins", or "hours").  Each resource may
be specified using a fixed value, or as a value that will be mulitplied by
the size of the input data (the sum of all relevant fastq.gz input files).
Resource names starting with an underscore are fixed, while those without an
underscore correspond to the resource units per Mb of input data.

The "queue" object specifies a Slurm queue to which the jobs will be submitted.
If none is specified, the default queue will be used.

```
    "queue": "batch",
    "resources": {
	    "01:fold": {"_mem":"16M", "_cpus":1, "_mins":10},
	    "02:map": {"_mem":"6G", "_cpus":4, "_mins":10},
	    "03:mergeruns": {"_mem":"128M", "_cpus":1, "_mins":10},
	    "04:dedup": {"_mem":"128M", "_cpus":1, "_mins":10},
	    "05:mergelibs": {"_mem":"128M", "_cpus":1, "_mins":10},
	    "06:realign1": {"_mem":"3G", "_cpus":2, "_mins":10},
	    "06:realign2": {"_mem":"3G", "_cpus":1, "_mins":10},
	    "06:realign3": {"_mem":"128M", "_cpus":1, "_mins":10},
	    "07:mark5mC": {"_mem":"128M", "_cpus":1, "_mins":10}
    },
```

Actually, the PP5mC pipeline does not require the full generality provided
by this framework.  Typically, the "_mem" object should be
used, as all pipeline stages have a fixed memory overhead, or overhead that
does not scale linearly with the job size.  The "_cpus" object should be used,
to specify a single cpu for all stages, except "03:map" and "06:realign1",
where "_cpus" should indicate the desired number of threads for `BWA` and
`$gatk -T RealignerTargetCreator` respectively.

We suggest that resource scaling be used only for the time objects ("secs" or
perhaps "mins").  See section below on running a test pipeline, in order to
obtain reasonable values, including a scaled value for "secs".

### A note regarding java and GATK ###

The PP5mC pipeline uses subcommands provided by the
[Genome Analysis ToolKit (GATK)](https://software.broadinstitute.org/gatk/),
which is implemented in java.  Hence, you will need java and GATK installed
to use the pipeline.  The pipeline uses the `$gatk` environment variable to
call GATK, which you should export into the environment from which you're
running PP5mC. I.e. before using the pipeline, run something like:
```
export gatk='java -Xmx2g -XX:ParallelGCThreads=2 -jar /path/to/GenomeAnalysisTK.jar'
```

The maximum heap size for a java application can be set using `-Xmx<M>`, where
`<M>` indicates the amount of memory to reserve for the process.  The command
above reserves 2 Gb of heap space.  The *heap* is just a name for the space
from which memory is drawn, as the application requests it.  If the java
application requests more memory than this, those requests will fail and the
application will likely terminate with an *out of memory* error.  In this
case, the maximum heap space must be increased.  The *stack* is another area
of memory, distinct from the heap, that is typically used for holding an
application's local variables.  The stack need not be as large as the heap (in
general), but each application thread is given its own stack space, whereas
application threads share the heap.  The maximum stack size can be altered,
but this is not usually necessary.  The default stack size is on the order of
1 Mb, and thus contributes little to overall memory usage unless the
application uses a very large number of threads.

It is important to be aware that the [java virtual machine consumes
resources](http://www.oracle.com/technetwork/java/javase/gc-tuning-6-140523.html)
separately to those used by the java application.
Perhaps unexpectedly, specifying the maximum heap space does not provide an
upper limit on the amount of memory that will be used when running a java
application.  The default garbage collector implemented in java
(`ParallelGCThreads`) will create 1 thread per CPU up to 8 CPUs + 5/8 per CPU
afterwards.  This means that on a 32 core cluster node, java will spawn 23
threads behind your back.  What's more, each garbage collection thread has a
non-negligible memory overhead---distinct from the heap space.
We strongly recommend limiting the maximum heap size (`-Xmx`) and the number
of garbage collector threads (`-XX:ParallelGCThreads`).  The "_mem" resource
object for the "06:realign1" and "06:realign2" pipeline stages should be set
*higher* than the value specified with `-Xmx` to account the java virtual
machine overhead.  An additional 1-2 Gb ought to be sufficient, depending on
the number of garbage collector threads.


Running the pipeline
--------------------

After writing a configuration file (e.g. `config.json`), first do a "dry run"
of the pipeline.  This will not submit any jobs to Slurm, but will perform
basic checks on the configuration, create a job dependency graph, and write job
scripts to the `acct/` folder.

```
pp5mC.py --dryrun config.json
```

We recommend looking at the `dependency_graph.pdf` output file to confirm the
configuration conforms to expectation, prior to launching the pipeline with:

```
pp5mC.py config.json
```

### Files/folders and hierarchy ###

The files created by the pipeline are:

dependency_graph.pdf
  ~ A graph showing which pipeline tasks depend on which
    other tasks, and the resources requested for each.

acct/
  ~ Folder containing job scripts and Slurm output.

Sample1/
  ~ Folder containing intermediate files pertaining to sample named "Sample1"

Sample1/Lib1
  ~ Folder containing intermediate files pertaining to library named "Lib1"
    belonging to sample "Sample1"

Sample1.UMD11.bam
  ~ Final alignment file for sample "Sample1", aligned to reference "UMD11".

Sample1.UMD11.methlist.txt.gz
  ~ Tab separated file containing counts of C/5mC calls for each cytosine
    in the Sample1.UMD11.bam alignment.

Sample1.UMD11.methylkit.{CpG,CHG,CHH}.txt.gz
  ~ [Methylkit](https://github.com/al2na/methylKit) format output files.

Sample1.UMD11.pileOmeth.{CpG,CHG,CHH}.txt.gz
  ~ [PileOmeth](https://github.com/dpryan79/MethylDackel) (now MethylDackel)
    format output files.

### Running a test pipeline ###

In order to best utilise cluster resources, run the pipeline on a subset
of your data first, using a relaxed set of resource constraints.
We provide a script (`res.py`) to be run on the output of 

```
```


Troubleshooting
---------------

This software has no bugs. :P
