# source this file, ie.
#   $ source phoenix_env.sh

module purge

module load BWA/0.7.15-foss-2016b
module load SAMtools/1.3.1-foss-2016b

module load GATK/3.7-Java-1.8.0_121
# NOTE: -Xmx parameter is per thread
java="java -Xmx2g"
gdir=$EBROOTGATK
export gatk="$java -jar $gdir/GenomeAnalysisTK.jar"

# pre-process 5mC pipeline
module load Python/2.7.13-foss-2016b
source /data/acad/programs/pp5mc-virtualenv/bin/activate
export PATH=/data/acad/programs/pp5mc:$PATH
export PATH=/data/acad/programs/pp5mc/slurm_pipeline:$PATH

# ensure files are group owned
chgrp acad_users .
chmod g+s .
umask 002
