#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -o sgeOut                        #-- output directory (fill in)
#$ -e sgeOut                        #-- error directory (fill in)
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l scratch=50G                   #-- SGE resources (home and scratch disks)
#$ -l h_rt=00:29:00                #-- runtime limit (see above; this requests 29 minutes)

#$ -l mem_free=2G                  #-- submits on nodes with enough free memory (required)  (will multiply by cores)
#$ -pe smp 4         #--four slots(cores) on a single machine


#$ -t 1-4    #-- array job.  This will create 4 jobs, each with a different $SGE_TASK_ID. Edit this with the number of files in fastqList.txt


# update this to match mem_free and num slots above
memoryForBB=$(expr 2 \* ${NSLOTS:-1})

fastqFullPath=`awk -v taskID=${SGE_TASK_ID:-13} 'NR==taskID {print $1}' fastqList.txt`
echo $fastqFullPath

genome=`awk -v taskID=${SGE_TASK_ID:-13} 'NR==taskID {print $2}' fastqList.txt`
# echo
# echo `echo $genome | wc`
# if edited on a mac, the file may have an \r end of line character
genome=`echo $genome | tr -d '\r'`
genome=`echo $genome | tr -d '\n'`
# echo `echo $genome | wc`
echo $genome
echo


## 0. In case TMPDIR is not set, e.g. on development nodes, set
##    it to local /scratch, if it exists, otherwise to /tmp
if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

echo using temp directory at $TMPDIR

## 1. Use a temporary working directory
startDir=`pwd`
cd "$TMPDIR"

fastq=`basename $fastqFullPath`


prefix=`basename $fastqFullPath | cut -f1 -d"."`
date
echo $fastqFullPath
echo Finding leaders...
wLeader=`echo $fastq | sed s/.fastq.gz$/.wLeader.fastq.gz/`

echo cp $fastqFullPath $fastq
cp $fastqFullPath $fastq
/wynton/home/krogan/bpolacco/tools/bbmap/bbduk.sh -Xmx"$memoryForBB"g \
          threads="${NSLOTS:-1}" \
          literal=CTTTCGATCTCTTGTAGATCTGTTCTC \
          k=27 \
          maskmiddle=f \
          in=`basename $fastqFullPath` \
          outm=$wLeader

echo Trimming leaders...
trimmed=`echo $fastq | sed s/.fastq.gz$/.wLeader.trimmed.fastq.gz/`

/wynton/home/krogan/bpolacco/tools/bbmap/bbduk.sh -Xmx"$memoryForBB"g \
          threads=${NSLOTS:-1} \
          literal=CTTTCGATCTCTTGTAGATCTGTTCTC \
          k=27 \
          in=$prefix.wLeader.fastq.gz \
          ktrim=l \
          out=$trimmed

echo Matching leader-trimmed reads...
matched=`echo $fastq | sed s/.fastq.gz$/.wLeader.trimmed.mapped.sam/`
echo output: $matched

/wynton/home/krogan/bpolacco/tools/bbmap/bbmap.sh -Xmx"$memoryForBB"g \
    threads=${NSLOTS:-1} \
    32bit=t \
    ref=$genome \
    rebuild=t \
    in=$trimmed \
    out=$matched \
    mappedonly=t \
    strandedcov=t \
    maxindel=100 \
    strictmaxindel=t \
    local=t


echo Converting sam to bam...
bam=`echo $fastq | sed s/.fastq.gz$/.wLeader.trimmed.mapped.bam/`

/wynton/home/krogan/bpolacco/tools/samtools-1.15.1/samtools view -S -b $matched > $bam


echo copying results...
echo cp $wLeader "$startDir"/withLeaders/$wLeader
echo cp $trimmed "$startDir"/withLeaders/"$trimmed"
echo cp $matched "$startDir"/withLeaders/"$matched"
echo cp $bam "$startDir"/withLeaders/"$bam"

cp $wLeader "$startDir"/withLeaders/$wLeader
cp $trimmed "$startDir"/withLeaders/"$trimmed"
cp $matched "$startDir"/withLeaders/"$matched"
cp $bam "$startDir"/withLeaders/"$bam"



date
hostname


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

