#!/bin/bash


if [ -z $1 ]; then
        echo 'Parameter 1 should be a directory containing *.fastq.gz files'
        exit 0
else
    fastqDir="$1" #/wynton/group/krogan/apelin/RNA-seq/For_others/Exp08_NY3/AdaptTrim
    echo "Will process files in $fastqDir"
fi



# copy the sge submission scripts, R script, reference fasta from source
# the sge sumbission script  will need some manual editing of the task ID line, but otherwise is ready
if [ -f "leaderReadsProcess.sh" ]; then
    echo "working with files in current directory"
else 
    echo "copying files from ../countSARS2Leaders/"
    cp -r ../countSARS2Leaders/* ./
fi

# the script depends upon a file which lists the fastq files, make that:
ls  $fastqDir/*fastq.gz > fastqList.txt
echo
echo __IMPORTANT__
echo fastqList.txt file created. You must edit it with a second column with genome paths before proceeding.

fastqCount=`grep -c . fastqList.txt`
echo 
echo __IMPORTANT__
echo $fastqCount fastq files found. Edit the file leaderReadsProcess.sh with  this as the upper bound of tasks
echo example: $'\n\t' '#$ -t 1-'$fastqCount 
echo 

# make a couple of directories that the sge script needs
mkdir -p sgeOut
mkdir -p withLeaders


echo 
echo after editing leaderReadsProcess.sh and fastqList.txt, process leaders with two sge jobs
echo first submit:
echo    qsub leaderReadsProcess.sh
echo
echo then when the above is complete, submit: 
echo    qsub ProcessBams.R.sh
echo

