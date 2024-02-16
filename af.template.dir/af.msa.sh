#!/bin/env bash
#
##$ -q gpu.q
#$ -N af.msa 
#$ -cwd
#$ -l h_rt=24:00:00
##$ -l h_rt=48:00:00
#$ -l mem_free=60G
#$ -l scratch=50G
##$ -l compute_cap=80,gpu_mem=40G

#$ -t 1-100             ## job array with xx tasks

#$ -j y
#$ -o jobLogs/$JOB_NAME-$JOB_ID-$TASK_ID.log


# if not running with sge task array, set to 5
taskID=${SGE_TASK_ID:-5}



#
# Compute cap for A100 GPU is 8.0 (40 or 80 GB), for A40 GPU is 8.6 (48 GB).
#

t0=$(date --rfc-3339=seconds)

echo "QUEUE: $QUEUE"
echo "SGE_GPU: $SGE_GPU"
export CUDA_VISIBLE_DEVICES=$SGE_GPU



./AF_saveMSAS.231.py --model_preset=multimer --job_id=$taskID \
        --master_fasta=masterFasta.fasta \
        --jobTable=AlphaFoldJobList.csv \
        --prevent_alphafold_output=True

t1=$(date --rfc-3339=seconds)
echo "Duration: $t0 -- $t1"
