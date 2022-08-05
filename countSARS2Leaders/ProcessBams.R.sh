#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -o sgeOut                        #-- output directory (fill in)
#$ -e sgeOut                        #-- error directory (fill in)
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l scratch=50G                   #-- SGE resources (home and scratch disks)
#$ -l h_rt=02:00:00                #-- runtime limit (see above; this requests 2 hours. With many genomes and/or bam files, this may not be enough)

#$ -l mem_free=16G                  #-- submits on nodes with enough free memory (required)  (will multiply by cores)




module load CBI
module load r

Rscript ProcessBams.R






## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"




