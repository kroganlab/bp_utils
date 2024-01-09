# AlphaFold Jobs Template Directory

## Simple Usage
1. Make a copy of this directory
2. Make a new `masterFasta.fasta` file (match format to `masterFasta.fasta.example`)
3. Make a new `AlphaFoldJobList.csv` file (match format to `AlphaFoldJobList.example.csv`)
4. Edit submission script `af.jobs.sh`, most importantly the number of tasks, but also new file names or job names if desired
5. Submit job with `qsub af.jobs.sh`
6. View queue and job status with `qstat`

## Advanced Usage (to optimize usage of cluster)
****

  
