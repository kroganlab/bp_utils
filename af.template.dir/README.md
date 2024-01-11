# AlphaFold Jobs Template Directory

## Simple Usage
1. Make a copy of this directory
2. Make a new `masterFasta.fasta` file (match format to `masterFasta.fasta.example`)
3. Make a new `AlphaFoldJobList.csv` file (match format to `AlphaFoldJobList.example.csv`)
4. Edit submission script `af.jobs.sh`, most importantly the number of tasks, but also new file names or job names if desired
5. Submit job with `qsub af.jobs.sh`
6. View queue and job status with `qstat`

## Advanced Usage (to optimize usage of cluster)

After setting up the directory (through step 3 above), you will submit several jobs which make better use of non-GPU and shorter-time queues to optimize run time (at the expense of your time).

### 1) MSA management

MSA creation is included in the main alphafold python script, but it does not actually benefit from access to a GPU. MSA generation can take hours, so we want the MSA step to run on the more plentiful CPU nodes. To limit the jobs to just MSA generation, the script creates a special locked file that will cause alphafold to crash when it tries to write to it upon MSA completion.  This option is activated in job script `af.msa.sh`.  As above, edit the task numbers in the .sh, then submit with

```
qsub af.msa.sh
```

  
### 2) Jobs that can finish in 2 hours

Wynton GPU jobs that take at most 2 hours (set in the shell script) are allowed to run on more GPU nodes than longer jobs.  This is sufficient time for many AF jobs, especially with the MSA step already completed.  Use `af.smallJobs.sh`, edit the task array field and:

```
qsub af.smallJobs.sh
```

### 3) Submit jobs with (mostly) no time limit
Now we can submit `af.jobs.sh` as in the simple usage above.  Remember to edit the task array number.
```
qsub af.jobs.sh
```

### 4) For any incomplete jobs (any time after 2 has started), get scores and save MSAs
If all goes well, this is not necessary, but sometimes it is good to run separately to extract scores and save MSAs if things go wrong, or to look at already-completed scores before tasks finish for all 5 models.  Edit task array number, then

```
qsub af.postProcess.sh
```

## AlphaFold Output

Handy command to get all the values from the output/*/scores.csv files into a single file:

```
grep r output/*/scores.csv > allScores.csv
```

Then load into R with

```
library (data.table)
scores <- fread ("allScores.csv")
setnames(scores, c("path", "ptm", "iptm"))
scores[, pair := tstrsplit(path, "/")[[2]]]
scores[, meanIPTM := mean(iptm), by = pair]
scores[, pair := factor(pair, levels = unique(pair[order(-meanIPTM)]))]


```

