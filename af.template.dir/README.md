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

MSA creation is included in the main alphafold python script, but it does benefit from access to a GPU. Two approaches are possible here.  Submit regular jobs to non-GPU nodes and kill the jobs manually or let them die when time limit is reached (edit and submit job script `af.noGPU.sh`).  Or **recommended**, submit jobs with a special locked file that will cause alphafold to crash when it tries to write to it.  This is activated in job script `af.msa.sh`.  As above, edit the task numbers, then submit with

```
qsub af.msa.sh
```

If instead you opt for the `af.noGPU.sh` script, this command is handy to kill jobs that have already passed the MSA step based on the log files. Note that `af.noaf` should be replaced if you changed the job name in the submission script:

```
qdel  `grep  "Running model" af.noaf.* | cut -f1 -d":" | sed 's/af.noaf.o//'`
```

Jobs killed manually will not get a chance to copy their MSA to the MSA repository, so you may want to run `af.postProcess.sh` after "af.noaf.sh" if you don't use `af.msa.sh`

  
### 2) Jobs that can finish in 2 hours

Wynton GPU jobs that take at most 2 hours (set in the shell script) are allowed to run on more GPU nodes than longer jobs.  This is sufficient time for many AF jobs.  Use `af.smallJobs.sh`, edit the task array field and:

```
qsub af.smallJobs.sh
```

### 3) Submit jobs with (mostly) no time limit
Now we can submit `af.jobs.sh` as in the simple usage above.  Remember to edit the task array number.
```
qsub af.jobs.sh
```

### 4) For any incomplete jobs (any time after 2 has started), get scores and save MSAs
If all goes well, this is not necessary, but some times it is good to run separately to extract scores and save MSAs if things go wrong, or to look at already-completed scores before tasks finish for all 5 models.  Edit task array number, then

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

