# AlphaFold Jobs Template Directory

## Simple Usage
1. Make a copy of this directory. This copy is your working directory
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

### NEW 2.a) 2hr incremental jobs start where 2 hour jobs left off (mostly)

With a modified singularity image, jobs can build starting from the first missing pdb. To do this use the setup\_singularity\_sandbox.sh file to create a new singularity image, then use af.incremental.sh.

```
cd af.template.dir
bash setup_singularity_sandbox.sh
mv ../alphafoldSandbox <someplacePermanent>
```


```
cd workDir
# edit af.incremental.sh with path to your alphafoldSandbox
qsub af.incremental.sh
```

repeat the above as needed, monitoring how much work is actually getting done. For those jobs where 2 hours is not enough to produce a model (and thus never advance), you will still have to run (3)

A handy way to submit jobs on wynton is to use job dependencies. In one session you can submit the msa job array and then 5 copies of the incremental af array, each can depend on completion of the corresponding task in the previous job.  Like so:

```
cd workDir
qsub af.msa.sh
# pay attention to the above job id assigned for each qsub
# and replace each 99999* below with the preceding job ID
qsub -hold_jid_ad 999991 af.incremental.sh  
qsub -hold_jid_ad 999992 af.incremental.sh  
qsub -hold_jid_ad 999993 af.incremental.sh  
qsub -hold_jid_ad 999994 af.incremental.sh  
qsub -hold_jid_ad 999995 af.incremental.sh  

```
`hold_jid_ad` says to hold individual tasks until completion of the corresponding task of the specified job_id.  If instead you use `hold_jid` your whole task array will wait for completion of all tasks of the specified job_id.

### 3) Submit jobs with (mostly) no time limit
Now we can submit `af.jobs.sh` as in the simple usage above.  Remember to edit the task array number.
```
qsub af.jobs.sh
```
Or better, copy and modify the `af.incremental.sh` file to run just those jobs that are incomplete. This command will find the directories in `output` that do not have the 5th relaxed pdb, and write a new `incompleteJobList.csv`:

```
find output  -mindepth 1 -maxdepth 1 -type d '!' -exec test -e "{}/relaxed_model_5_multimer_v3_pred_0.pdb" ';' -print | sed 's/__/\//' | awk -F"/"   'BEGIN{OFS = ","}{print NR,$2,$3;}' > incompleteJobList.csv
```
Then copy and update af.incremental.sh with the new number of tasks (the highest row number in `incompleteJobList.csv`), and make sure to point the af python command at the new `incompleteJobList.csv` file.

```
./AF_saveMSAS.231.py --model_preset=multimer --job_id=$taskID \
        --singularity_image=$singularitySandbox \
        --check_models_complete=True \
        --master_fasta=masterFasta.fasta \
        --jobTable=incompleteJobList.csv

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

Work in progress here, but this is handy for pulling out the useful information from the run logs:
```
grep -E "(Running model)|(seed)|(Total JAX model)|(Running monomer pipeline)|(MSA size)|(Total number of templates)"  jobLogs/*
```

