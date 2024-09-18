# Reproducible workflow to run ATACseq pipeline on hipergator

This workflow is to practice running an ATACseq workflow (sequence QC, trimming, alignment, peak calling) on a single pair of read files. 

Anywhere that you see text inside ```<>```, it must be replaced by a value specific to you (e.g. your username, filepath) **without** the ```<>```

### The ATACseq workflow can be replicated by running the scripts in scripts/ sequentially

## Step 1: Log onto hipergator, go to your main /blue directory, and copy this github repository

```ssh <username>@hpg2.rc.ufl.edu```

```cd </path/to/your/main/dir>``` # For you this is ```/blue/zhangw/<username>```

```git clone https://github.com/HeatherKates/TanziaPractice.git```

## Step 2: Change to the directory where the github repository is copied and see what's there

```cd TanziaPractice```

```tree```

## Step 3: Change to the directory where the scripts are so you can run the  scripts from there

```cd scripts/```

## Step 4: Run the first script to copy a single mate pair from /orange

Before you run it, look at the script to see what it will do

```more 1_copy_raw_reads.sbatch```

To submit the job and copy the reads, run:

 ```sbatch 1_copy_raw_reads.sbatch```

After you've submtited the job, check the status at any point:

```squeue -u <username>```

When your job is complete, check that the reads were copied successfully:

 ```ls ../raw_reads```

## Step 5: (ATACseq step 1) Run trimgalore on the raw reads to detect and trim contaminating adapters from the reads

Look at	the script before running to see what it is doing:

```more 2_trim_galore.sbatch```

Note the ```fastqc``` option means that fastqc will be run before and after trimming to generate the fastQC reports, so we do not need to do that step separately. We also included `multiqc` in the ```*sbatch``` script to summarize the results of any fastqc reports.

Submit the job to run trimgalore on your reads:

```sbatch 2_trim_galore.sbatch```

Check the status of your job at any point

```squeue -u <username>```

When your job is complete, check that the expected output was generated:

```ls ../trimgalore```

To view the fastqc output from trimgalore, use ```scp``` to transfer the file to your computer:

```scp <username>@hpg2.rc.ufl.edu:</path/to/multiqc_report.html> </path/to/destination>```

## Step 6: For ATACseq steps 2-4, follow the Step 5 steps for the remaining scripts in scripts/ one at a time 
### (do not start the next job before the previous is completed successfully)

To find where each step prints output, look in the ```*sbatch``` script as described above (hint, using ```more```)

When all steps are complete, check these output files to assess the success of alignment:

```more ../flagstat/aligned.bam.flagstat.log```

```more ../flagstat/aligned_filtered_sorted_duprmv.bam.flagstat.log```

Use ```more 4_samtools.sbatch``` to see how these files were created and what information they contain.

* ```flagstats/aligned.bam.flagstat.log``` is alignment stats of the inital alignment output by bowtie2.

* The stats of the final sam file output by ```4_samtools.sbatch``` are recorded in ```flagstats/aligned_filtered_sorted_duprmv.bam.flagstat.log```

Compare the results in these two files to assess the overall initial mapping of trimmed reads to the reference as well as how subsequent filtering affects reads present in the final alignment to be used for peak calling.
