## Reproducible workflow to run ATACseq pipeline on a single mate pair of read files

# The ATACseq workflow can be replicated by running the scripts in scripts/ sequentially

# Chage to the directory where the scripts are and execute scripts from there

```cd scripts/```

# Copy a single mate pair from /orange

Look at the script before running to see what it is doing using ```more 1_copy_raw_reads.sbatch```

To submit the job and copy the reads, run ```sbatch 1_copy_raw_reads.sbatch```

Check the status of your job using ```squeue -u <username>```

When your job is complete, check that the reads were copied using ```ls ../raw_reads```

# Run trimgalore on the raw reads to detect and trim contaminating adapters from the reads

Look at	the script before running to see what it is doing using	```more 2_trim_galore.sbatch```

To submit the job and copy the reads, run ```sbatch 2_trim_galore.sbatch```

Check the status of your job using ```squeue -u	<username>```

When your job is complete, check that the reads	were copied using ```ls ../trimgalore```

View the fastqc output from trimgalore that will provide information about QC metrics before and after trimming:
* use ```scp``` to transfer the multiqc_report.html file to your computer: 
* ```scp <username>@hpg2.rc.ufl.edu:/path/to/multiqc_report.html /path/to/destination```

# Follow the four steps above for the remaining scripts in scripts/ one at a time (do not start the next job before the previous is completed successfully)

To find where each step prints output, look in the *sbatch script as described above (hint, using ```more```)

When all steps are complete, check these files to assess the success of alignment:

```more ../flagstat/aligned.bam.flagstat.log```

```more ../flagstat/aligned_filtered_sorted_duprmv.bam.flagstat.log```

Use ```more 4_samtools.sbatch``` to see how these files were creaed and what information they contain.

aligned.bam.flagstat.log is alignment stats of the inital alignment output by bowtie2.

Many processing steps are executed by 4_samtools.sbatch, and the stats of the final sam file are recored in aligned_filtered_sorted_duprmv.bam.flagstat.lo

Compare the results in these two files to assess the overall initial mapping of trimmed reads to the reference as well as how subsequent filtering affects reads present in the final alignment to be used for peak calling.
