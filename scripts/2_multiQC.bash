#!/bin/bash
module load multiqc
# Run multiqc to summarize the fastqc results
multiqc 2b_trimgalore_fastqc/ -o 2b_trimgalore_fastqc/
