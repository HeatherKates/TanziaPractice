#!/bin/bash
#SBATCH --job-name=trimaglore.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb
#SBATCH --qos=zhangw
#SBATCH --account=zhangw
#SBATCH -t 24:00:00
#SBATCH --output=logs/trimaglore.out
#SBATCH --error=logs/trimgalore.err

module load trim_galore
#move to the parent directory so relative filepaths are accurate
cd ../
#trim_galore --phred33 --paired -fastqc --cores 2 --output_dir . 231-MSH2R4-1_S10_L001_R1_001.fastq.gz 231-MSH2R4-1_S10_L001_R2_001.fastq.gz 

#Run trimgalore 
trim_galore --phred33 --paired --fastqc --cores 2 --output_dir  trimgalore/ --basename MSH2R4-1 raw_reads/231-MSH2R4-1_S10_L001_R1_001.fastq.gz raw_reads/231-MSH2R4-1_S10_L001_R2_001.fastq.gz

#move the fastqc output to a different directory
mv *qc* trimgalore_fastqc

#run multiqc on fastqc files
module load multiqc
multiqc trimgalore_fastqc/
