#!/bin/bash
#SBATCH --job-name=trimaglore.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb
#SBATCH --qos=zhangw
#SBATCH --account=zhangw
#SBATCH -t 24:00:00
#SBATCH --output=../logs/trimaglore.out
#SBATCH --error=../logs/trimgalore.err

module load trim_galore
module load fastqc
#move to the parent directory so relative filepaths are accurate
cd ../

#Run trimgalore 
trim_galore --phred33 --paired --fastqc --cores 2 --output_dir  2a_trimgalore/ --basename MSH2KO-1 1_raw_reads/231-MSH2KO-1_S7_L001_R1_001.fastq.gz 1_raw_reads/231-MSH2KO-1_S7_L001_R2_001.fastq.gz

#Run fastqc on input files
fastqc -o 2b_trimgalore_fastqc/ 1_raw_reads/231-MSH2KO-1_S7_L001_R1_001.fastq.gz 1_raw_reads/231-MSH2KO-1_S7_L001_R2_001.fastq.gz

#move the fastqc output to a different directory
mv 2a_trimgalore/*qc* 2b_trimgalore_fastqc/

#run multiqc on fastqc files
module load multiqc
multiqc 2b_trimgalore_fastqc/ -o 2b_trimgalore_fastqc/
