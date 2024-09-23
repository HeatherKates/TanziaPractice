#!/bin/bash
#SBATCH --job-name=bowtie2.%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --time=24:00:00
#SBATCH --output=../logs/bowtie2_%A_%a.out
#SBATCH --error=../logs/bowtie2_%A_%a.err
#SBATCH --array=0-11  # Based on the sample list size

# Load Bowtie2
module load bowtie2/2.4.5

# Define the base path for trimmed data
base_trim_path="../2a_trimgalore"

# Read the sample name from the sample list
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" MSH2_sample_list.txt)

# Define the paths for trimmed files (L001 and L002)
R1_trim_L001="${base_trim_path}/${sample}_L001_val_1.fq.gz"
R2_trim_L001="${base_trim_path}/${sample}_L001_val_2.fq.gz"
R1_trim_L002="${base_trim_path}/${sample}_L002_val_1.fq.gz"
R2_trim_L002="${base_trim_path}/${sample}_L002_val_2.fq.gz"

# Concatenate the R1 and R2 files from both lanes
cat $R1_trim_L001 $R1_trim_L002 > ${base_trim_path}/${sample}_R1_combined.fq.gz
cat $R2_trim_L001 $R2_trim_L002 > ${base_trim_path}/${sample}_R2_combined.fq.gz

# Run Bowtie2 alignment
bowtie2 --phred33 --mm --maxins 500 --very-sensitive --threads 2 -x /blue/zhangw/hkates/reference_data/homo_sapiens.grch38.v104 \
-1 ${base_trim_path}/${sample}_R1_combined.fq.gz -2 ${base_trim_path}/${sample}_R2_combined.fq.gz \
-S ../3_bowtie/${sample}.sam

