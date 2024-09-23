#!/bin/bash
#SBATCH --job-name=samtools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 02:00:00
#SBATCH --output=../logs/postbigwig.out
#SBATCH --error=../logs/postbigwig.err

module load samtools
cd ../
#sort the bam file by 
samtools sort -@ 8 4a_samtools/aligned.bam -o 6_bigwig/aligned_sort_n.bam
