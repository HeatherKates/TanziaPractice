#!/bin/bash
#SBATCH --job-name=bigwig.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 02:00:00
#SBATCH --output=../logs/bigwig.out
#SBATCH --error=../logs/bigwig.err

module load deeptools/3.5.2
cd ../

#Step 5: Sort ethe filtered BAM by whatever else. This is for visualization
samtools sort 4a_samtools/aligned.sorted.bam > 6_bigwig/aligned_sort_n.bam

# Step 6: re-index the bam file just in case it's needed
samtools index 6_bigwig/aligned.sort_n.bam

bamCoverage --bam  6_bigwig/aligned_sort_n.bam -o 6_bigwig/aligned.sort_n.bw --binSize 10 --normalizeUsing RPKM -p=6

