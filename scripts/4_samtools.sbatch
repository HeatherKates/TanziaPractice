#!/bin/bash
#SBATCH --job-name=samtools.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb
#SBATCH --qos=timgarrett
#SBATCH --account=timgarrett
#SBATCH -t 24:00:00
#SBATCH --output=../logs/samtools.out
#SBATCH --error=../logs/samtools.err

cd ../
module load samtools/1.20
# Step 1: Convert the SAM file to BAM
samtools view -h -b 3_bowtie/MSH2KO-1.sam > 4a_samtools/aligned.bam

# Step 2: Apply filters and output filtered BAM
samtools view -h -b -F 1804 -f 2 -q 30 -@ 10 4a_samtools/aligned.bam > 4a_samtools/aligned.filt.bam 

# Step 4: Sort the filtered BAM by read name. This output file will be used for genrich peak calling
samtools sort -n -O BAM -@ 2 4a_samtools/aligned.filt.bam > 4a_samtools/aligned.sorted.bam

# Flagstat analysis
# Run flagstat on aligned.bam to get alignment statistics
samtools flagstat -@ 2 4a_samtools/aligned.bam > 4b_flagstat/aligned.bam.flagstat.log

