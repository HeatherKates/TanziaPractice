#!/bin/bash
#SBATCH --job-name=samtools.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb
#SBATCH --qos=zhangw
#SBATCH --account=zhangw
#SBATCH -t 24:00:00
#SBATCH --output=logs/samtools.out
#SBATCH --error=logs/samtools.err

cd ../
module load samtools/1.20
# Step 1: Convert the SAM file to BAM
samtools view -h -b bowtie/MSH2R4-1.sam > samtools/aligned.bam

# Step 2: Apply filters and output filtered BAM
samtools view -h -b -F 1804 -f 2 -q 30 -@ 10 samtools/aligned.bam > samtools/aligned_filtered.bam 

#Step 3: Remove mt contam (is this needed?)
samtools view -h -@ 15  samtools/aligned_filtered.bam | grep -v "chrM" | samtools view -h -b -@ 15 - > samtools/aligned_filt_noMT.bam

# Step 4: Sort the filtered BAM by read name
samtools sort -n -O BAM -@ 2 samtools/aligned_filt_noMT.bam > samtools/aligned_filt_noMT.sorted.bam

# Step 5: Fix mate information in the BAM file
samtools fixmate -m -@ 2 samtools/aligned_filt_noMT.sorted.bam samtools/aligned_filt_noMT_fixmate.bam

# Step 6: Sort the fixed BAM by coordinate
samtools sort -O BAM -@ 2 samtools/aligned_filt_noMT_fixmate.bam > samtools/aligned_filtered_sorted.bam

# Step 7: Mark and remove duplicates
samtools markdup -r -S -@ 2 samtools/aligned_filtered_sorted.bam samtools/aligned_filtered_sorted_duprmv.bam

# Step 8: re-index the bam file. Indexing helps in quickly accessing the reads.
samtools index samtools/aligned_filtered_sorted_duprmv.bam samtools/aligned_filtered_sorted_duprmv.bam.bam.bai


# Create a folder for your analysis
mkdir -p flagstat

# Flagstat analysis

# Run flagstat on aligned.bam to get alignment statistics
samtools flagstat -@ 2 samtools/aligned.bam > flagstat/aligned.bam.flagstat.log

# Run flagstat on the final BAM file after filtering, sorting, and duplicate removal
samtools flagstat -@ 2 samtools/aligned_filtered_sorted_duprmv.bam > flagstat/aligned_filtered_sorted_duprmv.bam.flagstat.log

