#!/bin/bash
#SBATCH --job-name=samtools.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb
#SBATCH --qos=zhangw
#SBATCH --account=zhangw
#SBATCH -t 24:00:00
#SBATCH --output=../logs/samtools.out
#SBATCH --error=../logs/samtools.err

cd ../
module load samtools/1.20
# Step 1: Convert the SAM file to BAM
samtools view -h -b 3_bowtie/MSH2KO-1.sam > 4a_samtools/aligned.bam

# Step 2: Apply filters and output filtered BAM
samtools view -h -b -F 1804 -f 2 -q 30 -@ 10 4a_samtools/aligned.bam > 4a_samtools/aligned_filtered.bam 

#Step 3: Remove mt contam (is this needed?)
samtools view -h -@ 15  4a_samtools/aligned_filtered.bam | grep -v "chrM" | samtools view -h -b -@ 15 - > 4a_samtools/aligned_filt_noMT.bam

# Step 4: Sort the filtered BAM by read name
samtools sort -n -O BAM -@ 2 4a_samtools/aligned_filt_noMT.bam > 4a_samtools/aligned_filt_noMT.sorted.bam

# Step 5: Fix mate information in the BAM file
samtools fixmate -m -@ 2 4a_samtools/aligned_filt_noMT.sorted.bam 4a_samtools/aligned_filt_noMT_fixmate.bam

# Step 6: Sort the fixed BAM by coordinate
samtools sort -O BAM -@ 2 4a_samtools/aligned_filt_noMT_fixmate.bam > 4a_samtools/aligned_filtered_sorted.bam

# Step 7: Mark and remove duplicates
samtools markdup -r -S -@ 2 4a_samtools/aligned_filtered_sorted.bam 4a_samtools/aligned_filtered_sorted_duprmv.bam

# Step 8: re-index the bam file. Indexing helps in quickly accessing the reads.
samtools index 4a_samtools/aligned_filtered_sorted_duprmv.bam 4a_samtools/aligned_filtered_sorted_duprmv.bam.bam.bai

# Flagstat analysis

# Run flagstat on aligned.bam to get alignment statistics
samtools flagstat -@ 2 4a_samtools/aligned.bam > 4b_flagstataligned.bam.flagstat.log

# Run flagstat on the final BAM file after filtering, sorting, and duplicate removal
samtools flagstat -@ 2 4a_samtools/aligned_filtered_sorted_duprmv.bam > 4b_flagstataligned_filtered_sorted_duprmv.bam.flagstat.log

