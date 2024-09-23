#!/bin/bash
#SBATCH --job-name=samtools.%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --time=24:00:00
#SBATCH --output=../logs/samtools_%A_%a.out
#SBATCH --error=../logs/samtools_%A_%a.err
#SBATCH --array=0-11  # Assuming 12 samples in the list

# Load the samtools module
module load samtools/1.20

# Move to the parent directory
cd ../

# Read the sample name from the sample list
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" MSH2_sample_list.txt)

# Step 1: Convert the SAM file to BAM
samtools view -h -b 3_bowtie/${sample}.sam > 4a_samtools/${sample}.bam

# Step 2: Apply filters and output filtered BAM
samtools view -h -b -F 1804 -f 2 -q 30 -@ 4 4a_samtools/${sample}.bam > 4a_samtools/${sample}.filt.bam 

# Step 3: Sort the filtered BAM by read name (used for peak calling)
samtools sort -n -O BAM -@ 4 4a_samtools/${sample}.filt.bam > 4a_samtools/${sample}.sorted.bam

# Step 4: Run flagstat on aligned BAM to get alignment statistics
mkdir -p 4b_flagstat
samtools flagstat -@ 4 4a_samtools/${sample}.bam > 4b_flagstat/${sample}.bam.flagstat.log

