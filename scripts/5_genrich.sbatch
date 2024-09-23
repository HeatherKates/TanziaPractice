#!/bin/bash
#SBATCH --job-name=genrich.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 04:00:00
#SBATCH --output=../logs/genrich.out
#SBATCH --error=../logs/genrich.err

module load genrich/0.5
cd ../

# Define the path where sorted BAM files are located
bam_path="4a_samtools"

# Find all the sorted BAM files
bam_files=$(ls $bam_path/*.sorted.bam)

# Separate KO files (controls) and non-KO files (experiments)
controls=$(echo "$bam_files" | grep "KO" | tr '\n' ',')
experiments=$(echo "$bam_files" | grep -v "KO" | tr '\n' ',')

# Remove trailing commas
controls=${controls%,}
experiments=${experiments%,}

# Run Genrich using the found experiment and control files
Genrich -t $experiments -c $controls -o 5_genrich/MSH2.narrowPeak -j -y -r -v -e MT 
