#!/bin/bash
#SBATCH --job-name=bigwig.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 02:00:00
#SBATCH --output=../logs/bigwig.out
#SBATCH --error=../logs/bigwig.err

module load deeptools/3.5.2
cd ../
bamCoverage --bam 4a_samtools/4a_samtools/aligned.coordsort.bam.bai -o 6_bigwig/aligned.coordsort.bw --binSize 10 --normalizeUsing RPKM -p=2
