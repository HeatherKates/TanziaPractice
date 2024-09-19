#!/bin/bash
#SBATCH --job-name=genrich.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 02:00:00
#SBATCH --output=../logs/genrich.out
#SBATCH --error=../logs/genrich.err

module load genrich/0.5
cd ../
Genrich -t 4a_samtools/aligned.sorted.bam -o 5_genrich/MSH2KO-1.narrowPeak -j -y -r -v -e MT 
#-c This is the option for control files when they are included
