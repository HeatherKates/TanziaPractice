#!/bin/bash
#SBATCH --job-name=bowtie2.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20gb
#SBATCH --qos=zhangw
#SBATCH --account=zhangw
#SBATCH -t 24:00:00
#SBATCH --output=bowtie2.out
#SBATCH --error=bowtie2.err

module load bowtie2/2.4.5
#module load samtools/0.1.16

bowtie2 \
--phred33 \
--local \
--mm \
--maxins 500 \
--very-sensitive \
--threads 6 \
-x /blue/zhangw/hkates/chipseq/reference_data/homo_sapiens.grch38.v104 \
-1 231-MSH2R4-1_S10_L001_R1_001_val_1.fq.gz \
-2 231-MSH2R4-1_S10_L001_R2_001_val_2.fq.gz \
-S 231-MSH2R4-1_S10_L001.v2.sam
