#!/bin/bash
#SBATCH --job-name=bowtie2.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb
#SBATCH --qos=zhangw
#SBATCH --account=zhangw
#SBATCH -t 24:00:00
#SBATCH --output=../logs/bowtie2.out
#SBATCH --error=../logs/bowtie2.err

module load bowtie2/2.4.5
cd ../
bowtie2 \
--phred33 \
--mm \
--maxins 500 \  # Adjusted based on the expected short insert size
--very-sensitive \
--threads 2 \
-x /blue/zhangw/hkates/chipseq/reference_data/homo_sapiens.grch38.v104 \
-1 2a_trimgalore/MSH2KO-1_val_1.fq.gz \
-2 2a_trimgalore/MSH2KO-1_val_2.fq.gz \
-S 3_bowtie/MSH2KO-1.sam
