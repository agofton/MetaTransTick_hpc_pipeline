#!/bin/bash

#SBATCH --job-name=k2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128GB
#SBATCH --time=48:00:00
#SBATCH --output=/datasets/work/hb-austicks/work/Project_Phoenix/data/pipeline_test_1/PP-K-p10/logs/k2_%A.slurm.log

date

module load kraken2/2.0.9b

db='/datasets/work/hb-austicks/work/Project_Phoenix/data/NGS/kraken2/maxikraken2_1903_140GB_March2019'
r1='/datasets/work/hb-austicks/work/Project_Phoenix/data/pipeline_test_1/PP-K-p10/PP-K-p10_R1.QC.paired.fastq.gz'
r2='/datasets/work/hb-austicks/work/Project_Phoenix/data/pipeline_test_1/PP-K-p10/PP-K-p10_R2.QC.paired.fastq.gz'
clasOut='/datasets/work/hb-austicks/work/Project_Phoenix/data/pipeline_test_1/PP-K-p10/PP-K-p10_R#.k2Class' 		#=1 or 2
unClasOut='/datasets/work/hb-austicks/work/Project_Phoenix/data/pipeline_test_1/PP-K-p10/PP-K-p10_R#.k2Unclass'		#=1 or 2
out='/datasets/work/hb-austicks/work/Project_Phoenix/data/pipeline_test_1/PP-K-p10/PP-K-p10_k2_out'
con=0.7

kraken2 --paired \
		--classified-out ${clasOut} \
		--unclassified-out ${unClasOut} \
		--db ${db} \
		--threads 20  \
		--output ${out} \
		--confidence ${con} \
		--use-names \
		${r1} ${r2}   
