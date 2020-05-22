#!/bin/bash

#SBATCH --job-name=QC_SAMPLEID
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB
#SBATCH --output=../logs/QC_SAMPLEID_%A.log
#SBATCH --time=00:10:00

module load fastqc/0.11.8
module load trimmomatic/0.38
module load python/3.6.1
module load bbtools/38.37

echo ""
date
echo ""

# 1. fastqc on raw data
in1=SCRATCHDIR/SAMPLEID_R1.fastq.gz
in2=SCRATCHDIR/SAMPLEID_R2.fastq.gz

mkdir SCRATCHDIR/fastqc_raw

fastqc \
	--outdir SCRATCHDIR/fastqc_raw \
	--format fastq \
	--threads ${SLURM_CPUS_PER_TASK} ${in1} ${in2}

# flow control
if [ $? -ne 0 ]; then
	echo ""; echo "FASTQC failed: ${in1}, ${in2}" ; date ; exit 1
else
	echo ""; echo "FastQC finished sucessfully: ${in1}, ${in2}" ; date; echo ""
fi

# 2. removing sequencing adapters and distal bases
cut_out1=SCRATCHDIR/SAMPLEID_R1.cutadapt.fastq.gz
cut_out2=SCRATCHDIR/SAMPLEID_R2.cutadapt.fastq.gz

cutadapt \
	-j ${SLURM_CPUS_PER_TASK} \
	-a AGATCGGAAGAG \
	-A AGATCGGAAGAG \
	--no-indels \
	-O 4 \
	-m 30 \
	-o ${cut_out1} \
	-p ${cut_out2} \
	${in1} ${in2}

# flow control
if [ $? -ne 0 ]; then
	echo ""; echo "CUTADAPT failed: ${in1}, ${in2}" ; date ; exit 1
else
	echo ""; echo "CUTADAPT finished sucessfully: ${in1}, ${in2}" ; date; echo ""
fi

# 3. trim low quality reads with a sliding window and remove short seqs < 30 bp
trim_sum=SCRATCHDIR/SAMPLEID.trimmomatic.summary
trim_paired_R1=SCRATCHDIR/SAMPLEID_R1.QC.paired.fastq.gz
trim_unpaired_R1=SCRATCHDIR/SAMPLEID_R1.QC.unpaired.fastq.gz
trim_paired_R2=SCRATCHDIR/SAMPLEID_R2.QC.paired.fastq.gz
trim_unpaired_R2=SCRATCHDIR/SAMPLEID_R2.QC.unpaired.fastq.gz

trimmomatic PE \
	-threads ${SLURM_CPUS_PER_TASK} \
	-summary ${trim_sum} \
	${cut_out1} ${cut_out2} \
	${trim_paired_R1} ${trim_unpaired_R1} ${trim_paired_R2} ${trim_unpaired_R2} \
	SLIDINGWINDOW:5:15 \
	MINLEN:30

# flow control
if [ $? -ne 0 ]; then
	echo ""; echo "TRIMMOMATIC failed: ${cut_out1}, ${cut_out2}" ; date ; exit 1
else
	echo ""; echo "TRIMMOMATIC finished sucessfully: ${cut_out1}, ${cut_out2}" ; date; echo ""
fi

# fastqc on QC data
mkdir SCRATCHDIR/fastqc_QC

fastqc \
	--outdir SCRATCHDIR/fastqc_QC \
	--format fastq \
	--threads ${SLURM_CPUS_PER_TASK} \
	${trim_paired_R1} ${trim_paired_R2} ${trim_unpaired_R1} ${trim_unpaired_R2}

# flow control
if [ $? -ne 0 ]; then
	echo ""; echo "FASTQC failed: ${trim_paired_R1}, ${trim_paired_R2}, ${trim_unpaired_R1}, ${trim_unpaired_R2}" ; date ; exit 1
else
	echo ""; echo "FastQC finished sucessfully: ${trim_paired_R1}, ${trim_paired_R2}, ${trim_unpaired_R1}, ${trim_unpaired_R2}" ; date; echo ""
fi

# Reporting number of reads at each step
count_in1=$(zcat ${in1} | grep -c "^+")
count_in2=$(zcat ${in2} | grep -c "^+")
count_cut_out1=$(zcat ${cut_out1} | grep -c "^+")
count_cut_out2=$(zcat ${cut_out2} | grep -c "^+")
count_trim_paired_R1=$(zcat ${trim_paired_R1} | grep -c "^+")
count_trim_paired_R2=$(zcat ${trim_paired_R2} | grep -c "^+")
count_trim_unpaired_R1=$(zcat ${trim_unpaired_R1} | grep -c "^+")
count_trim_unpaired_R2=$(zcat ${trim_unpaired_R2} | grep -c "^+")

echo "# reads in inputs and outputs."
echo "${in1}: ${count_in1}"
echo "${in2}: ${count_in2}"
echo "${cut_out1}: ${count_cut_out1}"
echo "${cut_out2}: ${count_cut_out2}"
echo "${trim_paired_R1}: ${count_trim_paired_R1}"
echo "${trim_paired_R2}: ${count_trim_paired_R2}"
echo "${trim_unpaired_R1}: ${count_trim_unpaired_R1}"
echo "${trim_unpaired_R2}: ${count_trim_unpaired_R2}"

# moving output back to hb-auticks
cp -r SCRATCHDIR/fastqc_raw OUTDIR/fastqc_raw
cp -r SCRATCHDIR/fastqc_QC OUTDIR/fastqc_QC
cp ${cut_out1} OUTDIR/
cp ${cut_out2} OUTDIR/
cp ${trim_sum} OUTDIR/
cp ${trim_paired_R1} OUTDIR/
cp ${trim_paired_R2} OUTDIR/
cp ${trim_unpaired_R1} OUTDIR/
cp ${trim_unpaired_R2} OUTDIR/

# print end date
echo ""
date
echo ""
