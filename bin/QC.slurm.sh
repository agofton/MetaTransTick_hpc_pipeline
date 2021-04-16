#!/bin/bash
date
source slurmParams.txt
source script_vars.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	fi
}

module load fastqc/0.11.8
module load trimmomatic/0.38
module load python/3.6.1
module load bbtools/38.37

mkdir ${qcDir}

# 1. FastQC on raw data.
mkdir ${qcDir}/fastqc_raw
fastqc \
	--outdir ${qcDir}/fastqc_raw \
	--format fastq \
	--threads ${SLURM_CPUS_PER_TASK} \
	${in1} ${in2}

errorExit "FASTQC failed: ${in1}, ${in2}"
echo "FastQC finished sucessfully: ${in1}, ${in2}"

# 2. Remove sequencing adapters & distal bases, & and reads under 30 nt
cutadapt \
	-j ${SLURM_CPUS_PER_TASK} \
	-a AGATCGGAAGAG \
	-A AGATCGGAAGAG \
	--no-indels -O 4 -m 30 \
	-o ${cutOut1} -p ${cutOut2} \
	${in1} ${in2}

errorExit "CUTADAPT failed: ${in1}, ${in2}"
echo "CUTADAPT finished sucessfully: ${in1}, ${in2}"

# 3. Trim low quality reads with a sliding window and remove short seqs < 30 bp
trimmomatic PE \
	-threads ${SLURM_CPUS_PER_TASK} \
	-summary ${trimSum} \
	${cutOut1} ${cutOut2} \
	${trimPE_R1} ${trimUP_R1} ${trimPE_R2} ${trimUP_P2} \
	SLIDINGWINDOW:5:15 MINLEN:30

errorExit "TRIMMOMATIC failed: ${cutOut1}, ${cutOut2}"
echo "TRIMMOMATIC finished sucessfully: ${cutOut1}, ${cutOut2}"

# 4. Fastqc on QC data.
mkdir ${qcDir}/fastqc_QC
fastqc \
	--outdir ${qcDir}/fastqc_QC \
	--format fastq \
	--threads ${SLURM_CPUS_PER_TASK} \
	${trimPE_R1} ${trimPE_R2} ${trimUP_R1} ${trimUP_P2}

errorExit "FASTQC failed: ${trimPE_R1}, ${trimPE_R2}, ${trimUP_R1}, ${trimUP_P2}"
echo "FastQC finished sucessfully: ${trimPE_R1}, ${trimPE_R2}, ${trimUP_R1}, ${trimUP_P2}"

# Reporting number of reads at each step
count_in1=$(zcat ${in1} | grep -c "^+")
count_in2=$(zcat ${in2} | grep -c "^+")
count_cutOut1=$(zcat ${cutOut1} | grep -c "^+")
count_cutOut2=$(zcat ${cutOut2} | grep -c "^+")
count_trimPE_R1=$(zcat ${trimPE_R1} | grep -c "^+")
count_trimPE_R2=$(zcat ${trimPE_R2} | grep -c "^+")
count_trimUP_R1=$(zcat ${trimUP_R1} | grep -c "^+")
count_trimUP_R2=$(zcat ${trimUP_P2} | grep -c "^+")

echo "# reads in inputs and outputs."
echo "${in1}: ${count_in1}"
echo "${in2}: ${count_in2}"
echo "${cutOut1}: ${count_cutOut1}"
echo "${cutOut2}: ${count_cutOut2}"
echo "${trimPE_R1}: ${count_trimPE_R1}"
echo "${trimPE_R2}: ${count_trimPE_R2}"
echo "${trimUP_R1}: ${count_trimUP_R1}"
echo "${trimUP_P2}: ${count_trimUP_R2}"

date
