#!/bin/bash
date

flowControl() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	else
		echo $2; date; echo "";
	fi
}

module load fastqc/0.11.8
module load trimmomatic/0.38
module load python/3.6.1
module load bbtools/38.37

# Fastqc on raw data
in1=../SAMPLEID_R1.fastq.gz
in2=../SAMPLEID_R2.fastq.gz

mkdir ../fastqc_raw
fastqc --outdir ../fastqc_raw --format fastq --threads ${SLURM_CPUS_PER_TASK} ${in1} ${in2}

flowControl "FASTQC failed: ${in1}, ${in2}" "FastQC finished sucessfully: ${in1}, ${in2}" 

# Removing sequencing adapters and distal bases
cutOut1=../SAMPLEID_R1.cutadapt.fastq.gz
cutOut2=../SAMPLEID_R2.cutadapt.fastq.gz

cutadapt -j ${SLURM_CPUS_PER_TASK} -a AGATCGGAAGAG -A AGATCGGAAGAG --no-indels -O 4 -m 30 -o ${cutOut1} -p ${cutOut2} ${in1} ${in2}
flowControl "CUTADAPT failed: ${in1}, ${in2}" "CUTADAPT finished sucessfully: ${in1}, ${in2}"

# Trim low quality reads with a sliding window and remove short seqs < 30 bp
trimSum=../SAMPLEID.trimmomatic.summary
trimPEr1=../SAMPLEID_R1.QC.paired.fastq.gz
trimUPr1=../SAMPLEID_R1.QC.unpaired.fastq.gz
trimPEr2=../SAMPLEID_R2.QC.paired.fastq.gz
trimUPr2=../SAMPLEID_R2.QC.unpaired.fastq.gz

trimmomatic PE -threads ${SLURM_CPUS_PER_TASK} -summary ${trimSum} \
	${cutOut1} ${cutOut2} \
	${trimPEr1} ${trimUPr1} ${trimPEr2} ${trimUPr2} \
	SLIDINGWINDOW:5:15 MINLEN:30
flowControl "TRIMMOMATIC failed: ${cutOut1}, ${cutOut2}" "TRIMMOMATIC finished sucessfully: ${cutOut1}, ${cutOut2}"

# Fastqc on QC data
mkdir ../fastqc_QC
fastqc --outdir ../fastqc_QC --format fastq --threads ${SLURM_CPUS_PER_TASK} \
	${trimPEr1} ${trimPEr2} ${trimUPr1} ${trimUPr2}
flowControl "FASTQC failed: ${trimPEr1}, ${trimPEr2}, ${trimUPr1}, ${trimUPr2}" "FastQC finished sucessfully: ${trimPEr1}, ${trimPEr2}, ${trimUPr1}, ${trimUPr2}"

# Reporting number of reads at each step
countIn1=$(zcat ${in1} | grep -c "^+")
countIn2=$(zcat ${in2} | grep -c "^+")
countCutOut1=$(zcat ${cutOut1} | grep -c "^+")
countCutOut2=$(zcat ${cutOut2} | grep -c "^+")
countTrimPEr1=$(zcat ${trimPEr1} | grep -c "^+")
countTrimPEr2=$(zcat ${trimPEr2} | grep -c "^+")
countTrimUPr1=$(zcat ${trimUPr1} | grep -c "^+")
countTrimUPr2=$(zcat ${trimUPr2} | grep -c "^+")

echo "# reads in inputs and outputs."
echo "${in1}: ${countIn1}"
echo "${in2}: ${countIn2}"
echo "${cutOut1}: ${countCutOut1}"
echo "${cutOut2}: ${countCutOut2}"
echo "${trimPEr1}: ${countTrimPEr1}"
echo "${trimPEr2}: ${countTrimPEr2}"
echo "${trimUPr1}: ${countTrimUPr1}"
echo "${trimUPr2}: ${countTrimUPr2}"

date
