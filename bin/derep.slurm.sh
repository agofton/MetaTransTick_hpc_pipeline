#!/bin/bash
date

flowControl() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	else
		echo $2
	fi
}

module load bbtools/38.37
module load fastqc/0.11.8

# input & output variables
in1=../SAMPLEID_R1.QC.paired.fastq.gz
in2=../SAMPLEID_R2.QC.paired.fastq.gz
int=../SAMPLEID.derep.interleaved.fastq.gz
out1=../SAMPLEID_R1.derep.fastq.gz
out2=../SAMPLEID_R2.derep.fastq.gz

# removes only optical duplicates when both R1 & R2 are identical. output is interleaved
dedupe.sh in1=${in1} in2=${in2} out=${int} ac=f threads=${SLURM_CPUS_PER_TASK}
flowControl "DEDUPE failed: ${in1}, ${in2}" "DEDUPE finished sucessfully: ${in1}, ${in2}"

# converts intleaved to paired
reformat.sh in=${int} out1=${out1} out2=${out2}
flowControl "REFORMAT failed: ${int}" "REFORMAT finished sucessfully: ${int}"

# cleanup
rm -f ${int}

# fastqc on derep files
mkdir ../fastqc_derep
fastqc --outdir ../fastqc_derep --format fastq --threads ${SLURM_CPUS_PER_TASK} \
	   ${out1} ${out2}
flowControl "FASTQC failed: ${out1}, ${out2}" "FastQC finished sucessfully: ${out1}, ${out2}"

# Reporting number of reads at each step
countIn1=$(zcat ${in1} | grep -c "^+")
countIn2=$(zcat ${in2} | grep -c "^+")
countOut1=$(zcat ${out1} | grep -c "^+")
countOut2=$(zcat ${out2} | grep -c "^+")

echo ""
echo "# Reads in input and output"
echo "${in1}: ${countIn1}"
echo "${in2}: ${countIn2}"
echo "${out1}: ${countOut1}"
echo "${out2}: ${countOut2}"

date
