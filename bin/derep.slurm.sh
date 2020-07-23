#!/bin/bash
date

source slurmParams.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	else
		echo $2; date; echo ""
	fi
}

module load bbtools/38.37
module load fastqc/0.11.8

# removes only optical duplicates when both R1 & R2 are identical.
dedupe.sh \
	in1=${trimPEr1} \
	in2=${trimPEr2} \
	out=${int} \
	ac=f \
	threads=${SLURM_CPUS_PER_TASK}

errorExit \
	"dedupe.sh failed: ${trimPEr1}, ${trimPEr2}" \
	"dedupe.sh finished sucessfully: ${trimPEr1}, ${trimPEr2}"

# converts intleaved to paired
reformat.sh \
	in=${int} \
	out1=${derepOut1} \
	out2=${derepOut2}

errorExit \
	"reformat.sh failed: ${int}" \
	"reformat.sh finished sucessfully: ${int}"

# cleanup
rm -f ${int}

# fastqc on derep files
mkdir ../fastqc_derep

fastqc \
	--outdir ../fastqc_derep \
	--format fastq \
	--threads ${SLURM_CPUS_PER_TASK} \
	${derepOut1} ${derepOut2}

errorExit \
	"FASTQC failed: ${derepOut1}, ${derepOut2}" \
	"FastQC finished sucessfully: ${derepOut1}, ${derepOut2}"

# Reporting number of reads at each step
countIn1=$(zcat ${trimPEr1} | grep -c "^+")
countIn2=$(zcat ${trimPEr2} | grep -c "^+")
countOut1=$(zcat ${derepOut1} | grep -c "^+")
countOut2=$(zcat ${derepOut2} | grep -c "^+")

echo ""
echo "# Reads in input and output"
echo "${trimPEr1}: ${countIn1}"
echo "${trimPEr2}: ${countIn2}"
echo "${derepOut1}: ${countOut1}"
echo "${derepOut2}: ${countOut2}"

date
