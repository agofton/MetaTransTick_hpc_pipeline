#!/bin/bash
date
source slurmParams.txt
source script_vars.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	fi
}

module load bbtools/38.37
module load fastqc/0.11.8

# 1. Removes only optical duplicates when both R1 & R2 are identical.
dedupe.sh \
	in1=${trimPE_R1} \
	in2=${trimPE_R2} \
	out=${interleaved} \
	ac=f \
	threads=${SLURM_CPUS_PER_TASK}

errorExit "dedupe.sh failed: ${trimPE_R1}, ${trimPE_R2}"
echo "dedupe.sh finished sucessfully: ${trimPE_R1}, ${trimPE_R2}"

# 2. Converts interleavedleaved to paired
reformat.sh \
	in=${interleaved} \
	out1=${derepOut1} \
	out2=${derepOut2}

errorExit "reformat.sh failed: ${interleaved}"
echo "reformat.sh finished sucessfully: ${interleaved}"

# cleanup
rm -f ${interleaved}

# 3. fastQC on derep files.
mkdir ${derep_fqc}
fastqc \
	--outdir ${derep_fqc} \
	--format fastq \
	--threads ${SLURM_CPUS_PER_TASK} \
	${derepOut1} ${derepOut2}

errorExit "FASTQC failed: ${derepOut1}, ${derepOut2}"
echo "FastQC finished sucessfully: ${derepOut1}, ${derepOut2}"

# Reporting number of reads at each step
countIn1=$(zcat ${trimPE_R1} | grep -c "^+")
countIn2=$(zcat ${trimPE_R2} | grep -c "^+")
countOut1=$(zcat ${derepOut1} | grep -c "^+")
countOut2=$(zcat ${derepOut2} | grep -c "^+")

echo ""
echo "# Reads in input and output"
echo "${trimPE_R1}: ${countIn1}"
echo "${trimPE_R2}: ${countIn2}"
echo "${derepOut1}: ${countOut1}"
echo "${derepOut2}: ${countOut2}"

date
