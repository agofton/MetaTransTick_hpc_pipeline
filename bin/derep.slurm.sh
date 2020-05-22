#!/bin/bash

#SBATCH --job-name=DR_SAMPLEID
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128GB
#SBATCH --output=../logs/derep_SAMPLEID_%A.log
#SBATCH --time=02:00:00

module load bbtools/38.37
module load fastqc/0.11.8

echo ""
date
echo ""

# input & output variables
in1=SCRATCHDIR/SAMPLEID_R1.QC.paired.fastq.gz
in2=SCRATCHDIR/SAMPLEID_R2.QC.paired.fastq.gz
int=SCRATCHDIR/SAMPLEID.derep.interleaved.fastq.gz
out1=SCRATCHDIR/SAMPLEID_R1.derep.fastq.gz
out2=SCRATCHDIR/SAMPLEID_R2.derep.fastq.gz

# removes only optical duplicates when both R1 & R2 are identical. output is interleaved
dedupe.sh \
	in1=${in1} \
	in2=${in2} \
	out1=${out1} \
	out2=${out2} \
	ac=f \
	threads=20

# flow control
if [ $? -ne 0 ]; then
	echo "DEDUPE failed: ${in1}, ${in2}"; date; exit 1
else
	echo "DEDUPE finished sucessfully: ${in1}, ${in2}"; date
fi

# converts intleaved to paired
reformat.sh \
	in=${int} \
	out1=${out1} \
	out2=${out2}

# cleanup
rm -f ${int}

# fastqc on derep files
mkdir SCRATCHDIR/fastqc_derep

fastqc --outdir SCRATCHDIR/fastqc_derep --format fastq --threads ${SLURM_CPUS_PER_TASK} ${out1} ${out2}

# flow control
if [ $? -ne 0 ]; then
	echo "FASTQC failed: ${out1}, ${out2}" ; date ; exit 1
else
	echo "FastQC finished sucessfully: ${out1}, ${out2}" ; date
fi

# Reporting number of reads at each step
count_in1=$(zcat ${in1} | grep -c "^+")
connt_in2=$(zcat ${in2} | grep -c "^+")

count_out1=$(zcat ${out1} | grep -c "^+")
connt_out2=$(zcat ${out2} | grep -c "^+")

echo "# Reads in input and output"
echo "${in1}: ${count_in1}"
echo "${in2}: ${count_in2}"
echo "${out1}: ${count_out1}"
echo "${out2}: ${count_out2}"

# copying outputs back to hb-austicks
cp -r SCRATCHDIR/fastqc_derep OUTDIR/fastqc_raw
cp ${out1} OUTDIR/
cp ${out2} OUTDIR/

# print end date
echo ""
date
echo ""

	

