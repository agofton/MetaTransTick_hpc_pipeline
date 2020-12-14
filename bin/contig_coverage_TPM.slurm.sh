#!/bin/bash
date
source slurmParams.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	fi
}

module load bwa/0.7.17
module load samtools/1.10.0
module load module load rstudio/1.2.1335

# Make output dir
mkdir -p ${conCovOutDir}

# Build bwa index of transcriptome assembly
bwa index -p ${bwaindex} ${trinFasta}

errorExit "BWA indexing failed ..."
echo "BWA indexing complete."

# Map reads using bwa mem
bwa mem \
	-t ${SLURM_CPUS_PER_TASK} \
	${conCovOutDir}/${bwaindex} \
	${derepOut1} ${derepOut2} | samtools sort > ${bamOut}

errorExit "Error mapping reads with bwa mem ..."
echo "bwa mem complete."

# generate stats
samtools coverage ${bamOut} > ${covSum}

errorExit "Error generating coverage summary with samtools coverage ..."
echo "Coverage summary complete."

# do TPM calcs
./tools/TPM_calc.R ${covSum}

errorExit "Error calculating TPM " 
echo "TPM calculation complete."

date
