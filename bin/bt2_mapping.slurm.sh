#!/bin/bash
date
source slurmParams.txt
source script_vars.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	fi
}

# Load software
module load bowtie/2.3.4
module load samtools/1.10.0
module load python/3.7.2

# Make output directory
mkdir -p ${bt2OutDir}

# Build bowtie2 index of whole transcriptome assembly
bowtie2-build --threads ${SLURM_CPUS_PER_TASK} \
			  ${trinFasta} ${bt2index}

errorExit "Bowtie2 indexing of trinity assembly failed ..."
echo "bowtie2 indexing complete. "

# Map all QCd reads to whole transcriptome with bowtie2
bowtie2 -p ${SLURM_CPUS_PER_TASK} \
		-x ${bt2index} \
		-1 ${derepOut1} \
		-2 ${derepOut2} \
		-S ${bt2_sam}
		--no-unal 2> ${bt2AlignStats}

errorExit "Bowtie2 failed :("
echo "Bowtie2 complete."

samtools view -@ ${SLURM_CPUS_PER_TASK} -Sb ${bt2_sam} |
	samtools sort -@ ${SLURM_CPUS_PER_TASK} > ${bt2SortedBam} && \
		rm -f ${bt2_sam}

# Generating mapping stats from each contig
samtools coverage -@ ${SLURM_CPUS_PER_TASK} ${bt2SortedBam} > ${bt2TransSum}

errorExit "Samtools coverage failed!"
echo "Samtools coverage complete. "

samtools flagstat -@ ${SLURM_CPUS_PER_TASK} ${bt2SortedBam} > ${flagstat}

errorExit "Samtools flagstat failed!"
echo "Samtools flagstat complete!"

date
