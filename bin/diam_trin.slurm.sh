#!/bin/bash
date

flowControl() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1
	else
		echo $2
	fi
}

module load bioref
module load diamond/0.9.28

database=/data/bioref/diamond_db/nr_191112_v0928.dmnd #updated Nov 12 2019
input=../SAMPLEID.trinity.fasta
output=../SAMPLEID.trinity.dmnd

diamond blastx --threads ${SLURM_CPUS_PER_TASK} --db ${database} --out ${output} \
	--outfmt 0 --query ${input} --strand both --max-target-seqs 100 --evalue 0.000000001 \
	--sensitive --min-orf 80
flowControli "Diamond failed: ${input}" "Diamond finished sucessfully: ${input}" 

# counting number of hits
numnohits=$(grep -c "No hits found" ${output})
numqueries=$(grep -c "^>" ${input})
num_hits=$[numqueries-numnohits]
echo "Number of reads with diamond hits"
echo "${out}: ${num_hits}"

date
