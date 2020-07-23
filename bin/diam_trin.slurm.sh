#!/bin/bash
date

source slurmParams.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1
	else
		echo $2
	fi
}

module load bioref
module load diamond/0.9.28

mkdir ${diamDir}

diamond blastx \
	--threads ${SLURM_CPUS_PER_TASK} \
	--db ${database} \
	--out ${diamOut} \
	--outfmt 0 \
	--query ${trinIn} \
	--strand both \
	--max-target-seqs 100 \
	--evalue 0.000000001 \
	--sensitive \
	--min-orf 80

errorExit \
	"diamond failed: ${trinIn}" \
	"diamond finished sucessfully: ${trinIn}"

diamond view \
	--daa ${diamOut} \
	--out ${diamTabOut} \
	--outfmt '6 qseqid sscinames sseqid stitle staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore'

errorExit \
	"diamond view failed: ${diamTabOut}" \
	"diamond view complete: ${diamTabOut}"

# counting number of hits
numnohits=$(grep -c "No hits found" ${diamOut})
numqueries=$(grep -c "^>" ${trinIn})
num_hits=$[numqueries-numnohits]
echo "Number of reads with diamond hits"
echo "${out}: ${num_hits}"

date
