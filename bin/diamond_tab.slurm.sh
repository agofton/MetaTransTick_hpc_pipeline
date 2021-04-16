#!/bin/bash
date
source slurmParams.txt
source script_vars.txt

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	fi
}

module load bioref
module load diamond/0.9.28

mkdir ${diamDir}

diamond blastx \
	--threads ${SLURM_CPUS_PER_TASK} \
	--db ${database_nr} \
	--out ${diamTabOut} \
	--outfmt 6 qseqid sscinames sseqid stitle staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore \
	--query ${trinFasta} \
	--strand both \
	--max-target-seqs 100 \
	--evalue 0.000000001 \
	--sensitive \
	--min-orf 80

errorExit "diamond failed: ${trinFasta}"
echo "diamond finished sucessfully: ${trinFasta}"

# take top hits only
awk '!a[$1]++' ${diamTabOut} > ${diamTopHits}
errorExit "Generating top hits filed ${diamTopHits}"

date
