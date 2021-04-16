#!/bin/bash
date
source slurmParams.txt
source script_vars.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1
	fi
}

module load bioref
module load diamond/0.9.28

mkdir ${diamDir}

diamond blastx \
	--threads ${SLURM_CPUS_PER_TASK} \
	--db ${database_nr} \
	--out ${diamOut} \
	--outfmt 100 \
	--query ${trinFasta} \
	--strand both \
	--max-target-seqs 100 \
	--evalue 0.000000001 \
	--sensitive \
	--min-orf 80

errorExit "diamond failed: ${trinFasta}"
echo "diamond finished sucessfully: ${trinFasta}"

date
