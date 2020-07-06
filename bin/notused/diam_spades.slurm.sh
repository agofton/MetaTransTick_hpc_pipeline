#!/bin/bash

#SBATCH --job-name=DiSp_SAMPLEID
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
#SBATCH --output=../logs/diamond_SAMPLEID_sapdes_%A.log
#SBATCH --mem=50GB

module load bioref
module load diamond/0.9.28

database=/data/bioref/diamond_db/nr_191112_v0928.dmnd #updated Nov 12 2019

input=SCRATCHDIR/SAMPLEID.spades.fasta
output=SCRATCHDIR/SAMPLEID.spades.dmnd

echo ""
date
echo ""

diamond blastx \
	--threads ${SLURM_CPUS_PER_TASK} \
	--db ${database} \
	--out ${output} \
	--outfmt 0 \
	--query ${input} \
	--strand both \
	--max-target-seqs 100 \
	--evalue 0.000000001 \
	--sensitive \
	--min-orf 80

# flow control
if [ $? -eq 0 ]; then
	echo "Diamond finished sucessfully: ${input}"; date
else
	echo "Diamond failed: ${input}"; date; exit 1
fi

# counting number of hits
numnohits=$(grep -c "No hits found" ${output})
numqueries=$(grep -c "^>" ${input})
num_hits=$[numqueries-numnohits]

echo "Number of reads with diamond hits"
echo "${out}: ${num_hits}"

# copying output back to hb-austicks
cp ${output} OUTDIR/

echo ""
date
echo ""
