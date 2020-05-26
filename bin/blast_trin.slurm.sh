#!/bin/bash

#SBATCH --job-name=BlTr_SAMPLEID
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=06:00:00
#SBATCH --output=../logs/blastN_trinity_SAMPLEID_%A.log
#SBATCH --mem=80GB

module load bioref
module load blast+/2.9.0

database="/data/bioref/blast/ncbi/nt"
in=SCRATCHDIR/SAMPLEID.trinity.fasta
out=SCRATCHDIR/SAMPLEID.trinity.blast

echo ""
date
echo ""

blastn \
	-task megablast \
	-query ${in} \
	-db ${database} \
	-out ${out} \
	-strand both \
	-num_threads ${SLURM_CPUS_PER_TASK} \
	-outfmt 0 \
	-num_alignments 100 \
	-num_descriptions 100 \
	-max_hsps 5 \
	-evalue 1e-10

# flow control
if [ $? -eq 0 ]; then
	echo ""; echo "Megablast finished sucessfully: ${in}"; date
else
	echo ""; echo "Megablast failed: ${in}"; date; exit 1
fi

# counting number of hits
numnohits=$(grep -c "No hits found" ${out})
numqueries=$(grep -c "^>" ${in})
num_hits=$[numqueries-numnohits]

echo "Number of reads with blast hits"
echo "${out}: ${num_hits}"

# copying output back to hb-austicks
cp ${out} OUTDIR/

echo ""
date
echo ""
