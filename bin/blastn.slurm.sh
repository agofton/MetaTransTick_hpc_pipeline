#!/bin/bash
date
source slurmParams.txt
source script_vars.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	fi
}

module load bioref
module load blast+/2.11.0

mkdir ${blastDir}

# Blast and write blast archive .ASN.1
blastn \
	-task megablast \
	-query ${trinFasta} \
	-db ${database_nt} \
	-out ${archive} \
	-strand both \
	-num_threads ${SLURM_CPUS_PER_TASK} \
	-outfmt 11 \
	-max_target_seqs 100 \
	-max_hsps 5 \
	-evalue 1e-10

errorExit "Megablast failed: ${trinIn}"
echo "Megablast finished sucessfully: ${trinIn}"

# Write pairwise output
blast_formatter \
	-archive ${archive} \
	-outfmt 0 \
	-out ${outfmt0}

errorExit "Writing pairwise output failed: ${outfmt0}"
echo "Pairwise output written: ${outfmt0}"

# Write tabular output
blast_formatter \
	-archive ${archive} \
	-out ${outfmt6} \
	-outfmt '6 qseqid ssciname saccver stitle staxid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

errorExit "Writing tabular output failed: ${outfmt6}"
echo "Tabular output written: ${outfmt6}"

# Writing top hits only
cat ${outfmt6} | awk '!a[$1]++' > ${TOPHITS}

# counting number of hits - this could also be done by counting number of uniq seqIDs in $outfmt6
NUM_NO_HITS=$(grep -c "No hits found" ${outfmt0})
NUM_QUERIES=$(grep -c "^>" ${trinFasta})
NUM_HITS=$[NUM_QUERIES-NUM_NO_HITS]
echo "Number of reads with blast hits:"
echo "${outfmt0}: ${NUM_HITS}"
numHits=$(awk '{print $1}' ${outfmt6} | sort -u | wc -l)
echo "${outfmt6}: ${NUM_HITS}"

date
