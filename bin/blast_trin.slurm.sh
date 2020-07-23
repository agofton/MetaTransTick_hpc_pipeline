#!/bin/bash
date

source slurmParams.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	else
		echo $2; date; echo "";
	fi
}

module load bioref
module load blast+/2.9.0

mkdir ${blastDir}

# Blast and write blast archive .ASN.1
blastn \
	-task megablast \
	-query ${trinIn} \
	-db ${database} \
	-out ${archive} \
	-strand both \
	-num_threads ${SLURM_CPUS_PER_TASK} \
	-outfmt 11 \
	-max_target_seqs 100 \
	-max_hsps 5 \
	-evalue 1e-10

errorExit \
	"Megablast failed: ${trinIn}" \
	"Megablast finished sucessfully: ${trinIn}"

# Write pairwise output
blast_formatter \
	-archive ${archive} \
	-outfmt 0 \
	-out ${outfmt0}

errorExit \
	"Writing pairwise output failed: ${outfmt0}" \
	"Pairwise output written: ${outfmt0}"

# Write tabular output
blast_formatter \
	-archive ${archive} \
	-out ${outfmt6} \
	-outfmt '6 qseqid ssciname saccver stitle staxid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

errorExit \
	"Writing tabular output failed: ${outfmt6}" \
	"Tabular output written: ${outfmt6}"

# Writing top hits only
cat ${outfmt6} | awk '!a[$1]++' > ${topHits}

# counting number of hits - this could also be done by counting number of uniq seqIDs in $outfmt6
numnohits=$(grep -c "No hits found" ${outfmt0})
numqueries=$(grep -c "^>" ${trinIn})
num_hits=$[numqueries-numnohits]
echo "Number of reads with blast hits:"
echo "${outfmt0}: ${num_hits}"
numHits=$(awk '{print $1}' ${outfmt6} | sort -u | wc -l)
echo "${outfmt6}: ${numHits}"

date
