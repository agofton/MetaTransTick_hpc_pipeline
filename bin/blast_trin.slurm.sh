#!/bin/bash
date

flowControl() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	else
		echo $2; date; echo "";
	fi
}

module load bioref
module load blast+/2.9.0

database="/data/bioref/blast/ncbi/nt"
in=../SAMPLEID.trinity.fasta
archive=../SAMPLEID.trinity.ASN.1
outfmt0=../SAMPLEID.trinity.blast
outfmt6=../SAMPLEID.trinity.blast.txt
topHits=../SAMPLEID.trinity.blast.topHits.txt

# Blast and write blast archive .ASN.1
blastn -task megablast -query ${in} -db ${database} -out ${archive} -strand both -num_threads ${SLURM_CPUS_PER_TASK} -outfmt 11 \
	-max_target_seqs 100 -max_hsps 5 -evalue 1e-10

flowControl "Megablast failed: ${in}" "Megablast finished sucessfully: ${in}"

# Write pairwise output
blast_formatter -archive ${archive} -outfmt 0 -out ${outfmt0}

flowControl "Writing pairwise output failed: ${outfmt0}" "Pairwise output written: ${outfmt0}"

# Write tabular output
blast_formatter -archive ${archive} -out ${outfmt6} \
	-outfmt '6 qseqid ssciname saccver stitle staxid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

flowControl "Writing tabular output failed: ${outfmt6}" "Tabular output written: ${outfmt6}"

# Writing top hits only
cat ${outfmt6} | awk '!a[$1]++' > ${topHits}

# counting number of hits - this could also be done by counting number of uniq seqIDs in $outfmt6
numnohits=$(grep -c "No hits found" ${outfmt0})
numqueries=$(grep -c "^>" ${in})
num_hits=$[numqueries-numnohits]
echo "Number of reads with blast hits:"
echo "${outfmt0}: ${num_hits}"
numHits=$(awk '{print $1}' ${outfmt6} | sort -u | wc -l)
echo "${outfmt6}: ${numHits}"

date
