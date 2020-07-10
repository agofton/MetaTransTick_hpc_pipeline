#!/bin/bash

#SBATCH --job-name=FqSpD
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80GB
#SBATCH --time=00:30:00
#SBATCH --output=/datasets/work/hb-austicks/work/Project_Phoenix/data/pipeline_test_1/PP-K-p10/species_summary.%A.slurm

date

module load bioref
module load blast+/2.9.0

blastDB='/data/bioref/blast/ncbi/nt'
blastQuery='/datasets/work/hb-austicks/work/Project_Phoenix/data/pipeline_test_1/PP-K-p10/test.fasta'
sumFile='/datasets/work/hb-austicks/work/Project_Phoenix/data/pipeline_test_1/PP-K-p10/speciesSummary.txt'

########################################################
blastOut=/scratch1/gof005/blastOut.tmp
topHits=/scratch1/gof005/topHits.tmp
taxList=/scratch1/gof005/taxlist.tmp

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	fi
}

# Perform blast
blastn -task megablast -query ${blastQuery} -db ${blastDB} -out ${blastOut} -outfmt '6 qseqid ssciname' -max_target_seqs 100 -max_hsps 5 -evalue 1e-10

	errorExit "blastN (megablast) failed ..."

# Write top hits only
awk '!a[$1]++' ${blastOut} > ${topHits}

	errorExit "Writing top hits failed ..."

# Print all unique taxon names
awk -F '\t' '{print $2}' ${topHits} | sort | sort -u > ${taxList}

	errorExit "Extracting uniq taxonIDs failed ..."

# Start output
date > $sumFile
echo Taxon$'\t'#Reads$'\t'%Reads >> $sumFile
echo '--------------------------------------------' >> $sumFile

# For each unique taxon name calculate the number and % of sequences assigned to that taxon and print to sdt out
while read -r line
do
	total=$(wc -l ${topHits} | awk '{print $1}')
	onePct=$(awk "BEGIN {print ${total}/100}")
	num=$(grep -c "${line}" ${topHits})
	pct=$(awk "BEGIN {print ${num}/${onePct}}")
	echo ${line}$'\t'${num}$'\t'${pct} >> $sumFile
done<${taxList}

# cleanup
rm -f ${blastOut}
rm -f ${topHits}
rm -f ${taxList}

