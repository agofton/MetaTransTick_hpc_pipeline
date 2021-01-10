#!/bin/bash

# Alexander W. Gofton, CSIRO, Jan 2021 alexander.gofton@csiro.au; alexander.gofton@gmail.com

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	else
		echo $2
	fi
}

hmessage="Script searches for a taxon key word (eg. Bacteria) in the lca_summary file (output of contig_lca_sum.py), extract contigs that match that taxon key word."
usage="Usage: $(basename "$0") -l lca_sum.txt -z taxon to search -s sample ID -t trinity.fasta -o output directory (output files will be names with -s sample ID)."


# Default params go here
THREADS=8

### Command line arguments ###
while getopts hl:z:s:t:o: option; do
        case "${option}" in
        	h) echo "$hmessage"
              echo "$usage"
              exit;;
        	l) LCA=$OPTARG;;
        	z) TAX=$OPTARG;;
        	s) SAMID=$OPTARG;;
        	t) FASTA=$OPTARG;;
			o) OUT=$OPTARG;;
        esac
done
shift $((OPTIND - 1))

echo ""
mkdir ${OUT}

# Find and extract contigs 
echo "Finding ${TAX} contigs in ${LCA} with grep..."

	RND=${RANDOM}
	TMP=tmp_${RND}_${SAMID}_${TAX}.txt
	awk -F "\t" '{print $1, $2, $3}' ${LCA} > ${TMP}
	errorExit "awk failed!" ""

	SEQIDS=${OUT}/${SAMID}_${TAX}_contigs.seqIDs.txt
	grep ${TAX} ${TMP} | awk -F " " '{print $1}' > ${SEQIDS}
	errorExit "grep failed!" ""

	# Check num contigs
	n_contigs=$(awk -F "\t" '{print $1, $2, $3}' ${LCA} | grep ${TAX} ${TMP} | wc -l)
	echo "Number of ${TAX} contigs: ${n_contigs}."
	echo ""
	
	if [[ ${n_contigs} = "0" ]]; then
		rm -f ${TMP}
		rm -f ${SEQIDS}
		exit 0
		date
	fi

	rm -f ${TMP}
	errorExit "Searching for ${TAX} contigs in ${LCA} failed!" "Searching for ${TAX} contigs in ${LCA} complete."
	
echo ""
echo "Extracting ${TAX} contigs from ${FASTA}..."

	TAXCONTIGS=${OUT}/${SAMID}_${TAX}_contigs.fasta
	usearch9.2_linux64 -fastx_getseqs ${FASTA} -labels ${SEQIDS} -fastaout ${TAXCONTIGS}
	errorExit "Extracting ${TAX} contigs from ${FASTA} failed!" "Extracting ${TAX} contigs from ${FASTA} complete."

date