#!/bin/bash

# Alexander W. Gofton, CSIRO, Jan 2021 alexander.gofton@csiro.au; alexander.gofton@gmail.com

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	fi
}

hmessage() {
	echo "Script searches for a taxon key word (eg. Bacteria) in a tab_separated lca summary file and extract contigs that match that taxon key word."
}

usage() {
	echo " Usage: $(basename "$0")" 
	echo "-l <lca_sum.txt>" 
	echo "-z <taxon to search>"
	echo "-t <trinity.fasta>"
	echo "-o <output .fasta file>"
	echo "-h [how this message]"
	echo "< > = required, [ ] = optional"
}

# Default params go here
threads=8
usearch=usearch9.2_linux64

### Command line arguments ###
while getopts hl:z:s:t:o: option; do
        case "${option}" in
        	h) hmessage ; usage ; exit 0 ;;
        	l) lca=$OPTARG;;
        	z) tax=$OPTARG;;
        	t) fasta=$OPTARG;;
			o) fastaOut=$OPTARG;;
        esac
done
shift $((OPTIND - 1))

# Find and extract contigs 
echo '###################################################'
echo "Finding ${tax} contigs in ${lca} with grep..."

	seqIDs=$(basename $fasta .fasta).seqIDs.txt
	grep ${tax} ${lca} | awk -F "\t" '{print $1}' > ${seqIDs}
	errorExit "grep failed!" ""

	# Check num contigs
	n_contigs=`wc -l ${seqIDs}`

#	if [[ "${n_contigs}" -eq "0" ]]; then
#		echo "Number of ${tax} contigs: 0, exiting script."
#		rm -f ${seqIDs}
#		exit 0
#	else
		echo "Number of ${tax} contigs: ${n_contigs}."
		echo ""
#	fi
	
echo ""
echo "Extracting ${tax} contigs from ${fasta}..."

	${usearch} -fastx_getseqs ${fasta} -labels ${seqIDs} -fastaout ${fastaOut}
	errorExit "Extracting ${tax} contigs from ${fasta} failed!" "Extracting ${tax} contigs from ${fasta} complete."




