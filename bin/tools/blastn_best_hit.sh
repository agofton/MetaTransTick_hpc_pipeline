#!/bin/bash

helpMessage() {
	echo "$0 usage: -i <input .blast or .diamond file> -r <input reads.fasta> -t <"nucl" or "prot"> -o <output.rma6> -l <output.lca.txt> -h print this help message";
	echo ""
	exit 1
}

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1
	else
		echo $2; date; echo ""
	fi
}


# Command line arguments
while getopts hi:o: option
do
	case "${option}"
		in
		h) helpMessage;;
	  	i) blast_in=${OPTARG};;
		o) output=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

#####

# Sort first by query name, then bitscore, then evalue, then %identity, and extract the best lines for each query.
# bitscore > evalue > %identiy
# k1=query
# k15=bitscore
# k14=evalue
# k6=%identity
# these colums will change if you alter the blast 6 oftfmt format.

sort -k1,1 -k15,15gr -k14,14g -k6,6gr ${blast_in} | sort -u -k 1,1 --merge > ${output}
