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
awk '!a[$1]++' ${blast_in} > ${output}

