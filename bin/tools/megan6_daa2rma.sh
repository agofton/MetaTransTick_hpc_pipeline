#!/bin/bash

helpMessage() {
	echo """
	$0 usage: 
	-i <input .daa file>  
	-o <output.rma6> 
	-h print this help message
	"""
	exit 1
}

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1
	fi
}

# Command line arguments
while getopts hi:o:t: option
do
	case "${option}"
		in
		h) helpMessage;;
        i) daa_in=${OPTARG};;
		o) rma_out=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

nr_db=/data/hb-austick/work/Project_Phoenix/data/megan_taxonomy_files/megan-map-Jan2021.db

daa2rma_exe=~/megan6/tools/daa2rma

${daa2rma_exe} \
	--in ${daa_in} \
	--out ${rma_out} \
	--useCompression \
	--maxMatchesPerRead 1000 \
	--classify \
	--minScore 100 \
	--maxExpected 0.0000000001 \
	--minPercentIdentity 60 \
	--topPercent 5 \
	--minSupport 2 \
	--lcaAlgorithm weighted \
	--lcaCoveragePercent 80 \
	--readAssignmentMode readCount \
	--mapDB ${nr_db} \
	--parseTaxonNames \
	--threads 8


errorExit "blast2rma failed ..." \


date
