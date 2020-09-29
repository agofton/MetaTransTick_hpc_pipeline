#!/bin/bash

helpMessage() {
	echo "$0 usage: -i <trinity output directory> -h [show this message]"
	echo ""
	exit 1
}

errorExit() {
	if [ $? -ne 0 ]; then
		echo $2; date; exit 1
	else
		echo $1; date
	fi
}

# Command line arguments
while getopts hi:s: option
do
	case "${option}"
		in
		h) helpMessage;; 
		i) inputDir=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

####
find ${inputDir}/read_partitions/ -name '*inity.fasta' | \
~/trinityrnaseq-Trinity-v2.8.4/util/support_scripts/partitioned_trinity_aggregator.pl \
		--token_prefix TRINITY_DN \
		--output_prefix ${inputDir}/Trinity.tmp

mv ${inputDir}/Trinity.tmp.fasta ${inputDir}/${inputDir}.Trinity.fasta

sed -i s/>TRINITY/>${sampleID}_TRINITY/g ${inputDir}/${inputDir}.Trinity.fasta

# Counting transcipts
count=$(grep -c "^>" ${inputDir}/${sampleID}.Trinity.fasta)
echo "Number of transcripts."
echo "${inputDir}/${sampleID}.Trinity.fasta: ${count}"
