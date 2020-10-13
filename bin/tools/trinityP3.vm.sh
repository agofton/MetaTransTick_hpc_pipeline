#!/bin/bash

helpMessage() {
	echo "$0 usage: -i <trinity output directory> -s <sampleID> -h [show this message]"
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
		s) sampleID=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

####
find ${inputDir}/read_partitions/ -name '*inity.fasta' | \
~/trinityrnaseq-Trinity-v2.8.4/util/support_scripts/partitioned_trinity_aggregator.pl \
		--token_prefix TRINITY_DN \
		--output_prefix ${inputDir}/Trinity.tmp

errorExit \
	"FIND complete, all trinity assemblies found and coalated :)" \
	"FIND failed :("

cd ${inputDir}
mv Trinity.tmp.fasta ../${sampleID}.trinity.fasta

errorExit \
	"Trinity.fasta moved to its final resting place :)" \
	"MV failed :("

sed -i "s@TRINITY@${sampleID}_TRINITY@g" ../${sampleID}.trinity.fasta

errorExit \
	"SED complete." \
	"SED failed, trinity.fasta does not have sampleIDs"
date
