#!/bin/bash

helpMessage() {
	echo "$0 usage: -i <trinity output directory> -s <SAMPLE_ID> -h [show this message]"
	echo "< > = required args, [ ] = optional args"
	exit 1
}

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	fi
}

# Command line arguments
while getopts hi:s: option
do
	case "${option}"
		in
		h) helpMessage;; 
		i) INPUT_DIR=${OPTARG};;
		s) SAMPLE_ID=${OPTARG};;	
	esac
done
shift $((OPTIND - 1))
######################
TRIN_PART_AGG=~/trinityrnaseq-v2.11.0/util/support_scripts/partitioned_trinity_aggregator.pl
TRIN_GSM=~/trinityrnaseq-v2.11.0/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py

find ${INPUT_DIR}/read_partitions/ -name '*inity.fasta' | ${TRIN_PART_AGG} --token_prefix TRINITY_DN --output_prefix ${INPUT_DIR}/Trinity.tmp
	errorExit "${TRIN_PART_AGG} failed :("

cp ${INPUT_DIR}/Trinity.tmp.fasta ${INPUT_DIR}/Trinity.fasta
	errorExit "mv ${INPUT_DIR}/Trinity.fasta failed :(" 

${TRIN_GSM} --trinity_fasta ${INPUT_DIR}/Trinity.fasta
	errorExit "${TRIN_GSM} failed :("

# copying output files and adding sample IDs
cd ${INPUT_DIR}

cp ${INPUT_DIR}/Trinity.fasta ../${SAMPLE_ID}.trinity.fasta
	errorExit "cp ${SAMPLE_ID}.trinity.fasta failed :("

cp ${INPUT_DIR}/Trinity.SuperTrans.fasta ../${SAMPLE_ID}.trinity.SuperTrans.fasta
	errorExit "cp ${SAMPLE_ID}.trinity.SuperTrans.fasta failed :("

sed -i "s@TRINITY@${SAMPLE_ID}_TRINITY@g" ../${SAMPLE_ID}.trinity.fasta
	errorExit "sed ${SAMPLE_ID}.trinity.SuperTrans.fasta failed :("

sed -i "s@TRINITY@${SAMPLE_ID}_TRINITY@g" ../${SAMPLE_ID}.trinity.SuperTrans.fasta
	errorExit "sed ${SAMPLE_ID}.trinity.SuperTrans.fasta failed :("	

echo ""
echo "Trinity phase 3 complete `date`"
echo ""
