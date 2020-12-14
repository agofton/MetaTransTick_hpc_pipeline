#!/bin/bash

helpMessage() {
	echo "$0 usage -f file.fasta -b file.blast.outfm6.txt -o output.blast.outfmt6.txt -h [show this message]"
	echo ""
	exit 1
}

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1
	fi
}

# cmd line args
while getopts hf:b:o: option
do
	case "${option}"
		in
		h) helpMessage;;
		f) fasta=${OPTARG};;
		b) blast=${OPTARG};;
		o) output=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

####################
RND=$RANDOM
seqIDs_tmp=seqId.${RND}.tmp

# Extract SeqIDs from fasta
echo "Pulling sequence IDs from ${fasta}..."
grep "^>" ${fasta} | sed 's/>//g' | sort -u > ${seqIDs_tmp}
errorExit "Extracting sequence IDs from fasta file failed"

# Search for these seqIDs in blast file
echo "Extracting blast hits from ${blast}..."
grep -F -w -f ${seqIDs_tmp} ${blast} > ${output}
errorExit "Searching (grep) for sequence IDs in blast file failed"

# Checking results
num_seqIDs=$(cat ${seqIDs_tmp} | wc -l)
echo "Number of seqIDs extracted from ${fasta}: ${num_seqIDs}"
num_blast_reads=$(awk '{print $1}' ${output} | sort -u | wc -l)
echo "Number of seqIDs extracted from ${blast}: ${num_blast_reads}"

# Cleanup
rm -f ${seqIDs_tmp}
errorExit "Failed to deleat tmp file (seqIDs)"

date
