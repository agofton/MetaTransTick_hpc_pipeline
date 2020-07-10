#!/bin/bash

helpMessage() {
	echo "$0 usage: -b <blast file in outfmt6 (tab separated; taxon sciname as column 2)> -o <output file> -h [show this message]>"
	exit 1
}

errorExit () {
	if [[ $? -ne 0 ]]; then
		echo $1
		date
		exit 1
	fi
}

# Default values
topHits=/scratch1/gof005/topHits.tmp
taxList=/scratch1/gof005/taxlist.tmp
output=./spp_summary.txt

# Command line args
while getopts hb:o: option
do
	case "${option}"
		in
			h) helpMessage;;
			b) blastFile;;
			o) output;;
	esac
done
shift $((OPTIND - 1))
#######################################
date

# Write top hits only
awk '!a[$1]++' ${blastFile} > ${topHits}

# Print all unique taxon names
awk -F '\t' '{print $2}' ${topHits} | sort | sort -u > ${taxList}

# Start output
echo Taxon$'\t'#Reads$'\t'%Reads > ${output}
echo '--------------------------------------------' >> ${output}

# For each unique taxon name calculate the number and % of sequences assigned to that taxon and print out file
while read -r line
do

	total=$(wc -l ${topHits} | awk '{print $1}')
	onePct=$(awk "BEGIN {print ${total}/100}")
	num=$(grep -c "${line}" $topHits)
	pct=$(awk "BEGIN {print ${num}/${onePct}}")
	echo ${line}$'\t'${num}$'\t'${pct} >> ${output}

done<${taxList}

# cleanup
rm -f ${topHits}
rm -f ${taxList}

date
