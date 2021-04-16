#!/bin/bash

# Converts lca (read name to taxon path output, tab separated) 
# from MEGAN6 into a true tsv, with each taxon level in a new column

helpMessage() {
	echo """$0 usage: 
	-i lca.txt - lca from from megan (read name to taxon path, tab separated) 
	-o output.txt 
	-h print this help message
	"""
	exit 0
}

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date
		exit 1
	fi
}

# Command line arguments
while getopts hi:o: option
do
	case "${option}"
		in
		h) helpMessage;;
	  	i) input=${OPTARG};;
		o) output=${OPTARG};;
	esac
done
shift $((OPTIND - 1))
######################
# 1. remove quotes
sed 's/"//g' ${input} > ${output}
	errorExit "sed 1 failed!"
# 2. replace spaces with underscores
sed -i 's/ /_/g' ${output}
	errorExit "sed 2 failed!"
# 3. remove ; at end of lines
sed -i 's/;$//g' ${output}
	errorExit "sed 2 failed!"
# 4. convert internal ; to tab
sed -i "s/;/\\t/g" ${output}
	errorExit "sed 2 failed!"

