#!/bin/bash

helpMessage() {
	echo "$0 usage: -i .txt lca from from megan (read name to taxon path) -h print this help message";
	echo "output = input_file.csv"
	echo ""
	exit 1
}

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date
		exit 1
	fi
}

# Command line arguments
while getopts hi: option
do
	case "${option}"
		in
		h) helpMessage;;
	  	i) input=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

#####################################
output=$(basename ${input} -ex.txt).lca.csv
assigned=$(basename ${input} -ex.txt).lca.assigned.csv

sed 's/,"/,/g' ${input} > ${output}
errorExit "sed #1 failed." ${output}

sed -i 's/;/,/g' ${output}
errorExit "sed 2 failed." ${output}

sed -i 's/,"//g' ${output}
errorExit "sed 3 failed" ${output}

rm -f ${input}
errorExit "rm ${input} failed."

grep -v -e "NCBI$" -e "cellular organisms$" -e "No hits$" -e "Not assigned$ "${output} > ${assigned}
errorExit "grep filtering assigned failed. "

rm -f ${output}
errorExit "rm ${output} failed."


