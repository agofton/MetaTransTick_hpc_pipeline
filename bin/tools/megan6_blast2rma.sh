#!/bin/bash

helpMessage() {
	echo """$0 usage:
		-i <input .blast or .diamond file>
		-r <input reads.fasta> 
		-t <nucl or prot> 
		-o <output.rma6> 
		-h [print this help message]
		< > = required, [ ] = optional
		"""
	exit 0
}

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	fi
}

# Command line arguments
while getopts hi:r:o:t: option
do
	case "${option}"
		in
		h) helpMessage;;
	    i) blast_in=${OPTARG};;
        r) reads_in=${OPTARG};;
		o) rma_out=${OPTARG};;
		t) typ=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

nt_db=/data/hb-austick/work/Project_Phoenix/data/megan_taxonomy_files/megan-nucl-Jan201.db
nr_db=/data/hb-austick/work/Project_Phoenix/data/megan_taxonomy_files/megan-map-Jan2021.db

# nucl or prot
if [[ "${typ}" == "nucl" ]]; then
	np='N'
	taxMap=${nt_db}
elif [[ "${typ}" == "prot" ]]; then
	np='X'
	taxMap=${nr_db}
else
	echo 'ERROR: -t must be "prot" or "nucl"'
fi

megan_exe=~/megan6/tools/blast2rma
bitscore=200

${megan_exe} --in ${blast_in} \
			 --reads ${reads_in} \
			 --format BlastText \
			 --out ${rma_out} \
			 --blastMode Blast${np} \
			 --useCompression \
			 --maxMatchesPerRead 1000 \
			 --minScore ${bitscore} \
			 --maxExpected 0.0000000001 \
			 --minPercentIdentity 70.0 \
			 --topPercent 5 \
			 --minSupport 2 \
			 --lcaAlgorithm weighted \
			 --lcaCoveragePercent 80 \
			 --readAssignmentMode readCount \
			 --mapDB ${taxMap} \
			 --threads 8

errorExit "blast2rma failed ..."
date
