#!/bin/bash

helpMessage() {
	echo "$0 usage: -i <input .blast or .diamond file> -r <input reads.fasta> -t <"nucl" or "prot"> -o <output.rma6> -h print this help message";
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

# lca params
maxMatchesPerRead='5000'
minScore='200' 					# bit score
maxExpected='0.0000000001' 		# e-value (1E-10)
minPercentIdentity='70'
topPercent='10'
#minSupportPercent='0' 			# 0 = off
minSupport='2' 					# no singletons
minPercentReadCover='0'
#minPercentReferenceCover='50.0'
lcaAlgorithm='weighted' 		# "naive", "weighted", "longReads"
readAssignmentMode='readCount'

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

# nucl or prot
if [ "${typ}" == "nucl" ]; then
	np='N'
	taxMap='/data/hb-austick/work/Project_Phoenix/data/megan_taxonomy_files/nucl_acc2tax-Jul2019.abin'
elif [ "${typ}" == "prot" ]; then
	np='X'
	taxMap='/data/hb-austick/work/Project_Phoenix/data/megan_taxonomy_files/prot_acc2tax-Jul2019X1.abin'
else
	echo 'ERROR: -t must be "prot" or "nucl"'
fi

date

~/megan/tools/blast2rma \
	-i ${blast_in} \
	-r ${reads_in} \
	-o ${rma_out} \
	-f BlastText \
	-bm Blast${np} \
	-c \
	-m ${maxMatchesPerRead} \
	-ms ${minScore} \
	-me ${maxExpected} \
	-mpi ${minPercentIdentity} \
	-top ${topPercent} \
	-sup ${minSupport} \
	-mrc ${minPercentReadCover} \
	-alg ${lcaAlgorithm} \
	-ram ${readAssignmentMode} \
	-a2t ${taxMap} \
	-v

errorExit \
	"blast2rma failed ..." \
	"blast2rma complete ..."
date
#~/megan/tools/blast2lca \
#	-i ${blast_in} \
#	-f BlastText \
#	-m Blast${np} \
#	-o ${lca_out} \
#	-sr \
#	-oro \
#	-ms ${minScore} \
#	-me ${maxExpected} \
#	-top ${topPercent} \
#	-mid ${minPercentIdentity} \
#	-tn \
#	-a2t ${taxMap} \
#	-v
#
#errorExit \
#	"blast2lca failed ..." \
#	"blast2lca complete ..."
