#!/bin/bash

helpMessage() {
	echo "$0 usage: -i <input .blast or .diamond file>  -t <"nucl" or "prot"> -o <output.lca.txt> -f input format <"DAA" or "BlastText"> -h print this help message";
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
while getopts hi:o:t:f: option
do
	case "${option}"
		in
		h) helpMessage;;
	    i) blast_in=${OPTARG};;
		o) lca_out=${OPTARG};;
		t) typ=${OPTARG};;
		f) format=${OPTARG};;
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

# test $format
if [ "${format}" == "BlastText" ]; then
	format='BlastText'
elif [ "${format}" == "DAA" ]; then
	format='DAA'
else
	echo 'ERROR -f <blast format> must be "BlastText" or "DAA". '
fi


date

~/megan/tools/blast2lca \
	-i ${blast_in} \
	-f ${format} \
	-m Blast${np} \
	-o ${lca_out} \
	-sr \
	-oro \
	-ms ${minScore} \
	-me ${maxExpected} \
	-top ${topPercent} \
	-mid ${minPercentIdentity} \
	-tn \
	-a2t ${taxMap} \
	-v

errorExit \
	"blast2lca failed ..." \
	"blast2lca complete ..."
