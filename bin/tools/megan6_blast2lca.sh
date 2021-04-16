#!/bin/bash

helpMessage() {
	echo """$0 usage: 
  			-i <input .blast or .diamond file>  
            -t <nucl or prot> 
            -o <output.lca.txt> 
            -f input format <DAA or BlastText> 
            -h print this help message
            """
	exit 1
}

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	fi
}

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

# test $format
if [[ "${format}" == "BlastText" ]]; then
	format='BlastText'
elif [[ "${format}" == "DAA" ]]; then
	format='DAA'
else
	echo 'ERROR -f <blast format> must be "BlastText" or "DAA". '
fi

blast2lca_exec=~/megan6/tools/blast2lca

${blast2lca_exec} \
	--input ${blast_in} \
	--format ${format} \
	--mode Blast${np} \
	--output ${lca_out} \
	--showranks \
	--officialRanksOnly \
	--minScore 200 \
	--maxExpected 0.0000000001 \
	--topPercent 10 \
	--minPercentIdentity 70 \
	--parseTaxonNames \
	--mapDB ${taxMap}

errorExit "blast2lca failed ..."
