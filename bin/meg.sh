#!/bin/bash

# script to be run on VM not pearcy - doesnt work becuse of some java error
hm='./meg.sh -i [input blast file] -r [reads.fasta] -o [output.rma6] -l [lca.txt] -t ["nucl" or "prot"] -h [print help]'

# default lca params
maxMatchesPerRead='5000'
minScore='100' 					# bit score
maxExpected='0.0000000001' 		# e-value (1E-10)
minPercentIdentity='50'
topPercent='10'
minSupportPercent='0' 			# 0 = off
minPercentReadCover='0'
minPercentReferenceCover='0'
lcaAlgorithm='naive' 			# "naive", "weighted", "longReads"
readAssignmentMode='readCount' 	

# Command line arguments
while getopts hi:r:o:t:l: option
do
	case "${option}"
		in
	        h) echo ${hm}
	           echo ""
	           exit;;
	        i) blast_in=${OPTARG};;
	        r) reads_in=${OPTARG};;
			o) rma_out=${OPTARG};;
			l) lca_out=${OPTARG};;
			t) typ=${OPTARG};;
	        :) printf "missing argument for  -%s\n" "$OPTARG" >&2
	           echo "$usage" >&2
	           exit 1;;
	   	   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
	           echo "$usage" >&2
	           exit 1;;
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

#b2r='~/megan/tools/blast2rma'
#
#echo ""
#date
#echo ""
#
#~/megan/tools/blast2rma -i ${blast_in} -r ${reads_in} -o ${rma_out} \
#	-f BlastText -bm Blast${np} \
#	-c -m ${maxMatchesPerRead} -ms ${minScore} \
#	-me ${maxExpected} -mpi ${minPercentIdentity} \
#	-top ${topPercent} -supp ${minSupportPercent} \
#	-mrc ${minPercentReadCover} -mrefc ${minPercentReferenceCover} \
#	-alg ${lcaAlgorithm} -ram ${readAssignmentMode} \
#	-a2t ${taxMap} -v


date
echo ""

~/megan/tools/blast2lca -i ${blast_in} -f BlastText -m Blast${np} -o ${lca_out} \
	  -sr -oro -ms ${minScore} -me ${maxExpected} -top ${topPercent} \
	 -mid ${minPercentIdentity} -tn -a2t ${taxMap} -v

echo ""
date
