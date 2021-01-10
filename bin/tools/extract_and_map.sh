#!/bin/bash

# Alexander W. Gofton, CSIRO, Jan 2021 alexander.gofton@csiro.au; alexander.gofton@gmail.com

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	else
		echo $2
	fi
}

hmessage="Script searches for a taxon key word (eg. Bacteria) in the lca_summary file (output of contig_lca_sum.py), extract contigs that match that taxon key word, mapps all reads back to those contigs with bwa mem, sorts the bam file, and creasts samtools flagstat.txt file."
usage="Usage: $(basename "$0") -l lca_sum.txt -z taxon to search -s sample ID -t trinity.fasta -f R1.fastx.gz -r R2.fastx.gz -o output directory (output files will be names with -s sample ID) -p threads"


# Default params go here
THREADS=8

### Command line arguments ###
while getopts hl:z:s:t:f:r:o:p: option; do
        case "${option}" in
                h) echo "$hmessage"
                   echo "$usage"
                   exit;;
                l) LCA=$OPTARG;;
                z) TAX=$OPTARG;;
                s) SAMID=$OPTARG;;
                t) FASTA=$OPTARG;;
                f) R1=$OPTARG;;
                r) R2=$OPTARG;;
                o) OUT=$OPTARG;;
				p) THREADS=$OPTARG;;
        esac
done
shift $((OPTIND - 1))

echo ""
mkdir ${OUT}


# Find and extract contigs 
echo "Finding ${TAX} contigs in ${LCA} with grep..."

	RND=${RANDOM}
	TMP=tmp_${RND}_${SAMID}_${TAX}.txt
	awk -F "\t" '{print $1, $2, $3}' ${LCA} > ${TMP}
	errorExit "awk failed!" ""

	SEQIDS=${OUT}/${SAMID}_${TAX}_contigs.seqIDs.txt
	grep ${TAX} ${TMP} | awk -F " " '{print $1}' > ${SEQIDS}
	errorExit "grep failed!" ""

	rm -f ${TMP}
	errorExit "Searching for ${TAX} contigs in ${LCA} failed!" "Searching for ${TAX} contigs in ${LCA} complete."
	
echo ""
echo "Extracting ${TAX} contigs from ${FASTA}..."

	TAXCONTIGS=${OUT}/${SAMID}_${TAX}_contigs.fasta
	usearch9.2_linux64 -fastx_getseqs ${FASTA} -labels ${SEQIDS} -fastaout ${TAXCONTIGS}
	errorExit "Extracting ${TAX} contigs from ${FASTA} failed!" "Extracting ${TAX} contigs from ${FASTA} complete."

echo ""
echo "Mapping reads against ${TAXCONTIGS} with bwa..."
	
	BWAINDEX=${OUT}/${SAMID}_${TAX}_contigs
	bwa index -p ${BWAINDEX} ${TAXCONTIGS}
	errorExit "BWA index failed!" "BWA index complete."

	echo ""
	BAM=${OUT}/${SAMID}_${TAX}_contig.bam
	bwa mem -t ${THREADS} ${BWAINDEX} ${R1} ${R2} | samtools sort > ${BAM}
	errorExit "BWA mem failed!" "BWA mem complete."

	echo ""
	FLAGSTAT=${OUT}/${SAMID}_${TAX}_contig.flagstat.txt
	samtools flagstat -@ ${THREADS} ${BAM} > ${FLAGSTAT} 
	errorExit "Samtools flagstat failed!" "Samtools flagstat complete."

echo "Script complete."
date