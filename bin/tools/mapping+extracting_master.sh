#!/bin/bash

# Alexander W. Gofton, CSIRO, Jan 2021 alexander.gofton@csiro.au; alexander.gofton@gmail.com

# Error control 
errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	else
		echo $2
	fi
}

# Default parameters
THREADS=8

# Command line args & help
hmessage=""
usage="Usage: $(basename "$0")"

while getopts ht:l:f:r:p:s:o: option; do
        case "${option}" in
        	h) echo "$hmessage"
              echo "$usage"
              exit;;  
        	t) TRIN_FASTA=$OPTARG;;
			l) LCA=$OPTARG;;
			f) R1=$OPTARG;;
			r) R2=$OPTARG;;
			p) THREADS=$OPTARG;;
			s) SAMPLE_ID=$OPTARG;;
			o) OUT_DIR=$OPTARG;;
        esac
done
shift $((OPTIND - 1))

# Set output
mkdir ${OUT_DIR}

# 1. Number of contigs
N_CONTIGS=$(grep -c "^>" ${TRIN_FASTA})
errorExit "Counting all contigs complete :)" "Counting all contigs failed :("

# 2. Number of MEGAN assigned contigs
N_ASS_CONTIGS=$(awk -F "\t" '{print $1}' ${LCA} | tail -n +2 | wc -l)
errorExit "Counting MEGAN assigned contigs complete :)" "Counting MEGAN assigned contigs failed :("

# 3. Map all reads to all contigs with bowtie2
INDEX=${OUT_DIR}/${SAMPLE_ID}_trinity_bt2_index
SAM=${OUT_DIR}/${SAMPLE_ID}_trinity_bt2.sam
TRIN_FLAGSTAT=${OUT_DIR}/${SAMPLE_ID}_trinity_bt2_flagstat.txt

~/bowtie2-2.4.2-linux-x86_64/bowtie2-build --threads ${THREADS} --quiet ${TRIN_FASTA} ${INDEX}
errorExit "bowtie2-build failed :(" "bowtie-build complete :)"

~/bowtie2-2.4.2-linux-x86_64/bowtie2 -p ${THREADS} -x ${INDEX_OUT} -1 ${R1} -2 ${R2} -S ${SAM} --fast
errorExit "bowtie2 failed :(" "bowtie2 complete :)"

samtools flagstat -@ ${THREADS} ${SAM} > ${TRIN_FLAGSTAT}
errorExit "samtools flagstat failed :(" "samtools flagstat complete :)"

# Extract tax contigs
extract_and_subset_sam() {
	SEQIDS=${OUT_DIR}/${SAMPLE_ID}_${1}_contigs.seqIDs.txt
	TAX_FASTA=${OUT_DIR}/${SAMPLE_ID}_${1}_contigs.fasta
	SUB_SAM=${OUT_DIR}/${SAMPLE_ID}_${1}_contigs.sam

	echo "Extracting $1 contigs..."

	awk -F "\t" '{print $1, $2, $3}' ${LCA} | grep ${1} | awk -F " " '{print $1}' > ${SEQIDS}
	errorExit "Finding $1 contigs in ${SAMPLE_ID} failed :(" "Finding $1 contigs in ${SAMPLE_ID} complete :)"

	usearch9.2_linux64 -fastx_getseqs ${TRIN_FASTA} -labels ${SEQIDS} -fastaout ${TAX_FASTA}
	errorExit "Extracting $1 contigs in ${SAMPLE_ID} failed :(" "Extracting $1 contigs in ${SAMPLE_ID} complete :)"

	grep -f ${SEQIDS} ${SAM} > ${SUB_SAM}
	errorExit "Subsetting $1 contigs in ${SAM} filed :(" "Subsetting $1 contigs in ${SAM} compelete :)"
}

extract_contigs "Ixodida"
extract_contigs "Viruses"
extract_contigs "Apicomplexa"
extract_contigs "Trypanosoma"
extract_contigs "Anaplasmataceae"
extract_contigs	"Francisella"
extract_contigs "Rickettsia"
extract_contigs "Borreliaceae"


## Calculation of mapping stats (python)
#calc_mapped_reads() {
#	python 3 - <<END
#import pandas as pd
#df = pd.read_csv("${TRIN_FLAGSTAT}", sep = '\s+', header = None, names=["col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11"])
#num_reads = df.iloc[0, 0] - df.iloc[1, 0] - df.iloc[2, 0]
#num_unmapped = df.iloc[0, 0] - df.iloc[4, 0]
#num_mapped = num_reads - num_unmapped 
#result = (num_mapped / num_reads) * 100
#rounded_result = round(result, 2)
#print(str(num_mapped) + " of " + str(num_reads) + " mapped.")
#print(str(rounded_result) + "% of reads mapped.")
#print()
#END
#}
#calc_mapped_reads
#errorExit "Calculaing mapping % failed :(" "Calculaing mapping % complete :)"




