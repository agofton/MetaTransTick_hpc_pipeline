#!/bin/bash

# Alexander W. Gofton, CSIRO, Jan 2021 alexander.gofton@csiro.au; alexander.gofton@gmail.com

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	fi
}

HMESSAGE="Script searches for a taxon key word (eg. Bacteria) in the lca_summary file (output of contig_lca_sum.py), extracts contigs that match that taxon key word from both trinity.fasta and bt2.sam. Then performs some basic mapping rate calculations which are printed to the screen."
USAGE="""Usage: $(basename "$0") 
		-l lca_sum.txt input file 
		-z taxon to search 
		-s seqIDs output file
		-f taxon contigs output.fasta 
		-t trinity.fasta input file 
		-m trinity.sam input file 
		-c tax.sam output file 
		-p tax.flagstat.txt outputfile 
		-y trinity.flagstat.txt input file
		-o mapping_stats_out.txt
"""

# Text colors
BLACK='\033[0;01m'
RED='\033[0;31m'
GREEN='\033[0;32m'
ORANGE='\033[0;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
GRAY='\033[0;37m'
NC='\033[0m' # No Color (White)

# Command line arguments
while getopts hl:z:s:f:t:m:c:p:y:o: option; do
    case "${option}" in
            h) printf "${CYAN}${HMESSAGE}"
			   printf "\n"
               printf "${ORANGE}${USAGE}"
               exit;;
            l) LCA=$OPTARG;;
            z) TAX=$OPTARG;;
            s) SEQIDS=$OPTARG;;
            f) CONTIGS=$OPTARG;;
            t) TRIN_FASTA=$OPTARG;;
			m) SAM=$OPTARG;;
			c) TAX_SAM=$OPTARG;;
			p) TAX_SAM_FLAGSTAT=$OPTARG;;
			y) TRIN_FLAGSTAT=$OPTARG;;
			o) OUT=$OPTARG;;
			/?) printf "${RED}Invalid option: -$OPTARG" 1>&2
    esac
done
shift $((OPTIND - 1))

printf "\n"
printf "${ORANGE}###############"

printf "${BLUE}Finding ${TAX} contigs IDs in ${LCA} ...\n"
awk -F "\t" '{print $1, $2, $3}' ${LCA} | grep ${TAX} | awk -F " " '{print $1}' > ${SEQIDS}
errorExit "Searching for ${TAX} contig IDs in ${LCA} failed!"
printf "${GREEN}Searching for ${TAX} contig IDs in ${LCA} complete... "
printf "${GREEN}${TAX} contig IDs written to ${SEQIDS} \n"

NSEQ=$(cat ${SEQIDS} | wc -l)
if [[ ${NSEQ} -eq 0 ]]; then
	rm -f ${SEQIDS} && printf "${RED}No ${TAX} contigs in ${TRIN_FASTA}!\n" && exit 0
else
	printf "${ORANGE}${TAX} contigs: ${NSEQ} \n"
fi

printf "${BLUE}Extracting ${TAX} contigs from ${TRIN_FASTA}...${NC} \n"
usearch9.2_linux64 -fastx_getseqs ${TRIN_FASTA} -labels ${SEQIDS} -fastaout ${CONTIGS}
printf "\n"
errorExit "Extracting ${TAX} contigs from ${TRIN_FASTA} failed!"
printf "${GREEN}Extracting ${TAX} contigs from ${TRIN_FASTA} complete... "
printf "${GREEN}${TAX} contigs writen to ${CONTIGS}. \n"

printf "${BLUE}Extracting contigs from ${SAM}... \n"
grep -f ${SEQIDS} ${SAM} > ${TAX_SAM}.tmp	
errorExit "Extracting ${TAX} contig alignments from ${SAM} failed!"
printf "${GREEN}${TAX} contigs extracting from ${SAM}."
printf "${GREEN}${TAX} contigs written to ${TAX_SAM}. \n"

printf "${BLUE}Sorting .sam file..."
samtools sort -@ 8 ${TAX_SAM}.tmp > ${TAX_SAM} && rm -f ${TAX_SAM}.tmp
errorExit "Samtools sort failed!"
printf "${GREEN}Samtools sort complete.\n"

printf "${BLUE}Calculating mapping stats with samtools and python ... \n"
samtools flagstat -@ 8 ${TAX_SAM} > ${TAX_SAM_FLAGSTAT}.tmp 2> /dev/null # Std error redirected to /dev/nul to silence verbose output to screen
errorExit "Samtools flagstat failed!"
printf "${GREEN}Samtools flagstat complete..."
printf "${GREEN}${TAX_SAM_FLAGSTAT} written."
printf "${NC}\n"

#Calcs in python
/data/hb-austick/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/tools/tax_mapping_pct.py \
	-t ${TRIN_FASTA} -l ${LCA} -f ${CONTIGS} -x ${TAX_SAM_FLAGSTAT} -y ${TRIN_FLAGSTAT} -z ${TAX} -o ${OUT}

printf "\n"
errorExit "Calculations in python3 failed!"
printf "${GREEN}Calculations in python complete."
date








