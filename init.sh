#!/bin/bash

# Alexander W. Gofton CSIRO, 2021, alexander.gofton@gmail.com; alexander.gofton@csiro.au

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

helpMessage() {
	printf "${ORANGE}./init.sh usage: \n"
	printf "${ORANGE}-f <sample_x_R1.fastq.gz> raw input R1 file \n"
	printf "${ORANGE}-r <sample_x_R2.fastq.gz> raw input R1 file \n" 
	printf "${ORANGE}-o <output_folder> all output will go here under a new directoy labeled by the sample ID \n"
	printf "${ORANGE}-h [show this message] \n"
	printf "${ORANGE}<> = manditory argument, [] = optional argument \n"
	printf "\n"
	printf "${ORANGE}This script initialises all downstream slurm processes and directories needed for trascriptome assembly and homology-based taxonomic assignment of transcripts. Assumes input files are in the format: sample_xxx_R1.fastq.gz & sample_xxx_R2.fastq.gz, where 'sample_xxx_' can be whatever samples ID you want (The sample ID is defined as anything preceeding '_R1' and '_R2'. This sample ID will be carried through the whole analysis and appended onto the transcript IDs - so make sure it is not too long. Best to use complete paths so that nothing stuffs up! \n"
	printf "${NC}"
	date
	exit 0
}

errorExit() {
	if [[ $? -ne 0 ]]; then
		printf "${RED}${1} /n"
		date
		exit 1
	fi
}

# Command line arguments
while getopts "hf:r:o:" OPTION
do
	case "${OPTION}"
		in
	        h) helpMessage;;
	        f) R1=${OPTARG};;
	        r) R2=${OPTARG};;
			o) OUTPUT_DIR=${OPTARG};;
			/?) printf "${RED}Invalid option: -$OPTARG" 1>&2

	esac
done
shift $((OPTIND - 1))

date

# Checking that R1 and R2 input files exist.
if [ ! -f "${R1}" ]; then
	errorExit "Error: Input R1 file does not exist."
fi
if [ ! -f "${R2}" ]; then
	errorExit "Error: Input R2 file does not exits."
fi

# Checking that R1 & R2 are not the same file.
if [ "${R1}" == "${R2}" ]; then
	errorExit "Input error: -f and -r cannot be the same."
fi

# Checking if OUTPUT_DIR already exists - will not overwrite existing directory.
if [ -d "${OUTPUT_DIR}" ]; then
	errorExit "Error: Output directory already exists. Cannot overwrite."
fi

# Get sampleID from R1 input file (Any string before 1st underscore).
SAMPLE_ID=$(basename ${R1} _R1.fastq.gz)
	printf "${GREEN}SampleID = ${SAMPLE_ID} \n"
	printf "${NC}"

# Setting up direcory structure.
OUT_DIR=${OUTPUT_DIR}/${SAMPLE_ID} 	# All output goes here
SCRIPTS_DIR=${OUT_DIR}/scripts 		# Scripts are stored and launched from here
mkdir -p ${OUT_DIR}
mkdir ${SCRIPTS_DIR}
mkdir ${OUT_DIR}/logs 				# All slurm log files will go here

# Copy raw data to working directory.
cp ${R1} ${OUT_DIR} 
	
	errorExit "Copying ${R1} to ${OUT_DIR} failed."
	printf "${GREEN}${R1} copied to working dir: ${OUT_DIR} \n"

cp ${R2} ${OUT_DIR}
	
	errorExit "Copying ${R2} to ${OUT_DIR} failed."
	printf "${GREEN}${R2} copied to working dir: ${OUT_DIR} \n"

SLURM_SCRIPTS=/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/*.slurm.sh
# Copy analysis scripts to working directory.
for FILE in ${SLURM_SCRIPTS}
do
	NEW_FILE_NAME=$(basename $FILE .slurm.sh).${SAMPLE_ID}.slurm.sh
	cp ${FILE} ${SCRIPTS_DIR}/${NEW_FILE_NAME}
		
		errorExit "Error: copying .slurm scripts to working directory failed."
		printf "${GREEN}${NEW_FILE_NAME} copied to ${SCRIPTS_DIR}. \n"
done

# Copy run.sh to working dir
RUN_SCRIPT='/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/run.sh'
cp ${RUN_SCRIPT} ${SCRIPTS_DIR}/$(basename ${RUN_SCRIPT} .sh).${SAMPLE_ID}.sh

	errorExit "Error: copying run.sh to working directory failed."
	printf "${GREEN}run.sh copied to ${SCRIPTS_DIR} \n"

# Copy slurmParams.txt to working dir and add universal variables
SLURM_PARAMS='/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/slurmParams.txt'
cp ${SLURM_PARAMS} ${SCRIPTS_DIR}/slurmParams.txt
	
	errorExit "Error: copying slurmParams.txt to working directory failed."
	printf "${GREEN}slurmParams.txt copied to ${SCRIPTS_DIR} \n."

sed -i s@SAMPLEID@${SAMPLE_ID}@g ${SCRIPTS_DIR}/slurmParams.txt

	errorExit "Error: writing SAMPLEID variable to slurmParams.txt failed."
	printf "${GREEN}${SAMPLE_ID} written to ${SCRIPTS_DIR}/slurmParams.txt \n"

# Copy tools dir to working direcory.
TOOLS_DIR='/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/tools/*'
mkdir ${SCRIPTS_DIR}/tools
cp ${TOOLS_DIR} ${SCRIPTS_DIR}/tools/
	
	errorExit "Error: copying bin/tools/ to working directory failed."
	printf "${GREEN}bin/tools copies to working directory. \n"

# Final comments to stdout.
printf "${BLUE}All output directories and analysis scripts created... \n"
printf "${BLUE}To launch analysis run: \n"
printf "${NC}cd ${OUT_DIR}/scripts && ./run.${SAMPLE_ID}.sh -m all \n"
date
