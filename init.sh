#!/bin/bash
date

helpMessage="
Usage: ./init.sh [OPTIONS]
Description:	This script initialises all downstream slurm processes and directories needed for trascriptome assembly and homology-based 
		taxonomic assignment of transcripts. Assumes input files are in the format: sample_xxx_R1.fastq.gz & sample_xxx_R2.fastq.gz, 
		where 'sample_xxx_' can be whatever samples ID you want (The sample ID is defined as anything preceeding '_R1' and '_R2'. 
		This sample ID will be carried through the whole analysis and appended onto the transcript IDs, so make sure it is not too long.
Options:
-f <sample_x_R1.fastq.gz> raw input R1 file
-r <sample_x_R2.fastq.gz> raw input R1 file
-o <output_folder> all output will go here under a new directoy labeled by the sample ID
-h [show this message]
< > = manditory argument, [  ] = optional argument
"

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1 ; date ; exit 1
	fi
}

# Command line arguments
while getopts "hf:r:o:" OPTION
do
	case "${OPTION}"
		in
	        h) echo "${helpMessage}
		   date
		   exit 0;;
	        f) R1=${OPTARG};;
	        r) R2=${OPTARG};;
		o) outDir=${OPTARG};;
		/?) echo "Invalid option: -$OPTARG" 1>&2
	esac
done
shift $((OPTIND - 1))

### Checking that R1 and R2 input files exist.
if [ ! -f "${R1}" ]; then
	errorExit "Error: Input R1 file does not exist."
elif [ ! -f "${R2}" ]; then
	errorExit "Error: Input R2 file does not exits."
fi

### Checking that R1 & R2 are not the same file.
if [ "${R1}" == "${R2}" ]; then
	errorExit "Input error: -f and -r cannot be the same."
fi

### Checking if OUTPUT_DIR already exists - will not overwrite existing directory.ß
if [ -d "${outDir}" ]; then
	errorExit "Error: Output directory already exists. Cannot overwrite."
fi

### Get sampleID from R1 input file (Any string before 1st underscore).
sampleID=$(basename ${R1} _R1.fastq.gz)
echo "SampleID = ${sampleID}"

### Setting up direcory structure.
wrkDir=${outDir}/${sampleID} 			# All output goes here
scriptsDir=${wrkDir}/scripts 		# Scripts are stored and launched from here
mkdir -p ${wrkDir}
mkdir ${scriptsDir}
mkdir ${wrkDir}/logs 				# All slurm log files will go here

### Copy raw data to working directory.
cp ${R1} ${wrkDir} 
errorExit "Copying ${R1} to ${wrkDir} failed." ; echo "${R1} copied to working dir: ${wrkDir}"

cp ${R2} ${wrkDir}
errorExit "Copying ${R2} to ${wrkDir} failed." ; echo "${R2} copied to working dir: ${wrkDir}"

### Copy analysis scripts to working directory.
masterScripts=/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/*.slurm.sh

for file in ${masterScripts}; do
	newFileName=$(basename $file .slurm.sh).${sampleID}.slurm.sh
	cp ${file} ${scriptsDir}/${newFileName}
	errorExit "Error: copying .slurm scripts to working directory failed." ; echo "${newFileName} copied to ${scriptsDir}."
done

### Copy run.sh to working dir
runScripts='/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/run.sh'

cp ${runScripts} ${scriptsDir}/$(basename ${runScripts} .sh).${sampleID}.sh
errorExit "Error: copying run.sh to working directory failed." ; echo "run.sh copied to ${scriptsDir}"

### Copy slurmParams.txt and script_vars.txt to working dir and add universal variables
slurmParams='/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/slurmParams.txt'
cp ${slurmParams} ${scriptsDir}/slurmParams.txt
errorExit "Error: copying slurmParams.txt to working directory failed." ; echo "slurmParams.txt copied to ${scriptsDir}"
sed -i s@SAMPLEID@${sampleID}@g ${scriptsDir}/slurmParams.txt
errorExit "Error: writing SAMPLEID variable to slurmParams.txt failed." ; echo "${sampleID} written to ${scriptsDir}/slurmParams.txt"

scriptVars='/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/script_vars.txt'
cp ${scriptVars} ${scriptsDir}/script_vars.txt
errorExit "Error: copying script_vars.txt to working directory failed." ; echo "script_vars.txt copied to ${scriptsDir}"
sed -i s@SAMPLEID@${sampleID}@g ${scriptsDir}/script_vars.txt
errorExit "Error: writing SAMPLEID variable to script_vars.txt failed." ; echo "${sampleID} written to ${scriptsDir}/script_vars.txt"

### Copy tools dir to working direcory.
toolsDir='/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/tools/*'

mkdir ${scriptsDir}/tools
cp ${toolsDir} ${scriptsDir}/tools/
errorExit "Error: copying bin/tools/ to working directory failed." ; echo "bin/tools copies to working directory."

# Final comments to stdout.
echo "All output directories and analysis scripts created..."
echo "To launch analysis run:"
echo "cd ${wrkDir}/scripts && ./run.${sampleID}.sh -m all"
date
