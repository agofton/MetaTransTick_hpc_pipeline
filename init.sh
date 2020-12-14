#!/bin/bash
date

helpMessage() {
	echo "$0 usage: -f sample_x_R1.fastq.gz <raw input R1 file> -r sample_x_R2.fastq.gz <raw input R1 file> -o <output folder> -h [show this message]";
	echo "<> = manditory argument, [] = optional argument"
	echo "";
	echo "Initialises all folders and scripts needed for trascriptome assembly and analysis. Assumes input files are in the format: sample_xxx_R1.fastq.gz & sample_xxx_R2.fastq.gz, where 'sample_xxx_' can be whatever samples ID you want. This sample ID will be carried through the whole analysis and appended onto the transcript IDs - so make sure it is not too long.";
	echo "";
	echo "Best to use complete paths!";
	echo "";
	exit 1;
}

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	fi
}

# Command line arguments
while getopts "hf:r:o:" OPT
do
	case "${OPT}"
		in
	        h) helpMessage;;
	        f) R1=${OPTARG};;
	        r) R2=${OPTARG};;
			o) outputDir=${OPTARG};;
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

# Checking if outputDir already exists - will not overwrite existing directory.
if [ -d "${outputDir}" ]; then
	errorExit "Error: Output directory already exists. Cannot overwrite."
fi

# Get sampleID from R1 input file (Any string before 1st underscore).
sampleID=$(basename ${R1} _R1.fastq.gz)
echo "SampleID = ${sampleID}"

# Setting up direcory structure.
outDir=${outputDir}/${sampleID} 	# All output goes here
scriptsDir=${outDir}/scripts 		# Scripts are stored and launched from here
mkdir -p ${outDir}
mkdir ${scriptsDir}
mkdir ${outDir}/logs 				# All slurm log files will do here

# Copy raw data to working directory.
cp ${R1} ${outDir} 
errorExit "Copying ${R1} to ${outDir} failed."
echo "${R1} copied to working dir: ${outDir}"

cp ${R2} ${outDir}
errorExit "Copying ${R2} to ${outDir} failed."
echo "${R2} copied to working dir: ${outDir}"

# Copy analysis scripts to working directory.
for x in /datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/*.slurm.sh
do
	newFilename=$(basename $x .slurm.sh).${sampleID}.slurm.sh
	cp ${x} ${scriptsDir}/${newFilename}
	errorExit "Error: copying .slurm scripts to working directory failed."
	echo "${newFilename} copied to ${scriptsDir}."
done

# Copy run.sh to working dir
runScript=/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/run.sh
cp ${runScript} ${scriptsDir}/$(basename ${runScript} .sh).${sampleID}.sh
errorExit "Error: copying run.sh to working directory failed."
echo "run.sh copied to ${scriptsDir}"

# Copy slurmParams.txt to working dir and add universal variables
slurmParams=/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/slurmParams.txt
cp ${slurmParams} ${scriptsDir}/slurmParams.txt
errorExit "Error: copying slurmParams.txt to working directory failed."
echo "slurmParams.txt copied to ${scriptsDir}."

sed -i s@SAMPLEID@${sampleID}@g ${scriptsDir}/slurmParams.txt
errorExit "Error: writing SAMPLEID variable to slurmParams.txt failed."
echo "${sampleID} written to ${scriptsDir}/slurmParams.txt"

# Copy tools dir to working direcory.
mkdir ${scriptsDir}/tools
cp /datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/tools/* ${scriptsDir}/tools/
errorExit "Error: copying bin/tools/ to working directory failed."
echo "bin/tools copies to working directory."

# Final comments to stdout.
echo "All output directories and analysis scripts created..."
echo ""
echo "To launch analysis:"
echo "cd ${outDir}/scripts && ./run.${sampleID}.sh -m all"
echo ""
date
