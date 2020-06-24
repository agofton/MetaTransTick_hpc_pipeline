#!/bin/bash
date

helpMessage() {
	echo "$0 usage: -f sample_x_R1.fastq.gz <raw input R1 file> -r sample_x_R2.fastq.gz <raw input R1 file> -o <output folder> -h show this message";
	echo "";
	echo "Initialises all folders and scripts needed for analysis. Assumes input files are in the format: sample_xxx_R1.fastq.gz & sample_xxx_R2.fastq.gz, where 'sample_xxx_' can be whatever samples ID you want. This sample ID will be carried through the whole analysis and appended onto the transcript IDs - so make sure it is not too long.";
	echo "";
	echo "USE COMPLETE PATH TO AVOID FUCK UPS!";
	echo "";
	exit 1;
}

flowControl() {							# $1 is error message
	if [ $? -ne 0 ]
	then
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

# Check R1 exists
if [ ! -f "${R1}" ]; then
	echo "Input R1 file does not exist."; date; exit 1
fi
# Check R2 exists
if [ ! -f "${R2}" ]; then
	echo "Input R2 file does not exits."; date; exit 1
fi
# Check R1 & R2 are not the same file
if [ "${R1}" == "${R2}" ]; then
	echo "Input error: -f and -r cannot be the same."; date; exit 1
fi
# Check if outputDir already exists - cancel overwrite
if [ -d "${outputDir}" ]; then
	echo "Output directory already exists. Cannot overwrite."; date; exit 1
fi

# Get sample ID from R1 input. Assumes file names are in format whatever_sampleX_R1.fastq.gz & whatever_sampleX_R2.fastq.gz
sampleID=$(basename ${R1} _R1.fastq.gz)
echo "SampleID = ${sampleID}"

# setting up direcory structure
outDir=${outputDir}/${sampleID} 	# All output goes here
scriptsDir=${outDir}/scripts 			# Scripts are stored and launched from here
mkdir -p ${outDir}
mkdir ${scriptsDir}
mkdir ${outDir}/logs 							# All slurm log files will do here

# copy raw data to sample working directory
cp ${R1} ${outDir}/ && flowControl "Copying ${R1} to ${outDir} failed."
cp ${R2} ${outDir}/ && flowControl "Copying ${R2} to ${outDir} failed."

# copy analysis scripts to working directory & editing scripts with universal variables
for x in /datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/*.slurm.sh
do
	newFilename=$(basename $x .slurm.sh).${sampleID}.slurm.sh
	cp ${x} ${scriptsDir}/${newFilename} && flowControl "Copying scripts failed." "Copying scripts..."
	sed -i s@SAMPLEID@${sampleID}@g ${scriptsDir}/${newFilename} && flowControl "Writing SAMPLEID variable failed."
	sed -i s@OUTDIR@${outDir}@g ${scriptsDir}/${newFilename} && flowControl "Writing OUTDIR variable failed."
done

# copy run.sh to working dir and add universal variables
runScript=/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/run.sh
newRunScript=$(basename $runScript .sh).${sampleID}.sh
cp ${runScript} ${scriptsDir}/${newRunScript} && flowControl "Copying run.sh failed." "Copying run.sh complete."
sed -i s/SAMPLEID/${sampleID}/g ${scriptsDir}/${newRunScript} && flowControl "Writing SAMPLEID variable to run.sh failed."

# copy slurmParams.txt to working dir and add universal variables
slurmParams=/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/slurmParams.txt
cp ${slurmParams} ${scriptsDir}/slurmParams.txt && flowControl "Copying slurmParams.txt failed." "Copying slurmParams.txt complete."
sed -i s@SAMPLEID@${sampleID}@g ${scriptsDir}/slurmParams.txt && flowControl "Writing SAMPLEID variable to slurmParams.txt failed."

# final comments to stdout
echo "Output and tmp directories and analysis scripts created..."
echo ""
echo "Run commands:"
echo "cd ${outDir}/scripts"
echo "nohup run.${sampleID}.sh &"
echo ""
date
