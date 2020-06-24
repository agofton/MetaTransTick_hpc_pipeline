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

flowControl() {
	if [ $? -ne 0 ]
	then
		echo $1; date; exit 1;
	else
		echo $2;
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

# check input files
if [ ! -f "${R1}" ]
then
	echo "Input R1 file does not exist."; date; exit 1
fi
##
if [ ! -f "${R2}" ]
then
	echo "Input R2 file does not exits."; date; exit 1
fi
##
if [ "${R1}" == "${R2}" ]
then
	echo "Input error: -f and -r cannot be the same."; date; exit 1
fi
#
if [ -d "${outputDir}" ]
then
	echo "Output directory already exists. Cannot overwrite."; date; exit 1
fi

# Get sample ID from R1 input. Assumes file names are in format whatever_sampleX_R1.fastq.gz & whatever_sampleX_R2.fastq.gz
sampleID=$(basename ${R1} _R1.fastq.gz)
echo "SampleID = ${sampleID}"

# setting up direcory structure
outDir=${outputDir}/${sampleID} 	# All output goes here
scriptsDir=${outDir}/scripts 		# Scripts are stored and launched from here
mkdir -p ${outDir}
mkdir ${scriptsDir}
mkdir ${outDir}/logs 				# All slurm log files will do here

# copy raw data to sample working directory
cp ${R1} ${outDir}/
flowControl "Copying raw R1 failed." "${R1} copied to ${outDir}."
cp ${R2} ${outDir}/
flowControl "Copying raw R2 failed." "${R2} copied to ${outDir}." 

# copy analysis scripts to sample scripts directory & editing scripts with universal variables
for x in /datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/*.slurm.sh
do
	newFilename=$(basename $x .slurm.sh).${sampleID}.slurm.sh
	cp ${x} ${scriptsDir}/${newFilename}
	flowControl "Copying scripts failed." "Copying scripts..."
	sed -i s@SAMPLEID@${sampleID}@g ${scriptsDir}/${newFilename} 	# SAMPLEID is the sample ID placeholder in the master scripts.
	flowControl "Writing SAMPLEID variable failed."
	sed -i s@OUTDIR@${outDir}@g ${scriptsDir}/${newFilename} 	# OUTDIR is the output directory name placeholder in the master scripts.
	flowControl "Writing OUTDIR variable failed."
done

# copy run script
runScript=/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/run.sh
newRunScript=$(basename $runScript .sh).${sampleID}.sh
cp ${runScript} ${scriptsDir}/${newRunScript}
flowControl "Copying run.sh failed." "Copying run.sh complete."
sed -i s/SAMPLEID/${sampleID}/g ${scriptsDir}/${newRunScript} 	# SAMPLEID is the sample placeholder in the master scripts
flowControl "Writing SAMPLEID variable to run.sh failed."

# copy slurmParams.txt
slurmParams=/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/slurmParams.txt
cp ${slurmParams} ${scriptsDir}/slurmParams.txt
flowControl "Copying slurmParams.txt failed." "Copying slurmParams.txt complete."
sed -i s@SAMPLEID@${sampleID}@g ${scriptsDir}/slurmParams.txt
flowControl "Writing SAMPLEID variable to slurmParams.txt failed."

echo "Output and tmp directories and analysis scripts created..."
echo ""
echo "Run commands:"
echo "cd ${outDir}/scripts"
echo "nohup run.${sampleID}.sh &"
echo ""
date

