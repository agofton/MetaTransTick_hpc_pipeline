#!/bin/bash

helpMessage() {
	echo "$0 usage: -f sample_x_R1.fastq.gz <raw input R1 file> -r sample_x_R2.fastq.gz <raw input R1 file> -o <output folder> -h show this message";
	echo "";
	echo "Initialises all folders and scripts needed for analysis. Assumes input files are in the format: sample_xxx_R1.fastq.gz & sample_xxx_R2.fastq.gz, where 'sample_xxx_' can be whatever samples ID you want. This sample ID will be carried through the whole analysis and appended onto the transcript IDs - so make sure it is not too long.";
	echo "";
	echo "USE COMPLETE PATH TO AVOID FUCK UPS!";
	echo "";
	exit 1;
}

errorExit() {
	msg=$1
	exitCode=${2:-1}
	echo ${msg};
	exit ${exitCode};
}

# Command line arguments
while getopts "hf:r:o:" OPT
do
	case "${OPT}"
		in
	        h) helpMessage;;
	        f) R1=${OPTARG};;
	        r) R2=${OPTARG};;
			o) output_dir=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

# check input files
if [ ! -f "${R1}" ]; then
	errorExit "Input R1 file does not exist."
fi

if [ ! -f "${R2}" ]; then
	errorExit "Input R2 file does not exits."
fi

if [ "${R1}" == "${R2}" ]; then
	errorExit "Input error: -f and -r cannot be the same."
fi

if [ -d "${output_dir}" ]
	errorExit "Output directory already exists. Cannot overwrite."
fi

# Get sample ID from R1 and R2 input.
# Assumes file names are in format whatever_sampleX_R1.fastq.gz & whatever_sampleX_R2.fastq.gz
sampleID=$(basename ${R1} _R1.fastq.gz)

# setting up direcory structure and scripts
out_dir=${output_dir}/${sampleID}
scripts_dir=${out_dir}/scripts

# make directories
mkdir -p ${out_dir}
mkdir ${scripts_dir}
mkdir ${out_dir}/logs

echo ""

# copy raw data to sample working directory
cp ${R1} ${out_dir}/
	
if [ $? -eq 0 ]; then
	echo "${R1} copied to ${out_dir}."
else
	errorExit "Copying raw R1 failed."
fi

cp ${R2} ${out_dir}/

if [ $? -eq 0 ]; then
	echo "${R2} copied to ${out_dir}."
else
	errorExit "Copying raw R2 failed."
fi

echo ""

# copy analysis scripts to sample scripts directory & editing scripts with new sampleID and filename

for x in /datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/*.slurm.sh
do
	new_filename=$(basename $x .slurm.sh).${sampleID}.slurm.sh
	cp ${x} ${scripts_dir}/${new_filename}
	sed -i s@SAMPLEID@${sampleID}@g ${scripts_dir}/${new_filename} 	# SAMPLEID is the sample ID placeholder in the master scripts.
	sed -i s@OUTDIR@${out_dir}@g ${scripts_dir}/${new_filename} 	# OUTDIR is the output directory name placeholder in the master scripts.

	if [ $? -eq 0 ]; then
		echo "${new_filename} copied to working directory."
	else
		errorExit "Copying scripts failed."
	fi
done

# copy run script
run_script=/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/bin/run.sh
new_run_script=$(basename $run_script .sh).${sampleID}.sh

cp ${run_script} ${scripts_dir}/${new_run_script}
sed -i s/SAMPLEID/${sampleID}/g ${scripts_dir}/${new_run_script} 	# SAMPLEID is the sample placeholder in the master scripts
	
	if [ $? -eq 0 ]; then
		echo "${new_run_script} copied to working directory."
	else
		errorExit "Copying run.sh failed."
	fi

echo ""
echo "Output and tmp directories and analysis scripts created..."
echo ""
echo "Run ${out_dir}/scripts/run.${sampleID}.sh to start analysis."
echo ""
echo "Run commands:"
echo "cd ${out_dir}/scripts"
echo "nohup run.${sampleID}.sh &"
echo ""
date

