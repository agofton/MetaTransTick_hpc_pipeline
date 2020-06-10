#!/bin/bash

hm1="Initialised all folders and scripts needed for analysis. Assumes input files are in the format: sample_xxx_R1.fastq.gz & sample_xxx_R2.fastq.gz, where 'sample_xxx_' can be whatever samples ID you want. This sample ID will be carried through the whole analysis and appended onto the transcript IDs - so make sure it is not too long."
hm2="usage: ./init.sh -f fwd_file_R1.fastq.gz -r rev_file_R2.fastq.gz -o output_directory -t tmp_dir"
hm3="A folder will be created in [-o] from the sampleID from the input file names, all output will go there. [-t] is located in /scratch1/gof005/tmp_dir all compute comes from here."
hm4="USE COMPLETE PATH TO AVOID FUCK UPS!!!!!"
# Command line arguments
while getopts hf:r:o:t: option
do
	case "${option}"
		in
	        h) echo ${hm1}
			   echo ${hm2}
			   echo ${hm3}
			   echo ${hm4}
	           echo ""
	           exit;;
	        f) R1=${OPTARG};;
	        r) R2=${OPTARG};;
			o) output_dir=${OPTARG};;
			t) tmp=${OPTARG};;
	        :) printf "missing argument for  -%s\n" "$OPTARG" >&2
	           echo "$usage" >&2
	           exit 1;;
	   	   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
	           echo "$usage" >&2
	           exit 1;;
	esac
done
shift $((OPTIND - 1))

echo ""

# Get sample ID from R1 and R2 input.
# Assumes file names are in format whatever_sampleX_R1.fastq.gz & whatever_sampleX_R2.fastq.gz
sampleID=$(basename ${R1} _R1.fastq.gz)

# setting up direcory structure and scripts
out_dir=${output_dir}/${sampleID}
scripts_dir=${out_dir}/scripts
sf=/scratch1/gof005/${tmp}

# make directories
mkdir -p ${out_dir}
mkdir ${scripts_dir}
mkdir ${out_dir}/logs
mkdir ${sf}

echo ""
# move data to scratch1/
cp ${R1} ${sf}/

if [ $? -eq 0 ]; then
	echo "${R1} copied to ${sf}."
else
	exit
fi

cp ${R2} ${sf}/

if [ $? -eq 0 ]; then
	echo "${R2} copied to ${sf}."
else
	exit
fi

echo ""
# move input data to sample working directory
cp ${R1} ${out_dir}/
	
if [ $? -eq 0 ]; then
	echo "${R1} copied to ${out_dir}."
else
	exit
fi

cp ${R2} ${out_dir}/

if [ $? -eq 0 ]; then
	echo "${R2} copied to ${out_dir}."
else
	exit
fi

echo ""
# copy analysis scripts to sample scripts directory & editing scripts with new sampleID and filename
for x in /datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/trinity_dist_mode/bin/*.slurm.sh
do
	new_filename=$(basename $x .slurm.sh).${sampleID}.slurm.sh
	cp ${x} ${scripts_dir}/${new_filename}
	sed -i s@SAMPLEID@${sampleID}@g ${scripts_dir}/${new_filename} # SAMPLEID is the sample ID placeholder in the master scripts.
	sed -i s@SCRATCHDIR@${sf}@g ${scripts_dir}/${new_filename} # SCRATCHDIR is the /scratch1/gof005/ directory name placeholder in the master scripts.
	sed -i s@OUTDIR@${out_dir}@g ${scripts_dir}/${new_filename} # OUTDIR is the output directory name placeholder in the master scripts.

	if [ $? -eq 0 ]; then
		echo "${new_filename} copied to working directory."
	else
		exit
	fi
done

# copy run script
run_script=/datasets/work/hb-austicks/work/Project_Phoenix/data/MetaTransTick_hpc_pipeline/trinity_dist_mode/bin/run.sh
new_run_script=$(basename $run_script .sh).${sampleID}.sh

cp ${run_script} ${scripts_dir}/${new_run_script}
sed -i s/SAMPLEID/${sampleID}/g ${scripts_dir}/${new_run_script} # SAMPLEID is the sample placeholder in the master scripts
	
if [ $? -eq 0 ]; then
	echo "${new_run_script} copied to working directory."
else
	exit
fi

echo ""
echo "Output and tmp directories and analysis scripts created..."
echo ""
echo "Input files copied to ${sf} for analysis."
echo ""
echo "Run ${out_dir}/scripts/run.${sampleID}.sh to start analysis."
echo ""
echo "Run commands:"
echo "cd ${out_dir}/scripts"
echo "nohup run.${sampleID}.sh &"
echo ""
date

