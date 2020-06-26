#!/bin/bash

flowControl() {
	if [ $? -ne 0 ]
	then
		echo $1; date; exit 1;
	else
		echo $2
	fi
}

module load trinity/2.8.4
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Create array of input files
arrayFile=trinP2InputArray.txt
ls -1 ../SAMPLEID_trinity_out/recursive_trinity.cmds.chunk.?? > ${arrayFile}
inArray=( `cut -f 1 ${arrayFile}` );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
	i=${SLURM_ARRAY_TASK_ID}
	/apps/trinity/2.8.4/trinity-plugins/BIN/ParaFly -c ${inArray[$i]} -CPU ${SLURM_CPUS_PER_TASK} -v -shuffle
	flowControl "ParaFly failed: ${inArray[$i]}" "ParaFly Complete: ${inArray[$i]}"
else
	echo "Error: missing array index as SLURM_ARRAY_TASK_ID"
fi
