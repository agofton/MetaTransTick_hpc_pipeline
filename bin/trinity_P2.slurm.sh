#!/bin/bash

flowControl() {
	if [ $? -ne 0 ]
	then
		echo $1; date; exit 1;
	else
		echo $2
	fi
}

createArray() {
	array="(`for x in ${1}
	do
		echo -n '"'
		echo -n ${x}
		echo -n '"'
	done`);"
	array=`sed -E 's@""@" \\\\\\n"@g' <<< ${array}
}

module load trinity/2.8.4
export OMP_NUM_THREADS=20

createArray "../SAMPLEID_trinity_out/recursive_trinity.cmds.chunk.??"
inArray=${array}

flowControl "Error in creating input file array for slurm array script." " "

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
	i=${SLURM_ARRAY_TASK_ID}
	/apps/trinity/2.8.4/trinity-plugins/BIN/ParaFly -c ../SAMPLEID_trinity_out/${inArray[$i]} -CPU ${SLURM_CPUS_PER_TASK} -v -shuffle
	flowControl "ParaFly failed: ${inArray[$i]}" "ParaFly Complete: ${inArray[$i]}"
else
	echo "Error: missing array index as SLURM_ARRAY_TASK_ID"
fi

