#!/bin/bash

#SBATCH -J TrP2_SAMPLEID
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=06:00:00
#SBATCH --output=../logs/trinity_P2_SAMPLEID_%A_%a.out
#SBATCH --mem=128GB
#SBATCH --array=0-9

module load trinity/2.8.4

export OMP_NUM_THREADS=20

# array will need to be manually expanded if more chunks are used. def(10)
in_array=("recursive_trinity.cmds.chunk.00" \
	"recursive_trinity.cmds.chunk.01" \
	"recursive_trinity.cmds.chunk.02" \
	"recursive_trinity.cmds.chunk.03" \
	"recursive_trinity.cmds.chunk.04" \
	"recursive_trinity.cmds.chunk.05" \
	"recursive_trinity.cmds.chunk.06" \
	"recursive_trinity.cmds.chunk.07" \
	"recursive_trinity.cmds.chunk.08" \
	"recursive_trinity.cmds.chunk.09");

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
	i=${SLURM_ARRAY_TASK_ID}

	/apps/trinity/2.8.4/trinity-plugins/BIN/ParaFly \
		-c OUTDIR/SAMPLEID_trinity_out/${in_array[$i]} \
		-CPU 20 \
		-v \
		-shuffle

	if [ $? -ne 0 ]
	then
		echo ""; echo ParaFly failed: ${in_array[$i]}; date; exit 1
	else
		echo ""; echo date; echo ""
	fi

else
	echo "Error: missing array index as SLURM_ARRAY_TASK_ID"
fi

