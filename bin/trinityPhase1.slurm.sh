#!/bin/bash
date
source slurmParams.txt
source script_vars.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	fi
}

module load trinity/2.11.0
module load perl
module load bowtie/2.3.4

export OMP_NUM_THREADS=${SLURM_NTASKS}
maxMemGB=$(echo ${SLURM_MEM_PER_NODE} | cut -c1-3)G

# link read_partions to /dev/shm/gof005.${SLURM_ID}/read_partitions
mkdir ${trinOutDir}
rm -f ${trinOutDir}/read_partitions
mkdir ${MEMDIR}/read_partitions
ln -s ${MEMDIR}/read_partitions ${trinOutDir}

Trinity \
	--seqType fq \
	--max_memory ${maxMemGB} \
	--left ${derepOut1} \
	--right ${derepOut2} \
	--CPU ${SLURM_NTASKS} \
	--min_contig_length 200 \
	--output ${trinOutDir} \
	--no_distributed_trinity_exec \
	--no_normalize_reads

errorExit "Trinity failed: ${derepOut1}, ${derepOut2}" 
echo "Trinity finished sucessfully: ${derepOut1}, ${derepOut2}"

# break symlink to /dev/shm & move data back to working dir
mkdir ${trinOutDir}/rp_tmp
mv ${MEMDIR}/read_partitions ${trinOutDir}/rp_tmp
rm -f ${trinOutDir}/read_partitions
mv ${trinOutDir}/rp_tmp/read_partitions ${trinOutDir}
rm -rf ${trinOutDir}/rp_tmp
cp -r ${trinOutDir}/read_partitions ${trinOutDir}/read_partitions_copy

# split up recursive_trinity.cmds into 10 chunks (1 chunk per node)
split \
	-d -n l/${nChunks} \
	${trinOutDir}/recursive_trinity.cmds \
	${trinOutDir}/recursive_trinity.cmds.chunk.

errorExit "Splitting recursive_trinity.cmds failed."

#output will be:
# 	recursive_trinity.cmds.chunk.00
# 	recursive_trinity.cmds.chunk.01
# 	recursive_trinity.cmds.chunk.02
#   recursive_trinity.cmds.chunk.03
#   recursive_trinity.cmds.chunk.04
#   recursive_trinity.cmds.chunk.05
#   recursive_trinity.cmds.chunk.06
#   recursive_trinity.cmds.chunk.07
#   recursive_trinity.cmds.chunk.08
#   recursive_trinity.cmds.chunk.09
#	etc...

date
