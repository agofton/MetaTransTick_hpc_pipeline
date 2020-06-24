#!/bin/bash
date

flowControl() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	else
		echo $2
	fi
}

source slurmParams.txt
flowControl "Cannot source variables from ./slurmParams.txt. File does not exits." " "

module load trinity/2.8.4
module load perl
module load bowtie/2.2.9

export OMP_NUM_THREADS=${SLURM_NTASKS}

# inputs
R1=../SAMPLEID_R1.derep.fastq.gz
R2=../SAMPLEID_R2.derep.fastq.gz
out=../SAMPLEID_trinity_out
maxMemGB=$(echo ${SLURM_MEM_PER_NODE} | cut -c1-3)G

# link read_partions to /dev/shm/gof005.${SLURM_ID}/read_partitions
mkdir ${out}
rm -f ${out}/read_partitions
mkdir ${MEMDIR}/read_partitions
ln -s ${MEMDIR}/read_partitions ${out}

Trinity --seqType fq --max_memory ${maxMemGB} --left ${R1} --right ${R2} --CPU ${SLURM_NTASKS} --min_contig_length 200 --output ${out} --no_normalize_reads --no_distributed_trinity_exec 
flowControl "Trinity failed: ${R1}, ${R2}" "Trinity finished sucessfully: ${R1}, ${R2}" 

# break symlink to /dev/shm & move data back to working dir
mkdir ${out}/rp_tmp 										# creat new dir for data in working dir
mv ${MEMDIR}/read_partitions ${out}/rp_tmp 					# move data from/dev/shm to working directory
rm -f ${out}/read_partitions 								# remove symlink pointer
mv ${out}/rp_tmp/read_partitions ${out} 					# remane data to read_partitions
rm -rf ${out}/rp_tmp 										# remove tmp dir
cp -r ${out}/read_partitions ${out}/read_partitions_copy 	# make backup incase of failed jobs

# split up recursive_trinity.cmds into 10 chunks (1 chunk per node)
split -d -n l/${nChunks} \
	../SAMPLEID_trinity_out/recursive_trinity.cmds \
	../SAMPLEID_trinity_out/recursive_trinity.cmds.chunk.
flowControl "Splitting recursive_trinity.cmds failed." " "
#output will be:
# 	recursive_trinity.cmds.chunk.00
# 	recursive_trinity.cmds.chunk.01
# 	recursive_trinity.cmds.chunk.02
#   recursive_trinity.cmds.chunk.03
#   recursive_trinity.cmds.chunk.04
#   recursive_trinity.cmds.chunk.05
#	recursive_trinity.cmds.chunk.06
# 	recursive_trinity.cmds.chunk.07
# 	recursive_trinity.cmds.chunk.08
#   recursive_trinity.cmds.chunk.09
# 	etc...

date
