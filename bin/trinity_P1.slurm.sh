#!/bin/bash

#SBATCH --job-name=TrP1_SAMPLEID
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --output=../logs/trinity_P1_SAMPLEID_%A.log
#SBATCH --mem=512GB

# MEMORY IN GB IS HARDCODED BELOW BECAUSE SLURM ENV VAR PRINTS IN KB NOT GB

module load trinity/2.8.4
module load perl
module load bowtie/2.2.9

export OMP_NUM_THREADS=${SLURM_NTASKS}

# inputs
R1=SCRATCHDIR/SAMPLEID_R1.derep.fastq.gz
R2=SCRATCHDIR/SAMPLEID_R2.derep.fastq.gz
out=SCRATCHDIR/SAMPLEID_trinity_out

# link read_partions to /dev/shm/gof005.${SLURM_ID}/read_partitions
mkdir ${out}
rm -f ${out}/read_partitions
mkdir ${MEMDIR}/read_partitions
ln -s ${MEMDIR}/read_partitions ${out}

echo ""
date
echo ""

Trinity \
	--no_distributed_trinity_exec \
	--seqType fq \
	--max_memory 510G \
	--left ${R1} \
	--right ${R2} \
	--CPU ${SLURM_NTASKS} \
	--min_contig_length 200 \
	--output ${out} \
	--no_normalize_reads

# flow control
if [ $? -eq 0 ]; then
	echo "Trinity finished sucessfully: ${R1}, ${R2}" ; date
else
	echo "Trinity failed: ${R1}, ${R2}"; date; exit 1
fi

# break symlink to /dev/shm & move data back to working dir
mkdir ${out}/rp_tmp 	# creat new dir for data in working dir
mv ${MEMDIR}/read_partitions ${out}/rp_tmp #move data from/dev/shm to working directory
rm -f ${out}/read_partitions # remove symlink pointer
mv ${out}/rp_tmp/read_partitions ${out} # remane data to read_partitions
rm -rf ${out}/rp_tmp # remove tmp dir
cp -r ${out}/read_partitions ${out}/read_partitions_copy # make backup incase of failed jobs

# split up recursive_trinity.cmds into 10 chunks (1 chunk per node)
split -d -n l/10 \
	SCRATCHDIR/SAMPLEID_trinity_out/recursive_trinity.cmds \
	SCRATCHDIR/SAMPLEID_trinity_out/recursive_trinity.cmds.chunk.

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

echo ""
date
echo ""
