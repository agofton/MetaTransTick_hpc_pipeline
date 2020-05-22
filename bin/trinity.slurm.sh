#!/bin/bash

#SBATCH --job-name=Tr_SAMPLEID
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=7-00:00:00
#SBATCH --output=../logs/trinity_SAMPLEID_%A.log
#SBATCH --mem=512GB

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
	--seqType fq \
	--max_memory ${SLURM_MEM_PER_NODE}G \
	--left ${R1} \
	--right ${R2} \
	--CPU ${SLURM_NTASKS} \
	--min_contig_length 200 \
	--output ${out} \
	--no_normalize_reads

# --no_distributed_trinity_exec 	#this will only run phase 1.

# flow control
if [ $? -eq 0 ]; then
	echo "Trinity finished sucessfully: ${R1}, ${R2}" ; date
else
	echo "Trinity failed: ${R1}, ${R2}"; date; exit 1
fi

# manipulating output
mv SCRATCHDIR/SAMPLEID_trinity_out.Trinity.fasta SCRATCHDIR/SAMPLEID.trinity.fasta
sed -i 's/>TRINITY/>SAMPLEID_TRINITY/g' SCRATCHDIR/SAMPLEID.trinity.fasta

# flow control
if [ $? -ne 0 ]; then
	echo "SED failed: ${R1}, ${R2}"; date; exit 1
fi

# counting transcipts
count=$(grep -c "^>" SCRATCHDIR/SAMPLEID.trinity.fasta)

echo "Number of transcripts."
echo "SCRATCHDIR/SAMPLEID.trinity.fasta: ${count}"

# moving output back to hb-auticks
cp SCRATCHDIR/SAMPLID.trinity.fasta OUTDIR/

echo ""
date
echo ""
