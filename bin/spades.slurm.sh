#!/bin/bash

#SBATCH --job-name=Sp_SAMPLEID
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
#SBATCH --output=../logs/spades.SAMPLEID_%A.log
#SBATCH --mem=512GB

module load spades/3.12.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

in1=SCRATCHDIR/SAMPLEID_R1.derep.fastq.gz
in2=SCRATCHDIR/SAMPLEID_R2.derep.fastq.gz
out=SCRATCHDIR/SAMPLEID_spades_out

echo ""
date
echo ""

spades.py \
	--rna \
	--threads ${SLURM_CPUS_PER_TASK} \
	--memory ${SLURM_MEM_PER_NODE} \
	-1 ${in1} \
	-2 ${in2} \
	-o ${out}

# flow control
if [ $? -eq 0 ]; then
	echo "SPAdes-rna finished sucessfully: ${in1}, ${in2}"; date
else
	echo "SPAdes-rna failed: ${in1}, ${in2}"; date; exit 1
fi

# copying  transcripts.fasta to sample working dir
mv ${out}/transcripts.fasta SCRATCHDIR/SAMPLEID.spades.fasta

# adding sample ID to start of contig ID
sed -i 's/>NODE/>SAMPLEID_NODE/g' SCRATCHDIR/SAMPLID.spades.fasta

# flow control
if [ $? -ne 0 ]; then
	echo "SED failed: ${R1}, ${R2}"; date; exit 1
fi

# counting transcipts
count=$(grep -c "^>" SCRATCHDIR/SAMPLEID.spades.fasta)

echo "Number of transcripts."
echo "SCRATCHDIR/SAMPLEID.spades.fasta: ${count}"

# copying output back to hb-austicks
cp SCRATCHDIR/SAMPLEID.spades.fasta OUTDIR/

# print script end time
echo ""
date
echo ""
