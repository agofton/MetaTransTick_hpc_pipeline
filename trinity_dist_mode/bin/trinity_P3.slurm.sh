#!/bin/bash

#SBATCH -J TrP3_SAMPLEID
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=../logs/trinity_P3_SAMPLEID_%A.out
#SBATCH --mem=8GB

module load trinity/2.8.4

export OMP_NUM_THREADS=1

echo ""
date
echo ""

# final job putting together all assemblies into xxx.Trinity.fasta

find \
	OUTDIR/SAMPLEID_trinity_out/read_partitions/ \
	-name '*inity.fasta' | \
	/apps/trinity/2.8.4/util/support_scripts/partitioned_trinity_aggregator.pl \
	--token_prefix TRINITY_DN \
	--output_prefix OUTDIR/SAMPLEID_trinity_out/Trinity.tmp

# flow control
if [ $? -ne 0 ]; then
	echo ""; echo "find failed."; date; exit 1
else
	echo ""; echo "FIND complete."; date; echo ""
fi

mv OUTDIR/SAMPLEID_trinity_out/Trinity.tmp.fasta \
	OUTDIR/SAMPLEID_trinity_out.Trinity.fasta

# flow control
if [ $? -ne 0 ]; then
	echo ""; echo "mv failed."; date; exit 1
else
	echo ""; echo "mv complete."; date; echo ""
fi

/apps/trinity/2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl \
	OUTDIR/SAMPLEID_trinity_out.Trinity.fasta > \
	OUTDIR/SAMPLEID_trinity_out.Trinity.fasta.gene_trans_map

# flow control
if [ $? -ne 0 ]; then
	echo ""; echo "trans_map.pl failed."; date; exit 1
else
	echo ""; echo "trans_map.pl complete."; date; echo ""
fi

# manipulating output
mv OUTDIR/SAMPLEID_trinity_out.Trinity.fasta OUTDIR/SAMPLEID.trinity.fasta

# flow control
if [ $? -ne 0 ]; then
	echo ""; echo "mv failed."; date; exit 1
else
	echo ""; echo "mv complete."; date; echo ""
fi

sed -i 's/>TRINITY/>SAMPLEID_TRINITY/g' OUTDIR/SAMPLEID.trinity.fasta

# flow control
if [ $? -ne 0 ]; then
	echo ""; echo "SED failed."; date; exit 1
else
	echo ""; echo "SED complete."; date; echo ""
fi

# counting transcipts
count=$(grep -c "^>" OUTDIR/SAMPLEID.trinity.fasta)

echo "Number of transcripts."
echo "OUTDIR/SAMPLEID.trinity.fasta: ${count}"

# sign off
echo ""
echo "Trinity completed successfully."
date
echo ""

