#!/bin/bash
date

flowControl() {
	if [ $? -ne 0 ]
	then
		echo $1; date; exit 1;
	fi
}

module load trinity/2.8.4
export OMP_NUM_THREADS=1 	# Shouldn't ever need more than this

##################################################################################################
# These steps (except the last sed command) are taken directly from the normal Trinity workflow. #
##################################################################################################

# Find and agregate all final assemblies
find ../SAMPLEID_trinity_out/read_partitions/ -name '*inity.fasta' | \
	/apps/trinity/2.8.4/util/support_scripts/partitioned_trinity_aggregator.pl \
	--token_prefix TRINITY_DN --output_prefix ../SAMPLEID_trinity_out/Trinity.tmp

flowControl "find failed."

# Move file
mv ../SAMPLEID_trinity_out/Trinity.tmp.fasta ../SAMPLEID_trinity_out.Trinity.fasta

flowControl "mv failed."

# Create gene-transcript map
/apps/trinity/2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl \
	../SAMPLEID_trinity_out.Trinity.fasta > ../SAMPLEID_trinity_out.Trinity.fasta.gene_trans_map

flowControl "trans_map.pl failed."

# Move file to final destiation
mv ../SAMPLEID_trinity_out.Trinity.fasta ../SAMPLEID.trinity.fasta

flowControl "mv failed."

# Add sampleID preflix to each transcript
sed -i 's/>TRINITY/>SAMPLEID_TRINITY/g' ../SAMPLEID.trinity.fasta

flowControl "sed failed."

# Counting transcipts
count=$(grep -c "^>" ../SAMPLEID.trinity.fasta)
echo "Number of transcripts."
echo "../SAMPLEID.trinity.fasta: ${count}"

# sign off
echo "Trinity completed successfully."
date
