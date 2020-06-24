#!/bin/bash
date

flowControl() {
	if [ $? -ne 0 ]
	then
		echo $1; date; exit 1;
	fi
}

module load trinity/2.8.4
export OMP_NUM_THREADS=1

find ../SAMPLEID_trinity_out/read_partitions/ -name '*inity.fasta' | \
	/apps/trinity/2.8.4/util/support_scripts/partitioned_trinity_aggregator.pl \
	--token_prefix TRINITY_DN --output_prefix ../SAMPLEID_trinity_out/Trinity.tmp
flowControl "find failed."

mv ../SAMPLEID_trinity_out/Trinity.tmp.fasta ../SAMPLEID_trinity_out.Trinity.fasta
flowControl "mv failed."

/apps/trinity/2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl \
	../SAMPLEID_trinity_out.Trinity.fasta > ../SAMPLEID_trinity_out.Trinity.fasta.gene_trans_map
flowControl "trans_map.pl failed."

mv ../SAMPLEID_trinity_out.Trinity.fasta ../SAMPLEID.trinity.fasta
	flowControl "mv failed."

sed -i 's/>TRINITY/>SAMPLEID_TRINITY/g' ../SAMPLEID.trinity.fasta
flowControl "sed failed."

# counting transcipts
count=$(grep -c "^>" ../SAMPLEID.trinity.fasta)
echo "Number of transcripts."
echo "../SAMPLEID.trinity.fasta: ${count}"

# sign off
echo "Trinity completed successfully."
date

