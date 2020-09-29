#!/bin/bash
date

source slurmParams.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	fi
}

module load trinity/2.8.4

export OMP_NUM_THREADS=1 	# Shouldn't ever need more than this

######################################################################################## These steps (except the last sed command) are taken directly from the normal Trinity workflow. #
#######################################################################################

# Find and agregate all final assemblies
find ${trinOutDir}/read_partitions/ -name '*inity.fasta' | \
	/apps/trinity/2.8.4/util/support_scripts/partitioned_trinity_aggregator.pl \
	--token_prefix TRINITY_DN --output_prefix ${trinOutDir}/Trinity.tmp

errorExit "find failed."

# Move file
mv ${trinOutDir}/Trinity.tmp.fasta ${trinOutDir}.Trinity.fasta

errorExit "mv failed."

# Create gene-transcript map
/apps/trinity/2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl \
	${trinOutDir}.Trinity.fasta > ${trinOutDir}.Trinity.fasta.gene_trans_map

errorExit "trans_map.pl failed."

# Move file to final destiation
mv ${trinOutDir}.Trinity.fasta ${trinFasta}

errorExit "mv failed."

# Add sampleID preflix to each transcript
sed -i s/>TRINITY/>${id}_TRINITY/g ${trinFasta}

errorExit "sed failed."

# Counting transcipts
count=$(grep -c "^>" ${trinFasta})
echo "Number of transcripts."
echo "${trinFasta}: ${count}"

# sign off
echo "Trinity completed successfully."
date
