#!/bin/bash
date
source slurmParams.txt

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date; exit 1;
	fi
}

# Load software
module load bowtie/2.3.4
module load samtools/1.10.0
module load python/3.7.2

# Make output directory
mkdir -p ${conCovOutDir}

# Build bowtie2 index of whole transcriptome assembly
bowtie2-build --threads ${SLURM_CPUS_PER_TASK} --quiet ${trinFasta} ${bt2index}
errorExit "Bowtie2 indexing of trinity assembly failed ..."
echo "BWA indexing complete. "
echo -n $(date)

# Map all QCd reads to whole transcriptome with bowtie2
bowtie2 -p ${SLURM_CPUS_PER_TASK} -x ${bt2index} -1 ${derepOut1} -2 ${derepOut2} -S ${samOut} --fast
errorExit "Bowtie2 failed :("
echo "Bowtie2 complete."

# Generate mapping statistics
samtools sort -@ ${SLURM_CPUS_PER_TASK} ${samOut} > ${samSorted}
errorExit "Samtools sort failed."

samtools coverage ${samSorted} > ${covSum}
errorExit "Samtools coverage failed!"
echo "Samtools coverage complete. "
echo -n $(date)

samtools flagstat -@ ${SLURM_CPUS_PER_TASK} ${samSorted} > ${samFlagstat}

# Fuction to do tpm calculations in python
TPM_calc() {
	python3 - <<END
import pandas as pd
df = pd.read_csv("${covSum}", sep = "\t")
df["RPK"] = df["numreads"] / (df["endpos"] / 1000)
total_RPK = df["RPK"].sum()
df["TPM"] = df["RPK"] / total_RPK
df.to_csv("${covSum}", sep = "\t", index = False, header = True)
END
}
TPM_calc
errorExit "TPM calc. in python failed."
echo "TPM calc. complete."

# Functiong to do % mapping statistics 
mapping_rate() {
	python3 - <<END
import pandas as pd
import os

# Get numbers from flagstat.txt by converting to df
df = pd.read_csv("${samFlagstat}", sep = '\s+', header = None, names=["col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11"])
mapped = (df.iloc[0, 0] - df.iloc[1, 0] - df.iloc[2, 0]) - (df.iloc[0, 0] - df.iloc[4, 0])

# Count number of reads in trinity.fasta
trin_fasta = open(${trinFasta})
reads = 0
for line in trin_fasta:
	if line.startswith(">"):
		reads = reads + 1
trin_fasta.close()

# Do maths
pct_mapped = (mapped / reads) * 100
rounded_pct_mapped = round(pct_mapped, 5)

# Append results to bottom of df
list = [{"col1": "Total trinity contigs:", "col2": reads}. {"col1": "Number of mapped reads:", "col2": mapped} {"col1": "Pct reads mapped:", "col2": rounded_pct_mapped}]
list_df = pd.DataFrame(list)
df2 = df.append(list_df)

# Write df2 flagstat.txt
df2.to_csv("${samFlagstat}", index = False, header = False, sep = ' ')
END
}
mapping_rate
errorExit "Mapping pct calcs failed."
echo "Mapping pct calcs complete."



date
