#!/bin/bash
date

module load bowtie/2.3.4
module load vsearch/2.13.4

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit 1
	else
		echo $2
	fi
}

# Variables
topHits=../SAMPLEID.trinity.blast.topHits.txt
tickSeqIds=../SAMPLEID.trinity.tickSeqIds.txt
tickSeqsFasta=../SAMPLEID.trinity.tickSeqs.fasta
usearch9=/home/gof005/usearch9.2
trinityFasta=../SAMPLEID.trinity.fasta
bt2Threads=20
bt2OutDir=../SAMPLEID_mapToHost
bt2Index=${bt2OutDir}/$(basename ${tickSeqsFasta} .fasta)
bt2SamOut=${bt2OutDir}/SAMPLEID.mappedToHost.sam
bt2UnconcGz=${bt2OutDir}/SAMPLEID.mappedToHost.unaligned.R%.fasta.gz
bt2AlconcGz=${bt2OutDir}/SAMPLEID.mappedToHost.aligned.R%.fasta.gz

# Find host reads from blast output and extract from .fasta
grep 'Ixodes\|Amblyomma\|Haemaphysalis\|Bothriocroton\|Rhipicephalus\|Boophilus\|Dermacentor\|Hyalomma\|Ornithodoros\|Argas' ${topHits} | \
	awk -F '\t' '{print $1}' | sort | sort -u > ${tickSeqIds}

vsearch --fastx_getseqs ${trinityFasta} --fastaout ${tickSeqsFasta} --labels ${tickSeqIds}

# Build bt2 index and run bt2
mkdir ${bt2OutDir}

bowtie2-build --threads ${bt2Threads} ${tickSeqsFasta} ${bt2Index}

bowtie2 -x ${bt2Index} -1 ${trimPEr1} -2 ${trimPEr2} \
		-f -t --fast --un-unal --threads ${bt2Threads} \
		--un-conc-gz ${bt2UnconcGz} --al-conc-gz ${bt2AlconcGz} -S ${bt2SamOut}





