#!/bin/bash
date

module load bowtie/2.3.4

errorExit() {
	if [[ $? -ne 0 ]]; then
		echo $1; date; exit1
	fi
}

outDir=../transcript_coverage
indexDir=${outDir}/bt2Index
db=${indexDir}/SAMPLEID.trinity

trimPEr1=../SAMPLEID_R1.QC.paired.fastq.gz
trimPEr2=../SAMPLEID_R2.QC.paired.fastq.gz

unConc=
alConc=

mkdir -p ${indexDir}

bowtie2-build --threads 20 ../SAMPLEID.trinity.fasta ${db}

# allow no missmatches
bowtie2 -x ${db} \
  	-1 ${trimPEr1} \
	-2 ${trimPEr2} \
	-f -t --fast --un-unal --threads 20 \
	--un-conc-gz --al-conc-gz
