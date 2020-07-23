#!/bin/bash
date
source ../slurmParams.txt

helpMess() {

}

errorExit() {
	if [ $? -ne 0 ]; then
		echo $1; date, exit 1
	fi
}

outFile=../../summary_report.txt

# Headers are:
# sampleID #PEreads #PEreadsCutadapt %PEreadsCutadapt #PEreadsTrimmomatic %PEreadsTrimmomatic #UPreadsTrimmomatic %UPreadsTrimmomatic #PEreadsDerep %PEreadsDerep #Transcripts
#
#
