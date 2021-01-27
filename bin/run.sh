#!/bin/bash

source slurmParams.txt

# Text colors
BLACK='\033[0;01m'
RED='\033[0;31m'
GREEN='\033[0;32m'
ORANGE='\033[0;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
GRAY='\033[0;37m'
NC='\033[0m' # No Color (White)

helpMessage() {
	printf "${ORANGE}./run.sh usage: 
	-m <all, derep, trinityP1, trinityP2, trinityP3, blast+diamond, onlyQC, onlyDerep, onlyTrinity, onlyTrinityP1, onlyTrinityP2, onlyTrinityP3, onlyBlast, onlyDiamond, onlyConCov> 
	-h [show this message]\n"
	printf "${ORANGE}Running with -m flags "all", "derep",  "trinityP1", "trinityP2", "trinityP3", or "blast+diamond" will start the pipeline at that step and then run to 	the end. \n"
	printf "${ORANGE}Running with -m flags "onlyQC", "onlyDerep",  "onlyTrinityP1", "onlyTrinityP2", "onlyTrinityP3", "onlyBlast", "onlyDiamond", or "onlyConCov" will run only that module.\n"
}

# Command line arguments
while getopts "hm:" OPTION
do
	case "${OPTION}"
		in
	        h) helpMessage;;
	        m) module=${OPTARG};;
			/?) printf "${RED}Invalid option: -$OPTARG" 1>&2
	esac
done
shift $((OPTIND - 1))

# define job options
all() {
	printf "\n"
	printf "${BLUE}Running work flow steps QC, Dereplication, Trinity phase 1, Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage.\n"
	printf "\n"
	# QC - no dependency
	job1=$(sbatch -J ${qcJobName} -N ${qcNodes} -n ${qcTasks} -c ${qcCPUsPerTask} -t ${qcTime} --mem ${qcMem} -o ${qcLog} ${qc})
	job1ID=$(sed 's/Submitted batch job //g' <<< ${job1})
			printf "\n"
			printf "${GREEN}Job 1: Quality Trimming; ${qc} queued with jobid=${job1ID}.\n"
			printf "\n"
	# derep  - dependent on QC
	job2=$(sbatch --dependency=afterok:${job1ID} -J ${drJobName} -N ${drNodes} -n ${drTasks} -c ${drCPUsPerTask} -t ${drTime} --mem ${drMem} -o ${drLog} ${derep})
	job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
			printf "${GREEN}Job 2: Dereplication; ${derep} queued with jobid=${job2ID}.\n"
			printf "${GREEN}${derep} will begin after successfull completion of ${qc}.\n"
			printf "\n"
	# trinity phase 1 - dependent on derep
	job3=$(sbatch --dependency=afterok:${job2ID} -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
			printf "${GREEN}Job 3: Trinity part 1; ${trinP1} queued with jobid=${job3ID}.\n"
			printf "${GREEN}${trinP1} will begin after successfull completion of ${derep}.\n"
			printf "\n"
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			printf "${GREEN}Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}.\n"
			printf "${GREEN}${trinP2} will begin after successfull completion of ${trinP1}.\n"
			printf "\n"
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			printf "${GREEN}Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}.\n"
			printf "${GREEN}${trinP3} will begin after successfull completion of ${trinP2}.\n"
			printf "\n"
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			printf "${GREEN}Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}.\n"
			printf "${GREEN}${blastnTrin} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# diamond .daa - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}.\n"
			printf "${GREEN}${diamTrin} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# diamond .txt - dependent on trinity phase 3
	job8=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diam2Log} ${diamTrin2})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			printf "${GREEN}Job 8: Diamond blastX to nr; ${diamTrin2} queued with jobid=${job8ID}.\n"
			printf "${GREEN}${diamTrin2} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# calculate contig coverage
	job9=$(sbatch --dependency=afterok:${job5ID} -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
			printf "${GREEN}Job 9: Calc. contig coverage; ${concov} queued with jobid=${job9ID}.\n"
			printf "${GREEN}${concov} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"

	# outputting PIDs for easy reference
	printf "${ORANGE}${job1ID}"
	printf "${ORANGE}${job2ID}"	
	printf "${ORANGE}${job3ID}"		
	printf "${ORANGE}${job4ID}"	
	printf "${ORANGE}${job5ID}"		
	printf "${ORANGE}${job6ID}"		
	printf "${ORANGE}${job7ID}"
	printf "${ORANGE}${job8ID}"		
	printf "${ORANGE}${job9ID}"	
}
derep() {
	printf "\n"
	printf "${BLUE}Running work flow steps Dereplication, Trinity phase 1, Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage.\n"
	printf "\n"
	# derep  - dependent on QC
	job2=$(sbatch -J ${drJobName} -N ${drNodes} -n ${drTasks} -c ${drCPUsPerTask} -t ${drTime} --mem ${drMem} -o ${drLog} ${derep})
	job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
			printf "${GREEN}Job 2: Dereplication; ${derep} queued with jobid=${job2ID}.\n"
			printf "\n"
	# trinity phase 1 - dependent on derep
	job3=$(sbatch --dependency=afterok:${job2ID} -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
			printf "${GREEN}Job 3: Trinity part 1; ${trinP1} queued with jobid=${job3ID}.\n"
			printf "${GREEN}${trinP1} will begin after successfull completion of ${derep}.\n"
			printf "\n"
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			printf "${GREEN}Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}.\n"
			printf "${GREEN}${trinP2} will begin after successfull completion of ${trinP1}.\n"
			printf "\n"
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			printf "${GREEN}Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}.\n"
			printf "${GREEN}${trinP3} will begin after successfull completion of ${trinP2}.\n"
			printf "\n"
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			printf "${GREEN}Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}.\n"
			printf "${GREEN}${blastnTrin} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# diamond .daa - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}.\n"
			printf "${GREEN}${diamTrin} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# diamond .txt - dependent on trinity phase 3
	job8=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diam2Log} ${diamTrin2})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin2} queued with jobid=${job8ID}.\n"
			printf "${GREEN}${diamTrin2} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# calculate contig coverage
	job9=$(sbatch --dependency=afterok:${job5ID} -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
			printf "${GREEN}Job 8: Calc. contig coverage; ${concov} queued with jobid=${job9ID}.\n"
			printf "${GREEN}${concov} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"

	# outputting PIDs for easy reference
	printf "${ORANGE}${job2ID}"	
	printf "${ORANGE}${job3ID}"		
	printf "${ORANGE}${job4ID}"	
	printf "${ORANGE}${job5ID}"		
	printf "${ORANGE}${job6ID}"		
	printf "${ORANGE}${job7ID}"
	printf "${ORANGE}${job8ID}"		
	printf "${ORANGE}${job9ID}"	
}
trinityP1() {
	printf "\n"
	printf "${BLUE}Running work flow steps Trinity phase 1, Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage.\n"
	printf "\n"
	# trinity phase 1 - dependent on derep
	job3=$(sbatch -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
			printf "${GREEN}Job 3: Trinity part 1; ${trinP1} queued with jobid=${job3ID}.\n"
			printf "\n"
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			printf "${GREEN}Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}.\n"
			printf "${GREEN}${trinP2} will begin after successfull completion of ${trinP1}.\n"
			printf "\n"
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			printf "${GREEN}Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}.\n"
			printf "${GREEN}${trinP3} will begin after successfull completion of ${trinP2}.\n"
			printf "\n"
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			printf "${GREEN}Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}.\n"
			printf "${GREEN}${blastnTrin} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# diamond .daa - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}.\n"
			printf "${GREEN}${diamTrin} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# diamond .txt - dependent on trinity phase 3
	job8=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diam2Log} ${diamTrin2})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin2} queued with jobid=${job8ID}.\n"
			printf "${GREEN}${diamTrin2} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# calculate contig coverage
	job9=$(sbatch --dependency=afterok:${job5ID} -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
			printf "${GREEN}Job 8: Calc. contig coverage; ${concov} queued with jobid=${job9ID}.\n"
			printf "${GREEN}${concov} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"

	# outputting PIDs for easy reference	
	printf "${ORANGE}${job3ID}"		
	printf "${ORANGE}${job4ID}"	
	printf "${ORANGE}${job5ID}"		
	printf "${ORANGE}${job6ID}"		
	printf "${ORANGE}${job7ID}"
	printf "${ORANGE}${job8ID}"		
	printf "${ORANGE}${job9ID}"	
}
trinityP2() {
	printf "\n"
	printf "Running work flow steps Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage.\n"
	printf "\n"
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			printf "${GREEN}Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}.\n"
			printf "\n"
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			printf "${GREEN}Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}.\n"
			printf "${GREEN}${trinP3} will begin after successfull completion of ${trinP2}.\n"
			printf "\n"
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			printf "${GREEN}Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}.\n"
			printf "${GREEN}${blastnTrin} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# diamond .daa - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}.\n"
			printf "${GREEN}${diamTrin} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# diamond .txt - dependent on trinity phase 3
	job8=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diam2Log} ${diamTrin2})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin2} queued with jobid=${job8ID}.\n"
			printf "${GREEN}${diamTrin2} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# calculate contig coverage
	job9=$(sbatch --dependency=afterok:${job5ID} -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
			printf "${GREEN}Job 8: Calc. contig coverage; ${concov} queued with jobid=${job9ID}.\n"
			printf "${GREEN}${concov} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"

	# outputting PIDs for easy reference		
	printf "${ORANGE}${job4ID}"	
	printf "${ORANGE}${job5ID}"		
	printf "${ORANGE}${job6ID}"		
	printf "${ORANGE}${job7ID}"
	printf "${ORANGE}${job8ID}"		
	printf "${ORANGE}${job9ID}"	
}
trinityP3() {
	printf "\n"
	printf "${BLUE}unning work flow steps Trinity phase 3, blast, diamond, & contig coverage.\n"
	printf "\n"
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			printf "${GREEN}Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}.\n"
			printf "\n"
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			printf "${GREEN}Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}.\n"
			printf "${GREEN}${blastnTrin} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# diamond .daa - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}.\n"
			printf "${GREEN}${diamTrin} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# diamond .txt - dependent on trinity phase 3
	job8=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diam2Log} ${diamTrin2})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin2} queued with jobid=${job8ID}.\n"
			printf "${GREEN}${diamTrin2} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"
	# calculate contig coverage
	job9=$(sbatch --dependency=afterok:${job5ID} -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
			printf "${GREEN}Job 8: Calc. contig coverage; ${concov} queued with jobid=${job9ID}.\n"
			printf "${GREEN}${concov} will begin after successfull completion of ${trinP3}.\n"
			printf "\n"

	# outputting PIDs for easy reference	
	printf "${ORANGE}${job5ID}"		
	printf "${ORANGE}${job6ID}"		
	printf "${ORANGE}${job7ID}"
	printf "${ORANGE}${job8ID}"		
	printf "${ORANGE}${job9ID}"	
}
blast+diamond() {
	printf "\n"
	printf "${BLUE}Running work flow steps blast, diamond, & contig coverage.\n"
	printf "\n"
	# blast - dependent on trinity phase 3
	job6=$(sbatch -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			printf "${GREEN}Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}.\n"
			printf "\n"
	# diamond .daa - dependent on trinity phase 3
	job7=$(sbatch -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}.\n"
			printf "\n"
	# diamond .txt - dependent on trinity phase 3
	job8=$(sbatch -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diam2Log} ${diamTrin2})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin2} queued with jobid=${job8ID}.\n"
			printf "\n"
	# calculate contig coverage
	job9=$(sbatch -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
			printf "${GREEN}Job 8: Calc. contig coverage; ${concov} queued with jobid=${job9ID}.\n"
			printf "\n"

	# outputting PIDs for easy reference		
	printf "${ORANGE}${job6ID}"		
	printf "${ORANGE}${job7ID}"
	printf "${ORANGE}${job8ID}"		
	printf "${ORANGE}${job9ID}"	
}
onlyQC() {
	printf "\n"
	printf "${BLUE}Running work flow step: QC.\n"
	printf "\n"
	# QC - no dependency
	job1=$(sbatch -J ${qcJobName} -N ${qcNodes} -n ${qcTasks} -c ${qcCPUsPerTask} -t ${qcTime} --mem ${qcMem} -o ${qcLog} ${qc})
	job1ID=$(sed 's/Submitted batch job //g' <<< ${job1})
			printf "\n"
			printf "${GREEN}Job 1: Quality Trimming; ${qc} queued with jobid=${job1ID}.\n"
			printf "\n"
}
onlyDerep() {
	printf "\n"
	printf "${BLUE}Running work flow step: Dereplication.\n"
	printf "\n"
	# derep - no dependency
	job2=$(sbatch -J ${drJobName} -N ${drNodes} -n ${drTasks} -c ${drCPUsPerTask} -t ${drTime} --mem ${drMem} -o ${drLog} ${derep})
	job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
			printf "${GREEN}Job 2: Dereplication; ${derep} queued with jobid=${job2ID}.\n"
			printf "\n"
}
onlyTrinity() {
	printf "\n"
	printf "${BLUE}Running work flow steps: Trinity phase 1, Trinity phase 2, Trinity phase 3.\n"
	printf "\n"
	# trinity phase 1 - no dependency
	job3=$(sbatch -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
			printf "${GREEN}Job 3: Trinity part 1; ${trinP1} queued with jobid=${job3ID}.\n"
			printf "\n"
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			printf "Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}.\n"
			printf "${GREEN}${trinP2} will begin after successfull completion of ${trinP1}.\n"
			printf "\n"
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			printf "${GREEN}Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}.\n"
			printf "${trinP3} will begin after successfull completion of ${trinP2}.\n"
			printf "\n"
	# writing PIDs
	printf "${ORANGE}${job3ID}"
	printf "${ORANGE}${job4ID}"
	printf "${ORANGE}${job5ID}"
}
onlyTrinityP1() {
	printf "\n"
	printf "${BLUE}Running work flow step: Trinity phase 1.\n"
	printf "\n"
	# trinity phase 1 - no dependency
	job3=$(sbatch -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
			printf "${GREEN}${trinP1} will begin after successfull completion of ${derep}.\n"
			printf "\n"
}
onlyTrinityP2() {
	printf "\n"
	printf "${BLUE}Running work flow step: Trinity phase 2.\n"
	printf "\n"
	# trinity phase 2 - no dependency
	job4=$(sbatch -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			printf "${GREEN}${trinP2} will begin after successfull completion of ${trinP1}.\n"
			printf "\n"
}
onlyTrinityP3() {
	printf "\n"
	printf "${BLUE}Running work flow step: Trinity phase 3.\n"
	printf "\n"
	# trinity phase 3 - no dependency
	job5=$(sbatch -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			printf "${GREEN}Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}.\n"
			printf "\n"
}
onlyBlast() {
	printf "\n"
	printf "${BLUE}Running work flow step: BlastN.\n"
	printf "\n"
	# blast - no dependency
	job6=$(sbatch -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			printf "${GREEN}Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}.\n"
			printf "\n"
}
onlyDiamond() {
	printf "\n"
	printf "${BLUE}Running work flow step: Diamond (blastX).\n"
	printf "\n"
	# diamond - no dependency
	job7=$(sbatch -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
		printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}.\n"
		printf "\n"
	# diamond .txt - dependent on trinity phase 3
	job8=$(sbatch -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diam2Log} ${diamTrin2})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			printf "${GREEN}Job 7: Diamond blastX to nr; ${diamTrin2} queued with jobid=${job8ID}.\n"
			printf "\n"
	# writing PIDs
	printf "${ORANGE}${job7ID}"
	printf "${ORANGE}${job8ID}"

}
onlyConCov() {
	printf "\n"
	printf "${BLUE}Running work flow step: Contig Coverage (bwa + samtools).\n"
	printf "\n"
	job9=$(sbatch -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
		printf "${GREEN}Job 8: Calc. contig coverage; ${concov} queued with jobid=${job9ID}.\n"
		printf "\n"
}

# logic for running modules
if [[ "${module}" == "all" ]]; then
	all
elif [[ "${module}" == "derep" ]]; then
	derep
elif [[ "${module}" == "trinityP1" ]]; then
	trinityP1
elif [[ "${module}" == "trinityP2" ]]; then
	trinityP2
elif [[ "${module}" == "trinityP3" ]]; then
	trinityP3
elif [[ "${module}" == "blast+diamond" ]]; then
	blast+diamond
elif [[ "${module}" == "onlyQC" ]]; then
	onlyQC
elif [[ "${module}" == "onlyDerep" ]]; then
	onlyDerep
elif [[ "${module}" == "onlyTrinity" ]]; then
	onlyTrinity
elif [[ "${module}" == "onlyTrinityP1" ]]; then
	onlyTrinityP1
elif [[ "${module}" == "onlyTrinityP2" ]]; then
	onlyTrinityP2
elif [[ "${module}" == "onlyTrinityP3" ]]; then
	onlyTrinityP3
elif [[ "${module}" == "onlyBlast" ]]; then
	onlyBlast
elif [[ "${module}" == "onlyDiamond" ]]; then
	onlyDiamond
elif [[ "${module}" == "onlyConCov" ]]; then
	onlyConCov
else
	printf "${RED}-m <module> not an accepted value. Options are all, derep, trinityP1, trinityP2, trinityP3, blast+diamond, onlyQC, onlyDerep, onlyTrinirt, onlyTrinityP1, onlyTrinityP2, onlyTrinityP3, onlyBlast, onlyDiamond, onlyConCov"
	date
	exit 1
fi
date
