#!/bin/bash
date

source slurmParams.txt

helpMessage() {
	echo "$0 usage: -m <all, derep, trinityP1, trinityP2, trinityP3, blast+diamond, onlyQC, onlyDerep, onlyTrinity, onlyTrinityP1, onlyTrinityP2, onlyTrinityP3, onlyBlast, onlyDiamond, onlyConCov> -h [show this message]"
	echo 'Running with -m flags "all", "derep",  "trinityP1", "trinityP2", "trinityP3", or "blast+diamond" will start the pipeline at that step and then run to the end.'
	echo 'Running with -m flags "onlyQC", "onlyDerep",  "onlyTrinityP1", "onlyTrinityP2", "onlyTrinityP3", "onlyBlast", "onlyDiamond", or "onlyConCov" will run only that module.'
}

flowControl() {
	if [ $? -ne 0 ]
	then
		echo $1; date; exit 1
	fi
}

# Command line arguments
while getopts "hm:" OPTION
do
	case "${OPTION}"
		in
	        h) helpMessage;;
	        m) module=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

# define job options
all() {
	echo ""
	echo "Running work flow steps QC, Dereplication, Trinity phase 1, Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage."
	echo ""
	# QC - no dependency
	job1=$(sbatch -J ${qcJobName} -N ${qcNodes} -n ${qcTasks} -c ${qcCPUsPerTask} -t ${qcTime} --mem ${qcMem} -o ${qcLog} ${qc})
	job1ID=$(sed 's/Submitted batch job //g' <<< ${job1})
			echo ""
			echo "Job 1: Quality Trimming; ${qc} queued with jobid=${job1ID}."
			echo ""
	# derep  - dependent on QC
	job2=$(sbatch --dependency=afterok:${job1ID} -J ${drJobName} -N ${drNodes} -n ${drTasks} -c ${drCPUsPerTask} -t ${drTime} --mem ${drMem} -o ${drLog} ${derep})
	job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
			echo "Job 2: Dereplication; ${derep} queued with jobid=${job2ID}."
			echo "${derep} will begin after successfull completion of ${qc}."
			echo ""
	# trinity phase 1 - dependent on derep
	job3=$(sbatch --dependency=afterok:${job2ID} -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
			echo "Job 3: Trinity part 1; ${trinP1} queued with jobid=${job3ID}."
			echo "${trinP1} will begin after successfull completion of ${derep}."
			echo ""
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			echo "Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}."
			echo "${trinP2} will begin after successfull completion of ${trinP1}."
			echo ""
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			echo "Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}."
			echo "${trinP3} will begin after successfull completion of ${trinP2}."
			echo ""
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			echo "Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}."
			echo "${blastnTrin} will begin after successfull completion of ${trinP3}."
			echo ""
	# diamond - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			echo "Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}."
			echo "${diamTrin} will begin after successfull completion of ${trinP3}."
			echo ""
	# calculate contig coverage
	job8=$(sbatch --dependency=afterok:${job5ID} -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Calc. contig coverage; ${concov} queued with jobid=${job8ID}."
			echo "${concov} will begin after successfull completion of ${trinP3}."
			echo ""

	# outputting PIDs for easy reference
	echo ${job1ID}
	echo ${job2ID}	
	echo ${job3ID}		
	echo ${job4ID}	
	echo ${job5ID}		
	echo ${job6ID}		
	echo ${job7ID}		
	echo ${job8ID}	
}

derep() {
	echo ""
	echo "Running work flow steps Dereplication, Trinity phase 1, Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage."
	echo ""
	# derep - no dependency
	job2=$(sbatch -J ${drJobName} -N ${drNodes} -n ${drTasks} -c ${drCPUsPerTask} -t ${drTime} --mem ${drMem} -o ${drLog} ${derep})
	job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
			echo "Job 2: Dereplication; ${derep} queued with jobid=${job2ID}."
			echo ""
	# trinity phase 1 - dependent on derep
	job3=$(sbatch --dependency=afterok:${job2ID} -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
			echo "Job 3: Trinity part 1; ${trinP1} queued with jobid=${job3ID}."
			echo "${trinP1} will begin after successfull completion of ${derep}."
			echo ""
	# trinity phase 2  - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			echo "Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}."
			echo "${trinP2} will begin after successfull completion of ${trinP1}."
			echo ""
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			echo "Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}."
			echo "${trinP3} will begin after successfull completion of ${trinP2}."
			echo ""
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			echo "Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}."
			echo "${blastnTrin} will begin after successfull completion of ${trinP3}."
			echo ""
	# diamond  - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			echo "Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}."
			echo "${diamTrin} will begin after successfull completion of ${trinP3}."
			echo ""
	# calculate contig coverage
	job8=$(sbatch --dependency=afterok:${job5ID} -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Calc. contig coverage; ${concov} queued with jobid=${job8ID}."
			echo "${concov} will begin after successfull completion of ${trinP3}."
			echo ""
	# outputting PIDs for easy reference
	echo ${job2ID}	
	echo ${job3ID}		
	echo ${job4ID}	
	echo ${job5ID}		
	echo ${job6ID}		
	echo ${job7ID}		
	echo ${job8ID}	
}

trinityP1() {
	echo ""
	echo "Running work flow steps Trinity phase 1, Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage."
	echo ""
	# trinity phase 1 - no dependency
	job3=$(sbatch -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
			echo "Job 3: Trinity part 1; ${trinP1} queued with jobid=${job3ID}."
			echo ""
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			echo "Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}."
			echo "${trinP2} will begin after successfull completion of ${trinP1}."
			echo ""
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			echo "Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}."
			echo "${trinP3} will begin after successfull completion of ${trinP2}."
			echo ""
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			echo "Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}."
			echo "${blastnTrin} will begin after successfull completion of ${trinP3}."
			echo ""
	# diamond - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			echo "Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}."
			echo "${diamTrin} will begin after successfull completion of ${trinP3}."
			echo ""
	job8=$(sbatch --dependency=afterok:${job5ID} -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Calc. contig coverage; ${concov} queued with jobid=${job8ID}."
			echo "${concov} will begin after successfull completion of ${trinP3}."
			echo ""
	# outputting PIDs for easy reference
	echo ${job3ID}		
	echo ${job4ID}	
	echo ${job5ID}		
	echo ${job6ID}		
	echo ${job7ID}		
	echo ${job8ID}	
}

trinityP2() {
	echo ""
	echo "Running work flow steps Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage."
	echo ""
	# trinity phase 2 - no dependency
	job4=$(sbatch -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			echo "Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}."
			echo ""
	# trinity phase 3- dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			echo "Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}."
			echo "${trinP3} will begin after successfull completion of ${trinP2}."
			echo ""
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			echo "Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}."
			echo "${blastnTrin} will begin after successfull completion of ${trinP3}."
			echo ""
	# diamond - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			echo "Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}."
			echo "${diamTrin} will begin after successfull completion of ${trinP3}."
			echo ""
	job8=$(sbatch --dependency=afterok:${job5ID} -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Calc. contig coverage; ${concov} queued with jobid=${job8ID}."
			echo "${concov} will begin after successfull completion of ${trinP3}."
			echo ""
	# outputting PIDs for easy reference
	echo ${job4ID}	
	echo ${job5ID}		
	echo ${job6ID}		
	echo ${job7ID}		
	echo ${job8ID}	
}

trinityP3() {
	echo ""
	echo "Running work flow steps Trinity phase 3, blast, diamond, & contig coverage."
	echo ""
	# trinity phase 3 - no dependency
	job5=$(sbatch -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			echo "Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}."
			echo ""
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			echo "Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}."
			echo "${blastnTrin} will begin after successfull completion of ${trinP3}."
			echo ""
	# diamond - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			echo "Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}."
			echo "${diamTrin} will begin after successfull completion of ${trinP3}."
			echo ""
	job8=$(sbatch --dependency=afterok:${job5ID} -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Calc. contig coverage; ${concov} queued with jobid=${job8ID}."
			echo "${concov} will begin after successfull completion of ${trinP3}."
			echo ""
	# outputting PIDs for easy reference
	echo ${job5ID}		
	echo ${job6ID}		
	echo ${job7ID}		
	echo ${job8ID}	
}

blast+diamond() {
	echo ""
	echo "Running work flow steps blast, diamond, & contig coverage."
	echo ""
	# blast - no dependency
	job6=$(sbatch -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			echo "Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}."
			echo ""
	# diamond - no dependency
	job7=$(sbatch -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			echo "Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}."
			echo ""
	job8=$(sbatch -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
	job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Calc. contig coverage; ${concov} queued with jobid=${job8ID}."
			echo ""
	
	# outputting PIDs for easy reference
	echo ${job6ID}		
	echo ${job7ID}
	echo ${job8ID}	
}

onlyQC() {
	echo ""
	echo "Running work flow step: QC."
	echo ""
	# QC - no dependency
	job1=$(sbatch -J ${qcJobName} -N ${qcNodes} -n ${qcTasks} -c ${qcCPUsPerTask} -t ${qcTime} --mem ${qcMem} -o ${qcLog} ${qc})
	job1ID=$(sed 's/Submitted batch job //g' <<< ${job1})
			echo ""
			echo "Job 1: Quality Trimming; ${qc} queued with jobid=${job1ID}."
			echo ""
}

onlyDerep() {
	echo ""
	echo "Running work flow step: Dereplication."
	echo ""
	# derep - no dependency
	job2=$(sbatch -J ${drJobName} -N ${drNodes} -n ${drTasks} -c ${drCPUsPerTask} -t ${drTime} --mem ${drMem} -o ${drLog} ${derep})
	job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
			echo "Job 2: Dereplication; ${derep} queued with jobid=${job2ID}."
			echo ""
}

onlyTrinity() {
	echo ""
	echo "Running work flow steps: Trinity phase 1, Trinity phase 2, Trinity phase 3."
	echo ""
	# trinity phase 1 - no dependency
	job3=$(sbatch -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
			echo "Job 3: Trinity part 1; ${trinP1} queued with jobid=${job3ID}."
			echo ""
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			echo "Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}."
			echo "${trinP2} will begin after successfull completion of ${trinP1}."
			echo ""
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			echo "Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}."
			echo "${trinP3} will begin after successfull completion of ${trinP2}."
			echo ""
	# writing PIDs
	echo ${job3ID}
	echo ${job4ID}
	echo ${job5ID}
}

onlyTrinityP1() {
	echo ""
	echo "Running work flow step: Trinity phase 1."
	echo ""
	# trinity phase 1 - no dependency
	job3=$(sbatch -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
			echo "${trinP1} will begin after successfull completion of ${derep}."
			echo ""
}

onlyTrinityP2() {
	echo ""
	echo "Running work flow step: Trinity phase 2."
	echo ""
	# trinity phase 2 - no dependency
	job4=$(sbatch -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			echo "${trinP2} will begin after successfull completion of ${trinP1}."
			echo ""
}

onlyTrinityP3() {
	echo ""
	echo "Running work flow step: Trinity phase 3."
	echo ""
	# trinity phase 3 - no dependency
	job5=$(sbatch -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			echo "Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}."
			echo ""
}

onlyBlast() {
	echo ""
	echo "Running work flow step: BlastN."
	echo ""
	# blast - no dependency
	job6=$(sbatch -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			echo "Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}."
			echo ""
}

onlyDiamond() {
	echo ""
	echo "Running work flow step: Diamond (blastX)."
	echo ""
	# diamond - no dependency
	job7=$(sbatch -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
	echo "Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}."
	echo ""
}

onlyConCov() {
	echo ""
	echo "Running work flow step: Contig Coverage (bwa + samtools)."
	echo ""
job8=$(sbatch -J ${concovJobName} -N ${concovNodes} -n ${concovTasks} -c ${concovCPUsPerTask} -t ${concovTime} --mem ${concovMem} -o ${concovLog} ${concov})
job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
		echo "Job 8: Calc. contig coverage; ${concov} queued with jobid=${job8ID}."
		echo ""
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
	echo '-m <module> not an accepted value. Options are "all", "derep", "trinityP1", "trinityP2", "trinityP3", "blast+diamond", "onlyQC", "onlyDerep", "onlyTrinirt", "onlyTrinityP1", "onlyTrinityP2", "onlyTrinityP3", "onlyBlast", "onlyDiamond", "onlyConCov"'
	date
	exit 1
fi
date
