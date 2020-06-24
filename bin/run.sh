#!/bin/bash
date

helpMessage() {
	echo "$0 usage: -s <all, derep, trinityP1, trinityP2, trinityP3, blast+diamond, onlyQC, onlyDerep, onlyTrinity, onlyTrinityP1, onlyTrinityP2, onlyTrinityP3, onlyBlast, onlyDiamond> -h [show this message]"
	echo 'Running with -s flags "all", "derep",  "trinityP1", "trinityP2", "trinityP3", or "blast+diamond" will start the pipeline at that step and then run to the end.'
	echo 'Running with -s flags begining with "only" with run only that module.'
}

flowControl() {
	if [ $? -ne 0 ]
	then
		echo $1; date; exit 1
	fi
}

source slurmParams.txt && flowControl "Cannot source variables from ./slurmParams.txt. File does not exits."

# job listing
qc=QC.SAMPLEID.slurm.sh
derep=derep.SAMPLEID.slurm.sh
trinP1=trinity_P1.SAMPLEID.slurm.sh
trinP2=trinity_P2.SAMPLEID.slurm.sh
trinP3=trinity_P3.SAMPLEID.slurm.sh
blastnTrin=blast_trin.SAMPLEID.slurm.sh
diamTrin=diam_trin.SAMPLEID.slurm.sh

# Command line arguments
while getopts "hs:" OPT
do
	case "${OPT}"
		in
	        h) helpMessage;;
	        s) steps=${OPTARG};;
	esac
done
shift $((OPTIND - 1))

if [ "${step}" == "all" ]
then
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

elif [[ "${step}" == "derep" ]]; then
	# derep - no dependency
	job2=$(sbatch -J ${drJobName} -N ${drNodes} -n ${drTasks} -c ${drCPUsPerTask} -t ${drTime} --mem ${drMem} -o ${drLog} ${derep})
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

elif [[ "${step}" == "trinityP1" ]]; then
	# trinity phase 1 - no dependency
	job3=$(sbatch -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
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

elif [[ "${step}" == "trinityP2" ]]; then
	# trinity phase 2 - no dependency
	job4=$(sbatch -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
		job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			echo "Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}."
			echo "${trinP2} will begin after successfull completion of ${trinP1}."
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

elif [[ "${step}" == "trinityP3" ]]; then
	# trinity phase 3 - no dependency
	job5=$(sbatch -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
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

elif [[ "${step}" == "blast+diamond" ]]; then
	# blast - no dependency
	job6=$(sbatch -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
		job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			echo "Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}."
			echo "${blastnTrin} will begin after successfull completion of ${trinP3}."
			echo ""
	# diamond - no dependency
	job7=$(sbatch -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
		job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			echo "Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}."
			echo "${diamTrin} will begin after successfull completion of ${trinP3}."
			echo ""

elif [[ "${step}" == "onlyQC" ]]; then
	# QC - no dependency
	job1=$(sbatch -J ${qcJobName} -N ${qcNodes} -n ${qcTasks} -c ${qcCPUsPerTask} -t ${qcTime} --mem ${qcMem} -o ${qcLog} ${qc})
		job1ID=$(sed 's/Submitted batch job //g' <<< ${job1})
			echo ""
			echo "Job 1: Quality Trimming; ${qc} queued with jobid=${job1ID}."
			echo ""

elif [[ "${step}" == "onlyDerep" ]]; then
	# derep - no dependency
	job2=$(sbatch -J ${drJobName} -N ${drNodes} -n ${drTasks} -c ${drCPUsPerTask} -t ${drTime} --mem ${drMem} -o ${drLog} ${derep})
		job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
			echo "Job 2: Dereplication; ${derep} queued with jobid=${job2ID}."
			echo "${derep} will begin after successfull completion of ${qc}."
			echo ""

elif [[ "${step}" == "onlyTrinity" ]]; then
	# trinity phase 1 - no dependency
	job3=$(sbatch -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
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

elif [[ "${step}" == "onlyTrinityP1" ]]; then
	# trinity phase 1 - no dependency
	job3=$(sbatch -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})
		job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
			echo "Job 3: Trinity part 1; ${trinP1} queued with jobid=${job3ID}."
			echo "${trinP1} will begin after successfull completion of ${derep}."
			echo ""

elif [[ "${step}" == "onlyTrinityP2" ]]; then
	# trinity phase 2 - no dependency
	job4=$(sbatch -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})
		job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
			echo "Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}."
			echo "${trinP2} will begin after successfull completion of ${trinP1}."
			echo ""

elif [[ "${step}" == "onlyTrinityP3" ]]; then
	# trinity phase 3 - no dependency
	job5=$(sbatch -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})
		job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
			echo "Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}."
			echo "${trinP3} will begin after successfull completion of ${trinP2}."
			echo ""

elif [[ "${step}" == "onlyBlast" ]]; then
	# blast - no dependency
	job6=$(sbatch -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
		job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
			echo "Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}."
			echo "${blastnTrin} will begin after successfull completion of ${trinP3}."
			echo ""

elif [[ "${step}" == "onlyDiamond" ]]; then
	# diamond - no dependency
	job7=$(sbatch -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})
		job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
			echo "Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}."
			echo "${diamTrin} will begin after successfull completion of ${trinP3}."
			echo ""

else
	echo "-s <step> not an accepted value. Options are all, QC, derep, trinityP1, trinityP2, trinityP3, blast+diamond, onlyBlast, onlyDiamond"
	date
	exit 1
fi
date
