#!/bin/bash

flowControl() {
	if [ $? -ne 0 ]
	then
		echo $1; date; exit 1
	else
		echo $2
	fi
}

source slurmParams.txt
flowControl "Cannot source variables from ./slurmParams.txt. File does not exits." " "

echo ${qcJobName}


# job listing
qc=QC.SAMPLEID.slurm.sh
derep=derep.SAMPLEID.slurm.sh
trinP1=trinity_P1.SAMPLEID.slurm.sh
trinP2=trinity_P2.SAMPLEID.slurm.sh
trinP3=trinity_P3.SAMPLEID.slurm.sh
blastnTrin=blast_trin.SAMPLEID.slurm.sh
diamTrin=diam_trin.SAMPLEID.slurm.sh

# QC - no dependencies, start asap.
job1=$(sbatch -J ${qcJobName} -N ${qcNodes} -n ${qcTasks} -c ${qcCPUsPerTask} -t ${qcTime} --mem ${qcMem} -o ${qcLog} ${qc})

	job1ID=$(sed 's/Submitted batch job //g' <<< ${job1})
	
	echo ""
	echo "Job 1: Quality Trimming; ${qc} queued with jobid=${job1ID}."
	echo ""

# derep
job2=$(sbatch --dependency=afterok:${job1ID} -J ${drJobName} -N ${drNodes} -n ${drTasks} -c ${drCPUsPerTask} -t ${drTime} --mem ${drMem} -o ${drLog} ${derep})

	job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
	
	echo "Job 2: Dereplication; ${derep} queued with jobid=${job2ID}."
	echo "${derep} will begin after successfull completion of ${qc}."
	echo ""

# trinity phase 1
job3=$(sbatch --dependency=afterok:${job2ID} -J ${t1JobName} -N ${t1Nodes} -n ${t1Tasks} -t ${t1Time} --mem ${t1Mem} -o ${t1Log} ${trinP1})

	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
	
	echo "Job 3: Trinity part 1; ${trinP1} queued with jobid=${job3ID}."
	echo "${trinP1} will begin after successfull completion of ${derep}."
	echo ""

# trinity phase 2
job4=$(sbatch --dependency=afterok:${job3ID} -J ${t2JobName} -N ${t2Nodes} -n ${t2Tasks} -c ${t2CPUsPerTask} -t ${t2Time} --mem ${t2Mem} -o ${t2Log} -a ${t2Arrays} ${trinP2})

	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})

	echo "Job 4: Trinity part 2; ${trinP2} queued with jobid=${job4ID}."
	echo "${trinP2} will begin after successfull completion of ${trinP1}."
	echo ""
	
# trinity phase 3
job5=$(sbatch --dependency=afterok:${job4ID} -J ${t3JobName} -N ${t3Nodes} -n ${t3Tasks} -c ${t3CPUsPerTask} -t ${t3Time} --mem ${t3Mem} -o ${t3Log} ${trinP3})

	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})

	echo "Job 5: Trinity part 3 (final); ${trinP3} queued with jobid=${job5ID}."
	echo "${trinP3} will begin after successfull completion of ${trinP2}."
	echo ""

# blast
job6=$(sbatch --dependency=afterok:${job5ID} -J ${blastJobName} -N ${blastNodes} -n ${blastTasks} -c ${blastCPUsPerTask} -t ${blastTime} --mem ${blastMem} -o ${blastLog} ${blastnTrin})
	
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})		

	echo "Job 6: blastN to nt; ${blastnTrin} queued with jobid=${job6ID}."
	echo "${blastnTrin} will begin after successfull completion of ${trinP3}."
	echo ""

# diamond
job7=$(sbatch --dependency=afterok:${job5ID} -J ${diamJobName} -N ${diamNodes} -n ${diamTasks} -c ${diamCPUsPerTask} -t ${diamTime} --mem ${diamMem} -o ${diamLog} ${diamTrin})

	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})		
	
	echo "Job 7: Diamond blastX to nr; ${diamTrin} queued with jobid=${job7ID}."
	echo "${diamTrin} will begin after successfull completion of ${trinP3}."
	echo ""

