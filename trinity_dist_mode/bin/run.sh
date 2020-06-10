#!/bin/bash

# job listing
trin_P1=trinity_P1.SAMPLEID.slurm.sh
trin_P2=trinity_P2.SAMPLEID.slurm.sh
trin_P3=trinity_P3.SAMPLEID.slurm.sh

# trinity_P1
job1=$(sbatch ${trin_P1})
job1ID=$(sed 's/Submitted batch job //g' <<< ${job1})
	
	echo "Job 1: Trinity part 1; ${trin_P1} queued with jobid=${job1ID}."
	echo ""

# trinity_P2
job2=$(sbatch --dependency=afterok:${job1ID} ${trin_P2})
job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})

	echo "Job 2: Trinity part 2; ${trin_P2} queued with jobid=${job2ID}."
	echo "${trin_P2} will begin after successfull completion of ${trin_P1}."
	echo ""
	
# trinity_P3
job3=$(sbatch --dependency=afterok:${job2ID} ${trin_P3})
job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})

	echo "Job 3: Trinity part 3 (final); ${trin_P3} queued with jobid=${job3ID}."
	echo "${trin_P3} will begin after successfull completion of ${trin_P2}."
	echo ""


# last step, converting .blast and .diamond files to .rma6 files for megan is done on vm not pearcey.

