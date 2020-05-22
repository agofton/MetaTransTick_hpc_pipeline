#!/bin/bash

# job listing
qc=QC.SAMPLEID.slurm.sh
derep=derep.SAMPLEID.slurm.sh
trinity=trinity.SAMPLEID.slurm.sh
trin_P1=trinity_P1.SAMPLEID.slurm.sh
trin_P2=trinity_P2.SAMPLEID.slurm.sh
trin_P3=trinity_P3.SAMPLEID.slurm.sh
spades=spades.SAMPLEID.slurm.sh
blastn_trin=blast_trin.SAMPLEID.slurm.sh
blastn_spades=blast_spades.SAMPLEID.slurm.sh
diam_trin=diam_trin.SAMPLEID.slurm.sh
diam_spades=diam_spades.SAMPLEID.slurm.sh

# queue jobs with dependencies
# QC
job1=$(sbatch ${qc})
job1ID=$(sed 's/Submitted batch job //g' <<< ${job1})
	
	echo ""
	echo "Job 1: Quality Trimming; ${qc} queued with jobid=${Job1ID}."
	echo ""

# derep
job2=$(sbatch --dependency=afterok:${job1ID} ${derep})
job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
	
	echo "Job 2: Dereplication; ${derep} queued with jobid=${job2ID}."
	echo "${derep} will begin after successfull completion of ${qc}."
	echo ""

# trinity_P1
job3=$(sbatch --dependency=afterok:${job2ID} ${trin_P1})
job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
	
	echo "Job 3: Trinity part 1; ${trin_P1} queued with jobid=${job3ID}."
	echo "${trin_P1} will begin after successfull completion of ${derep}."
	echo ""

# trinity_P2
job4=$(sbatch --dependency=afterok:${job3ID} ${trin_P2})
job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})

	echo "Job 4: Trinity part 2; ${trin_P2} queued with jobid=${job4ID}."
	echo "${trin_P2} will begin after successfull completion of ${trin_P1}."
	echo ""
	
# trinity_P3
job5=$(sbatch --dependency=afterok:${job4ID} ${trin_P3})
job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})

echo "Job 5: Trinity part 3 (final); ${trin_P3} queued with jobid=${job5ID}."
	echo "${trin_P3} will begin after successfull completion of ${trin_P2}."
	echo ""

# blast
job6=$(sbatch --dependency=afterok:${job5ID} ${blastn_trin})
job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})		

	echo "Job 6: blastN to nt; ${blastn_trin} queued with jobid=${job6ID}."
	echo "${blastn_trin} will begin after successfull completion of ${trin_P3}."
	echo ""

# diamond
job7=$(sbatch --dependency=afterok:${job5ID} ${diam_trin})
job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})		
	
	echo "Job 7: Diamond blastX to nr; ${diam_trin} queued with jobid=${job7ID}."
	echo "${diam_trin} will begin after successfull completion of ${trin_P3}."
	echo ""

# last step, converting .blast and .diamond files to .rma6 files for megan is done on vm not pearcey.

