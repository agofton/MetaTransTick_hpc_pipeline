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
# derep
job2=$(sbatch --dependency=afterok:${job1ID} ${derep})
	job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
# trinity_P1
job3=$(sbatch --dependency=afterok:${job2ID} ${trin_P1})
	job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
# trinity_P2
job4=$(sbatch --dependency=afterok:${job3ID} ${trin_P2})
	job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
# trinity_P3
job5=$(sbatch --dependency=afterok:${job4ID} ${trin_P3})
	job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
# blast
job6=$(sbatch --dependency=afterok:${job5ID} ${blastn_trin})
	job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})		
# diamond
job7=$(sbatch --dependency=afterok:${job5ID} ${diam_trin})
	job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})		

	
# last step, converting .blast and .diamond files to .rma6 files for megan is done on vm not pearcey.

