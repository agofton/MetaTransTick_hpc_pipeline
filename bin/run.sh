#!/bin/bash

date
source slurmParams.txt
source script_vars.txt

helpMessage() {
	"""
	./run.sh usage: 
	-m <all, derep, trinityP1, trinityP2, trinityP3, blast+diamond, onlyqc, onlyDerep, onlyTrinity, onlyTrinityP1, onlyTrinityP2, onlyTrinityP3, onlyBlast, onlyDiamond, onlybt2_mapping> 
	-h [show this message]
	
	Running with -m flags "all", "derep",  "trinityP1", "trinityP2", "trinityP3", or "blast+diamond" will start the pipeline at that step and then run to the end.
	
	Running with -m flags "onlyqc", "onlyDerep",  "onlyTrinityP1", "onlyTrinityP2", "onlyTrinityP3", "onlyBlast", "onlyDiamond", or "onlybt2_mapping" will run only that module.
	"""
}

# Command line arguments
while getopts "hm:" OPTION
do
	case "${OPTION}"
		in
			h) helpMessage;;
	     	m) module=${OPTARG};;
			/?) echo "Invalid option: -$OPTARG" 1>&2
	esac
done
shift $((OPTIND - 1))

# define job options
all() {
	echo "Running work flow steps qc, Dereplication, Trinity phase 1, Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage." 
	# qc - no dependency
	job1=$(sbatch -J ${qc_jobName} \
				  -N ${qc_nodes} \
				  -n ${qc_tasks} \
				  -c ${qc_cpus} \
				  -t ${qc_time} \
				  --mem ${qc_mem} \
				  -o ${qc_log} \
				  ${qc})
		job1ID=$(sed 's/Submitted batch job //g' <<< ${job1})
		echo "Job 1: Quality Trimming ${qc} queued with jobid ${job1ID}."
		echo ""
	# derep  - dependent on qc
	job2=$(sbatch --dependency=afterok:${job1ID} \
				  -J ${derep_jobName} \
				  -N ${derep_nodes} \
				  -n ${derep_tasks} \
				  -c ${derep_cpus} \
				  -t ${derep_time} \
				  --mem ${derep_mem} \
				  -o ${derep_log} \
				  ${derep})
		job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
		echo "Job 2: Dereplication queued with jobid ${job2ID}."
		echo "${derep} will begin after successfull completion of ${qc}."
		echo ""
	# trinity phase 1 - dependent on derep
	job3=$(sbatch --dependency=afterok:${job2ID} \
				  -J ${trinP1_jobName} \
				  -N ${trinP1_nodes} \
				  -n ${trinP1_tasks} \
				  -t ${trinP1_time} \
				  --mem ${trinP1_mem} \
				  -o ${trinP1_log} \
				  ${trinP1})
		job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
		echo "Job 3: Trinity phase 1 queued with jobid ${job3ID}."
		echo "${trinP1} will begin after successfull completion of ${derep}."
		echo ""
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} \
				  -J ${trinP2_jobName} \
				  -N ${trinP2_nodes} \
				  -n ${trinP2_tasks} \
				  -c ${trinP2_cpus} \
				  -t ${trinP2_time} \
				  --mem ${trinP2_mem} \
				  -o ${trinP2_log} \
				  -a ${trinP2_array} \
				  ${trinP2})
		job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
		echo "Job 4: Trinity phase2 queued with jobid ${job4ID}."
		echo "${trinP2} will begin after successfull completion of ${trinP1}."
		echo ""
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} \
				  -J ${trinP3_jobName} \
				  -N ${trinP3_nodes} \
				  -n ${trinP3_tasks} \
				  -c ${trinP3_cpus} \
				  -t ${trinP3_time} \
				  --mem ${trinP3_mem} \
				  -o ${trinP3_log} \
				  ${trinP3})
		job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
		echo "Job 5: Trinity phase3 queued with jobid ${job5ID}."
		echo "${trinP3} will begin after successfull completion of ${trinP2}."
		echo ""
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${blastn_jobName} \
				  -N ${blastn_nodes} \
				  -n ${blastn_tasks} \
				  -c ${blastn_cpus} \
				  -t ${blastn_time} \
				  --mem ${blastn_mem} \
				  -o ${blastn_log} \
				  ${blastn})
		job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
		echo "Job 6: blastN to nt queued with jobid ${job6ID}."
		echo "${blastn} will begin after successfull completion of ${trinP3}."
		echo ""
	# diamond .daa - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${diamond_daa_jobName} \
				  -N ${diamond_daa_nodes} \
				  -n ${diamond_daa_tasks} \
				  -c ${diamond_daa_cpus} \
				  -t ${diamond_daa_time} \
				  --mem ${diamond_daa_mem} \
				  -o ${diamond_daa_log} \
				  ${diamond_daa})
		job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
		echo "Job 7: Diamond blastX to nr queued with jobid ${job7ID}."
		echo "${diamond_daa} will begin after successfull completion of ${trinP3}."
		echo ""
	# diamond .tab - dependent on trinity phase 3
	job8=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${diamond_tab_jobName} \
				  -N ${diamond_tab_nodes} \
				  -n ${diamond_tab_tasks} \
				  -c ${diamond_tab_cpus} \
				  -t ${diamond_tab_time} \
				  --mem ${diamond_tab_mem} \
				  -o ${diamond_tab_log} \
				  ${diamond_tab})
		job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Diamond blastX to nr queued with jobid ${job8ID}."
			echo "${diamond_tab} will begin after successfull completion of ${trinP3}."
			echo ""
	# calculate contig coverage
	job9=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${bt2_jobName} \
				  -N ${bt2_nodes} \
				  -n ${bt2_tasks} \
				  -c ${bt2_cpus} \
				  -t ${bt2_time} \
				  --mem ${bt2_mem} \
				  -o ${bt2_log} \
				  ${bt2_mapping})
		job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
		echo "Job 9: bowtie2 mapping queued with jobid ${job9ID}."
		echo "${bt2_mapping} will begin after successfull completion of ${trinP3}."
		echo""
	# outputting PIDs for easy reference
	echo "${job1ID}"
	echo "${job2ID}"	
	echo "${job3ID}"		
	echo "${job4ID}"	
	echo "${job5ID}"		
	echo "${job6ID}"		
	echo "${job7ID}"
	echo "${job8ID}"		
	echo "${job9ID}"	
}
################################################################################
derep() {
	echo "Running work flow steps Dereplication, Trinity phase 1, Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage."
	# derep 
	job2=$(sbatch -J ${derep_jobName} \
				  -N ${derep_nodes} \
				  -n ${derep_tasks} \
				  -c ${derep_cpus} \
				  -t ${derep_time} \
				  --mem ${derep_mem} \
				  -o ${derep_log} \
				  ${derep})
		job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
		echo "Job 2: Dereplication queued with jobid ${job2ID}."
		echo ""
	# trinity phase 1 - dependent on derep
	job3=$(sbatch --dependency=afterok:${job2ID} \
				  -J ${trinP1_jobName} \
				  -N ${trinP1_nodes} \
				  -n ${trinP1_tasks} \
				  -t ${trinP1_time} \
				  --mem ${trinP1_mem} \
				  -o ${trinP1_log} \
				  ${trinP1})
		job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
		echo "Job 3: Trinity phase 1 queued with jobid ${job3ID}."
		echo "${trinP1} will begin after successfull completion of ${derep}."
		echo ""
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} \
				  -J ${trinP2_jobName} \
				  -N ${trinP2_nodes} \
				  -n ${trinP2_tasks} \
				  -c ${trinP2_cpus} \
				  -t ${trinP2_time} \
				  --mem ${trinP2_mem} \
				  -o ${trinP2_log} \
				  -a ${trinP2_array} \
				  ${trinP2})
		job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
		echo "Job 4: Trinity phase2 queued with jobid ${job4ID}."
		echo "${trinP2} will begin after successfull completion of ${trinP1}."
		echo ""
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} \
				  -J ${trinP3_jobName} \
				  -N ${trinP3_nodes} \
				  -n ${trinP3_tasks} \
				  -c ${trinP3_cpus} \
				  -t ${trinP3_time} \
				  --mem ${trinP3_mem} \
				  -o ${trinP3_log} \
				  ${trinP3})
		job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
		echo "Job 5: Trinity phase3 queued with jobid ${job5ID}."
		echo "${trinP3} will begin after successfull completion of ${trinP2}."
		echo ""
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${blastn_jobName} \
				  -N ${blastn_nodes} \
				  -n ${blastn_tasks} \
				  -c ${blastn_cpus} \
				  -t ${blastn_time} \
				  --mem ${blastn_mem} \
				  -o ${blastn_log} \
				  ${blastn})
		job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
		echo "Job 6: blastN to nt queued with jobid ${job6ID}."
		echo "${blastn} will begin after successfull completion of ${trinP3}."
		echo ""
	# diamond .daa - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${diamond_daa_jobName} \
				  -N ${diamond_daa_nodes} \
				  -n ${diamond_daa_tasks} \
				  -c ${diamond_daa_cpus} \
				  -t ${diamond_daa_time} \
				  --mem ${diamond_daa_mem} \
				  -o ${diamond_daa_log} \
				  ${diamond_daa})
		job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
		echo "Job 7: Diamond blastX to nr queued with jobid ${job7ID}."
		echo "${diamond_daa} will begin after successfull completion of ${trinP3}."
		echo ""
	# diamond .tab - dependent on trinity phase 3
	job8=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${diamond_tab_jobName} \
				  -N ${diamond_tab_nodes} \
				  -n ${diamond_tab_tasks} \
				  -c ${diamond_tab_cpus} \
				  -t ${diamond_tab_time} \
				  --mem ${diamond_tab_mem} \
				  -o ${diamond_tab_log} \
				  ${diamond_tab})
		job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Diamond blastX to nr queued with jobid ${job8ID}."
			echo "${diamond_tab} will begin after successfull completion of ${trinP3}."
			echo ""
	# calculate contig coverage
	job9=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${bt2_jobName} \
				  -N ${bt2_nodes} \
				  -n ${bt2_tasks} \
				  -c ${bt2_cpus} \
				  -t ${bt2_time} \
				  --mem ${bt2_mem} \
				  -o ${bt2_log} \
				  ${bt2_mapping})
		job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
		echo "Job 9: bowtie2 mapping queued with jobid ${job9ID}."
		echo "${bt2_mapping} will begin after successfull completion of ${trinP3}."
		echo""	
	# outputting PIDs for easy reference
	echo "${job2ID}"	
	echo "${job3ID}"		
	echo "${job4ID}"	
	echo "${job5ID}"		
	echo "${job6ID}"		
	echo "${job7ID}"
	echo "${job8ID}"		
	echo "${job9ID}"
}
################################################################################
trinityP1() {
	echo "Running work flow steps Trinity phase 1, Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage."
	# trinity phase 1
	job3=$(sbatch -J ${trinP1_jobName} \
				  -N ${trinP1_nodes} \
				  -n ${trinP1_tasks} \
				  -t ${trinP1_time} \
				  --mem ${trinP1_mem} \
				  -o ${trinP1_log} \
				  ${trinP1})
		job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
		echo "Job 3: Trinity phase 1 queued with jobid ${job3ID}."
		echo ""
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} \
				  -J ${trinP2_jobName} \
				  -N ${trinP2_nodes} \
				  -n ${trinP2_tasks} \
				  -c ${trinP2_cpus} \
				  -t ${trinP2_time} \
				  --mem ${trinP2_mem} \
				  -o ${trinP2_log} \
				  -a ${trinP2_array} \
				  ${trinP2})
		job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
		echo "Job 4: Trinity phase2 queued with jobid ${job4ID}."
		echo "${trinP2} will begin after successfull completion of ${trinP1}."
		echo ""
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} \
				  -J ${trinP3_jobName} \
				  -N ${trinP3_nodes} \
				  -n ${trinP3_tasks} \
				  -c ${trinP3_cpus} \
				  -t ${trinP3_time} \
				  --mem ${trinP3_mem} \
				  -o ${trinP3_log} \
				  ${trinP3})
		job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
		echo "Job 5: Trinity phase3 queued with jobid ${job5ID}."
		echo "${trinP3} will begin after successfull completion of ${trinP2}."
		echo ""
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${blastn_jobName} \
				  -N ${blastn_nodes} \
				  -n ${blastn_tasks} \
				  -c ${blastn_cpus} \
				  -t ${blastn_time} \
				  --mem ${blastn_mem} \
				  -o ${blastn_log} \
				  ${blastn})
		job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
		echo "Job 6: blastN to nt queued with jobid ${job6ID}."
		echo "${blastn} will begin after successfull completion of ${trinP3}."
		echo ""
	# diamond .daa - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${diamond_daa_jobName} \
				  -N ${diamond_daa_nodes} \
				  -n ${diamond_daa_tasks} \
				  -c ${diamond_daa_cpus} \
				  -t ${diamond_daa_time} \
				  --mem ${diamond_daa_mem} \
				  -o ${diamond_daa_log} \
				  ${diamond_daa})
		job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
		echo "Job 7: Diamond blastX to nr queued with jobid ${job7ID}."
		echo "${diamond_daa} will begin after successfull completion of ${trinP3}."
		echo ""
	# diamond .tab - dependent on trinity phase 3
	job8=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${diamond_tab_jobName} \
				  -N ${diamond_tab_nodes} \
				  -n ${diamond_tab_tasks} \
				  -c ${diamond_tab_cpus} \
				  -t ${diamond_tab_time} \
				  --mem ${diamond_tab_mem} \
				  -o ${diamond_tab_log} \
				  ${diamond_tab})
		job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Diamond blastX to nr queued with jobid ${job8ID}."
			echo "${diamond_tab} will begin after successfull completion of ${trinP3}."
			echo ""
	# calculate contig coverage
	job9=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${bt2_jobName} \
				  -N ${bt2_nodes} \
				  -n ${bt2_tasks} \
				  -c ${bt2_cpus} \
				  -t ${bt2_time} \
				  --mem ${bt2_mem} \
				  -o ${bt2_log} \
				  ${bt2_mapping})
		job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
		echo "Job 9: bowtie2 mapping queued with jobid ${job9ID}."
		echo "${bt2_mapping} will begin after successfull completion of ${trinP3}."
		echo""
	# outputting PIDs for easy reference	
	echo "${job3ID}"		
	echo "${job4ID}"	
	echo "${job5ID}"		
	echo "${job6ID}"		
	echo "${job7ID}"
	echo "${job8ID}"		
	echo "${job9ID}"	
}
################################################################################
trinityP2() {
	echo "Running work flow steps Trinity phase 2, Trinity phase 3, blast, diamond, & contig coverage."
	# trinity phase 2
	job4=$(sbatch -J ${trinP2_jobName} \
				  -N ${trinP2_nodes} \
				  -n ${trinP2_tasks} \
				  -c ${trinP2_cpus} \
				  -t ${trinP2_time} \
				  --mem ${trinP2_mem} \
				  -o ${trinP2_log} \
				  -a ${trinP2_array} \
				  ${trinP2})
		job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
		echo "Job 4: Trinity phase2 queued with jobid ${job4ID}."
		echo ""
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} \
				  -J ${trinP3_jobName} \
				  -N ${trinP3_nodes} \
				  -n ${trinP3_tasks} \
				  -c ${trinP3_cpus} \
				  -t ${trinP3_time} \
				  --mem ${trinP3_mem} \
				  -o ${trinP3_log} \
				  ${trinP3})
		job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
		echo "Job 5: Trinity phase3 queued with jobid ${job5ID}."
		echo "${trinP3} will begin after successfull completion of ${trinP2}."
		echo ""
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${blastn_jobName} \
				  -N ${blastn_nodes} \
				  -n ${blastn_tasks} \
				  -c ${blastn_cpus} \
				  -t ${blastn_time} \
				  --mem ${blastn_mem} \
				  -o ${blastn_log} \
				  ${blastn})
		job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
		echo "Job 6: blastN to nt queued with jobid=${job6ID}."
		echo "${blastn} will begin after successfull completion of ${trinP3}."
		echo ""
	# diamond .daa - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${diamond_daa_jobName} \
				  -N ${diamond_daa_nodes} \
				  -n ${diamond_daa_tasks} \
				  -c ${diamond_daa_cpus} \
				  -t ${diamond_daa_time} \
				  --mem ${diamond_daa_mem} \
				  -o ${diamond_daa_log} \
				  ${diamond_daa})
		job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
		echo "Job 7: Diamond blastX to nr queued with jobid ${job7ID}."
		echo "${diamond_daa} will begin after successfull completion of ${trinP3}."
		echo ""
	# diamond .tab - dependent on trinity phase 3
	job8=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${diamond_tab_jobName} \
				  -N ${diamond_tab_nodes} \
				  -n ${diamond_tab_tasks} \
				  -c ${diamond_tab_cpus} \
				  -t ${diamond_tab_time} \
				  --mem ${diamond_tab_mem} \
				  -o ${diamond_tab_log} \
				  ${diamond_tab})
		job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Diamond blastX to nr queued with jobid ${job8ID}."
			echo "${diamond_tab} will begin after successfull completion of ${trinP3}."
			echo ""
	# calculate contig coverage
	job9=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${bt2_jobName} \
				  -N ${bt2_nodes} \
				  -n ${bt2_tasks} \
				  -c ${bt2_cpus} \
				  -t ${bt2_time} \
				  --mem ${bt2_mem} \
				  -o ${bt2_log} \
				  ${bt2_mapping})
		job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
		echo "Job 9: bowtie2 mapping queued with jobid ${job9ID}."
		echo "${bt2_mapping} will begin after successfull completion of ${trinP3}."
		echo""
	# outputting PIDs for easy reference			
	echo "${job4ID}"	
	echo "${job5ID}"		
	echo "${job6ID}"		
	echo "${job7ID}"
	echo "${job8ID}"		
	echo "${job9ID}"	
}
################################################################################
trinityP3() {
	echo "Running work flow steps Trinity phase 3, blast, diamond, & contig coverage."
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch -J ${trinP3_jobName} \
				  -N ${trinP3_nodes} \
				  -n ${trinP3_tasks} \
				  -c ${trinP3_cpus} \
				  -t ${trinP3_time} \
				  --mem ${trinP3_mem} \
				  -o ${trinP3_log} \
				  ${trinP3})
		job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
		echo "Job 5: Trinity phase3 queued with jobid ${job5ID}."
		echo ""
	# blast - dependent on trinity phase 3
	job6=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${blastn_jobName} \
				  -N ${blastn_nodes} \
				  -n ${blastn_tasks} \
				  -c ${blastn_cpus} \
				  -t ${blastn_time} \
				  --mem ${blastn_mem} \
				  -o ${blastn_log} \
				  ${blastn})
		job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
		echo "Job 6: blastN to nt queued with jobid ${job6ID}."
		echo "${blastn} will begin after successfull completion of ${trinP3}."
		echo ""
	# diamond .daa - dependent on trinity phase 3
	job7=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${diamond_daa_jobName} \
				  -N ${diamond_daa_nodes} \
				  -n ${diamond_daa_tasks} \
				  -c ${diamond_daa_cpus} \
				  -t ${diamond_daa_time} \
				  --mem ${diamond_daa_mem} \
				  -o ${diamond_daa_log} \
				  ${diamond_daa})
		job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
		echo "Job 7: Diamond blastX to nr queued with jobid ${job7ID}."
		echo "${diamond_daa} will begin after successfull completion of ${trinP3}."
		echo ""
	# diamond .tab - dependent on trinity phase 3
	job8=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${diamond_tab_jobName} \
				  -N ${diamond_tab_nodes} \
				  -n ${diamond_tab_tasks} \
				  -c ${diamond_tab_cpus} \
				  -t ${diamond_tab_time} \
				  --mem ${diamond_tab_mem} \
				  -o ${diamond_tab_log} \
				  ${diamond_tab})
		job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Diamond blastX to nr queued with jobid ${job8ID}."
			echo "${diamond_tab} will begin after successfull completion of ${trinP3}."
			echo ""
	# calculate contig coverage
	job9=$(sbatch --dependency=afterok:${job5ID} \
				  -J ${bt2_jobName} \
				  -N ${bt2_nodes} \
				  -n ${bt2_tasks} \
				  -c ${bt2_cpus} \
				  -t ${bt2_time} \
				  --mem ${bt2_mem} \
				  -o ${bt2_log} \
				  ${bt2_mapping})
		job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
		echo "Job 9: bowtie2 mapping queued with jobid ${job9ID}."
		echo "${bt2_mapping} will begin after successfull completion of ${trinP3}."
		echo""
	# outputting PIDs for easy reference				
	echo "${job5ID}"		
	echo "${job6ID}"		
	echo "${job7ID}"
	echo "${job8ID}"		
	echo "${job9ID}"
}
################################################################################
blast+diamond() {
	echo "Running work flow steps blast, diamond, & contig coverage."
	# blast
	job6=$(sbatch -J ${blastn_jobName} \
				  -N ${blastn_nodes} \
				  -n ${blastn_tasks} \
				  -c ${blastn_cpus} \
				  -t ${blastn_time} \
				  --mem ${blastn_mem} \
				  -o ${blastn_log} \
				  ${blastn})
		job6ID=$(sed 's/Submitted batch job //g' <<< ${job6})
		echo "Job 6: blastN to nt queued with jobid ${job6ID}."
		echo ""
	# diamond .daa
	job7=$(sbatch -J ${diamond_daa_jobName} \
				  -N ${diamond_daa_nodes} \
				  -n ${diamond_daa_tasks} \
				  -c ${diamond_daa_cpus} \
				  -t ${diamond_daa_time} \
				  --mem ${diamond_daa_mem} \
				  -o ${diamond_daa_log} \
				  ${diamond_daa})
		job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
		echo "Job 7: Diamond blastX to nr queued with jobid ${job7ID}."
		echo ""
	# diamond .tab - dependent on trinity phase 3
	job8=$(sbatch -J ${diamond_tab_jobName} \
				  -N ${diamond_tab_nodes} \
				  -n ${diamond_tab_tasks} \
				  -c ${diamond_tab_cpus} \
				  -t ${diamond_tab_time} \
				  --mem ${diamond_tab_mem} \
				  -o ${diamond_tab_log} \
				  ${diamond_tab})
		job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
			echo "Job 8: Diamond blastX to nr queued with jobid ${job8ID}."
			echo ""
	# calculate contig coverage
	job9=$(sbatch -J ${bt2_jobName} \
				  -N ${bt2_nodes} \
				  -n ${bt2_tasks} \
				  -c ${bt2_cpus} \
				  -t ${bt2_time} \
				  --mem ${bt2_mem} \
				  -o ${bt2_log} \
				  ${bt2_mapping})
		job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
		echo "Job 9: bowtie2 mapping queued with jobid ${job9ID}."
		echo""
	# outputting PIDs for easy reference						
	echo "${job6ID}"		
	echo "${job7ID}"
	echo "${job8ID}"		
	echo "${job9ID}"
}
################################################################################
onlyqc() {
	echo "Running work flow step: qc."
	# qc - no dependency
	job1=$(sbatch -J ${qc_jobName} \
				  -N ${qc_nodes} \
				  -n ${qc_tasks} \
				  -c ${qc_cpus} \
				  -t ${qc_time} \
				  --mem ${qc_mem} \
				  -o ${qc_log} \
				  ${qc})
		job1ID=$(sed 's/Submitted batch job //g' <<< ${job1})
		echo "Job 1: Quality Trimming queued with jobid ${job1ID}."
}
################################################################################
onlyDerep()	{
	echo "Running work flow step: Dereplication."
	# derep - no dependency
	job2=$(sbatch -J ${derep_jobName} \
				  -N ${derep_nodes} \
				  -n ${derep_tasks} \
				  -c ${derep_cpus} \
				  -t ${derep_time} \
				  --mem ${derep_mem} \
				  -o ${derep_log} \
				  ${derep})
		job2ID=$(sed 's/Submitted batch job //g' <<< ${job2})
		echo "Job 2: Dereplication queued with jobid ${job2ID}."	 
}
################################################################################
onlyTrinity() {
	echo "Running work flow steps: Trinity phase 1, Trinity phase 2, Trinity phase 3."
	# trinity phase 1 - no dependency
	job3=$(sbatch -J ${trinP1_jobName} \
				  -N ${trinP1_nodes} \
				  -n ${trinP1_tasks} \
				  -t ${trinP1_time} \
				  --mem ${trinP1_mem} \
				  -o ${trinP1_log} \
				  ${trinP1})
		job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
		echo "Job 3: Trinity phase 1 queued with jobid ${job3ID}."
	# trinity phase 2 - dependent on trinity phase 1
	job4=$(sbatch --dependency=afterok:${job3ID} -J ${trinP2_jobName} \
				  -N ${trinP2_nodes} \
				  -n ${trinP2_tasks} \
				  -c ${trinP2_cpus} \
				  -t ${trinP2_time} \
				  --mem ${trinP2_mem} \
				  -o ${trinP2_log} \
				  -a ${trinP2_array} \
				  ${trinP2})
		job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
		echo "Job 4: Trinity phase 2 queued with jobid ${job4ID}."
		echo "${trinP2} will begin after successfull completion of ${trinP1}."
	# trinity phase 3 - dependent on trinity phase 2
	job5=$(sbatch --dependency=afterok:${job4ID} -J ${trinP3_jobName} \
				  -N ${trinP3_nodes} \
				  -n ${trinP3_tasks} \
				  -c ${trinP3_cpus} \
				  -t ${trinP3_time} \
				  --mem ${trinP3_mem} \
				  -o ${trinP3_log} \
				  ${trinP3})
		job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
		echo "Job 5: Trinity phase 3 queued with jobid ${job5ID}."
		echo "${trinP3} will begin after successfull completion of ${trinP2}."
	# writing PIDs
	echo "${job3ID}"
	echo "${job4ID}"
	echo "${job5ID}"
}
################################################################################
onlyTrinityP1() {
	echo "Running work flow step: Trinity phase 1."
	# trinity phase 1 - no dependency
	job3=$(sbatch -J ${trinP1_jobName} \
				  -N ${trinP1_nodes} \
				  -n ${trinP1_tasks} \
				  -t ${trinP1_time} \
				  --mem ${trinP1_mem} \
				  -o ${trinP1_log} \
				  ${trinP1})
		job3ID=$(sed 's/Submitted batch job //g' <<< ${job3})
		echo "Trinity phase 1 queued with jobid ${job3ID}"
}
################################################################################
onlyTrinityP2() {
	echo "Running work flow step: Trinity phase 2."
	# trinity phase 2 - no dependency
	job4=$(sbatch -J ${trinP2_jobName} \
				  -N ${trinP2_nodes} \
				  -n ${trinP2_tasks} \
				  -c ${trinP2_cpus} \
				  -t ${trinP2_time} \
				  --mem ${trinP2_mem} \
				  -o ${trinP2_log} \
				  -a ${trinP2_array} \
				  ${trinP2})
		job4ID=$(sed 's/Submitted batch job //g' <<< ${job4})
		echo "Trinity phase 2 queued with jobid ${job4ID}"
}
################################################################################
onlyTrinityP3() {
	echo "Running work flow step: Trinity phase 3."
	# trinity phase 3 - no dependency
	job5=$(sbatch -J ${trinP3_jobName} \
				  -N ${trinP3_nodes} \
				  -n ${trinP3_tasks} \
				  -c ${trinP3_cpus} \
				  -t ${trinP3_time} \
				  --mem ${trinP3_mem} \
				  -o ${trinP3_log} \
				  ${trinP3})
		job5ID=$(sed 's/Submitted batch job //g' <<< ${job5})
		echo "Trinity phase 3 queued with jobid ${job5ID}."
}
################################################################################
onlyBlast() {
	echo "Running work flow step: BlastN."
	# blast - no dependency
	job6=$(sbatch -J ${blastn_jobName} \
				  -N ${blastn_nodes} \
				  -n ${blastn_tasks} \
				  -c ${blastn_cpus} \
				  -t ${blastn_time} \
				  --mem ${blastn_mem} \
				  -o ${blastn_log} \
				  ${blastn})
		job6ID=$(sed 's/Submitted batch job //g' <<< ${job6}) 
		echo "Job 6: blastN to nt queued with jobid ${job6ID}."
}
################################################################################
onlyDiamond() {
	echo "Running work flow step: Diamond (blastX)."
	# diamond - no dependency
	job7=$(sbatch -J ${diamond_daa_jobName} \
				  -N ${diamond_daa_nodes} \
				  -n ${diamond_daa_tasks} \
				  -c ${diamond_daa_cpus} \
				  -t ${diamond_daa_time} \
				  --mem ${diamond_daa_mem} \
				  -o ${diamond_daa_log} \
				  ${diamond_daa})
		job7ID=$(sed 's/Submitted batch job //g' <<< ${job7})
		echo "Job 7: Diamond blastX to nr queued with jobid ${job7ID}."
	# diamond .txt - dependent on trinity phase 3
	job8=$(sbatch -J ${diamond_tab_jobName} \
				  -N ${diamond_tab_nodes} \
				  -n ${diamond_tab_tasks} \
				  -c ${diamond_tab_cpus} \
				  -t ${diamond_tab_time} \
				  --mem ${diamond_tab_mem} \
				  -o ${diamond_tab_log} \
				  ${diamond_tab}) \
		job8ID=$(sed 's/Submitted batch job //g' <<< ${job8})
		echo "Job 7: Diamond blastX to nr queued with jobid ${job8ID}."
	# writing PIDs
	echo "${job7ID}"
	echo "${job8ID}"
}
################################################################################
onlybt2_mapping() {
	echo "Running work flow step: Contig Coverage (bowtie2 + samtools)."
	job9=$(sbatch -J ${bt2_jobName} \
				  -N ${bt2_nodes} \
				  -n ${bt2_tasks} \
				  -c ${bt2_cpus} \
				  -t ${bt2_time} \
				  --mem ${bt2_mem} \
				  -o ${bt2_log} \
				  ${bt2_mapping})
		job9ID=$(sed 's/Submitted batch job //g' <<< ${job9})
		echo "Job 8: bowtie2 mapping queued with jobid ${job9ID}."
}
################################################################################
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
elif [[ "${module}" == "onlyqc" ]]; then
	onlyqc
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
elif [[ "${module}" == "onlybt2_mapping" ]]; then
	onlybt2_mapping
else
	echo "-m <module> not an accepted value. Options are all, derep, trinityP1, trinityP2, trinityP3, blast+diamond, onlyqc, onlyDerep, onlyTrinirt, onlyTrinityP1, onlyTrinityP2, onlyTrinityP3, onlyBlast, onlyDiamond, onlybt2_mapping"
	date
	exit 1
fi
date
