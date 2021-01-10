#!/bin/bash

# Example metatrascriptomics workflow script with all parameters
# Actual analysis was performed on a HPC cluster through SLURM

# FastQC raw data - fastqc version 0.11.8
fastqc --outdir ${qcDir}/fastqc_raw --format fastq --threads ${SLURM_CPUS_PER_TASK} ${in1} ${in2}

# Remove adapters and distal bases - cutadapt version 2.8 in python version 3.6.1
cutadapt -j ${SLURM_CPUS_PER_TASK} -a AGATCGGAAGAG -A AGATCGGAAGAG --no-indels -O 4 -m 30 -o ${cutOut1} -p ${cutOut2} ${in1} ${in2}

# Trimmomatic version 0.38
trimmomatic PE -threads ${SLURM_CPUS_PER_TASK} -summary ${trimSum} ${cutOut1} ${cutOut2} ${trimPEr1} ${trimUPr1} ${trimPEr2} ${trimUPr2} SLIDINGWINDOW:5:15 MINLEN:30

# FastQC QC data - fastqc version 0.11.8
fastqc --outdir ${qcDir}/fastqc_QC --format fastq --threads ${SLURM_CPUS_PER_TASK} ${trimPEr1} ${trimPEr2} ${trimUPr1} ${trimUPr2}

# Removes only optical duplicates dedupe.sh in bbtools version 38.37
dedupe.sh in1=${trimPEr1} in2=${trimPEr2} out=${int} ac=f threads=${SLURM_CPUS_PER_TASK}
# Converts intleaved to paired
reformat.sh in=${int} out1=${derepOut1} out2=${derepOut2}

# FastQC dereplicated reads data - fastqc version 0.11.8
fastqc --outdir ${derepFQCoutDir} --format fastq --threads ${SLURM_CPUS_PER_TASK} ${derepOut1} ${derepOut2}

# Trinity denovo assembly with version 2.11.0
Trinity --seqType fq --max_memory ${maxMemGB} --left ${derepOut1} --right ${derepOut2} --CPU ${SLURM_NTASKS} --min_contig_length 200 --output ${trinOutDir} --no_normalize_reads

# Map reads back to assembly with bwa version 0.7.17 and samtools version 1.10.0
bwa index -p ${bwaindex} ${trinFasta}
bwa mem -t ${SLURM_CPUS_PER_TASK} ${conCovOutDir}/${bwaindex} ${derepOut1} ${derepOut2} | samtools sort > ${bamOut}
samtools coverage ${bamOut} > ${covSum}

# Calculate Transcripts Per Million (TPM) with custom R stript
./TPM_calc.R ${covSum}

# Calculate mapping rate with custom python script and samtools version 1.10
./mapping_percentage.py -i ${bamOut} -f {flagstat_out} -t ${threads}

# Blastn all transcripts to ncbi nt with megablast (blast+ version 2.11.0)
blastn -task megablast -query ${trinIn} -db ${database_nt} -out ${archive} -strand both -num_threads ${SLURM_CPUS_PER_TASK} -outfmt 11 -max_target_seqs 100 -max_hsps 5 -evalue 1e-10
blast_formatter -archive ${archive} -outfmt 0 -out ${outfmt0}
blast_formatter -archive ${archive} -out ${outfmt6} -outfmt '6 qseqid ssciname saccver stitle staxid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
awk '!a[$1]++' ${outfmt6} > ${blastntopHits} # Takes top hit only for each transcript (from tabular output)

# Diamond (blastx) all transcripts to ncbi nr (Diamond version 0.9.28)
diamond blastx --threads ${SLURM_CPUS_PER_TASK} --db ${database_nr} --out ${diamOut} --outfmt 100 --query ${trinIn} --strand both --max-target-seqs 100 --evalue 0.000000001 --sensitive --min-orf 80
diamond blastx --threads ${SLURM_CPUS_PER_TASK} --db ${database_nr} --out ${diamTabOut} --outfmt 6 qseqid sscinames sseqid stitle staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore --query ${trinIn} --strand both --max-target-seqs 100 --evalue 0.000000001 --sensitive --min-orf 80
awk '!a[$1]++' ${diamTabOut} > ${diamTopHits} # Takes top hit only for each transcript (from tabular output)

# MEGAN megan v - relaxed
~/megan/tools/blast2rma -i ${blast_in} -r ${reads_in} -o ${rma_out} -f BlastText -bm BlastN -c -m 5000 -ms 200 -me 0.0000000001 -mpi 70 -top 10 -sup 2 -mrc 50.0 -alg weighted -ram readCount -a2t nucl_acc2tax-Jul2019.abin -v
~/megan/tools/daa2rma -i ${blast_in} -o ${rma_out} -c -m 5000 -me 0.0000000001 -ms 200 -mpi 70 -top 10 -sup 2 -mrc 0 -alg weighted -ram readCount -a2t prot_acc2tax-Jul2019X1.abin -v

# .rma6 files were then opened in Megan 6.18.10, and LCA (read name to taxon path .csv) files were exported.

# Generating megan + blast + diamond comparison files with custom python script
./contig_lca_sum.py --lca_bn blastn_lca.csv --lca_bx blastx_lca.csv --bn ${blastntopHits} --bx ${diamTopHits} --tpm ${covSum} --out output.txt

# Extracting host (tick reads), read mapping, and generating stats
grep "Ixodida" output.txt | awk '{print $1}' > host_contigs.seqIDs
usearch9.2 -fastx_getseqs ${trinFasta} -labels host_contigs.seqIDs -fastaout host_contigs.fasta
bwa index -p host_contigs host_contigs.fasta 
bwa mem -t ${threads} host_contigs ${derepOut1} ${derepOut2} | samtools sort > ${host_bam}
./mapping_percentage.py -i ${host_bam} -f {flagstat_out} -t ${threads}

# Extracting viral (tick reads), read mapping, and generating stats
grep "Virus" output.txt | awk '{print $1}' > viral_contigs.seqIDs
usearch9.2 -fastx_getseqs ${trinFasta} -labels viral_contigs.seqIDs -fastaout viral_contigs.fasta
bwa index -p viral_contigs viral_contigs.fasta 
bwa mem -t ${threads} viral_contigs ${derepOut1} ${derepOut2} | samtools sort > ${viral_bam}
./mapping_percentage.py -i ${viral_bam} -f {flagstat_out} -t ${threads}
