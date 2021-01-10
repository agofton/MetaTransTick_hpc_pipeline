#!/usr/bin/python3

# Alexander W. Gofton, CSIRO, Jan 2021 alexander.gofton@csiro.au; alexander.gofton@gmail.com

import pandas as pd
import os
import sys
import argparse

parser = argparse.ArgumentParser(description = "Extracts contigs of a given taxon, maps all reads to those contigs using bwa mem, ")
parser.add_argument("-l", "--lca_sum_file", help = "lca summary file. Output of contig_lca_sum.py", type = str, required = True)
parser.add_argument("-t", "--taxon", help = "Taxon to filter contigs by.", type = str, required = True)
parser.add_argument("-s", "--sampleID", help = "sample ID prefix.", type = str, required = True)
parser.add_argument("-z", "--trin_fasta", help = "trinity fasta file.", type = str, required = True)
parser.add_argument("-f", "--R1", help = "R1.fastq for mapping.", type = str, required = True)
parser.add_argument("-r", "--R2", help = "R2.fastq for mapping.", type = str, required = True)
parser.add_argument("-b", "--trinity_bam", help = "bam file of whole trinity assembly", type = str, required = True)
parser.add_argument("-o", "--out_file", help = "output.txt", type = str, required = True)
parser.parse_args()
args = parser.parse_args()

# Filter and extract contigs and based on MEGAN LCA
# Get contig IDs
seqIDs = args.sampleID + '.' + args.taxon + '_contigs.seqIDs.txt'
fasta = args.sampleID + '.' + args.taxon + '_contigs.fasta'
prefix = args.sampleID + '.' + args.taxon + '_contigs'
bam = args.sampleID + '.' + args.taxon + '_contigs.bam'
flagstat = args.sampleID + '.' + args.taxon + '_contigs.flagstat.txt'

print("\n")

os.system('grep ' + '"' + args.taxon + '" ' + args.lca_sum_file + " | awk '{print $1}' > " + seqIDs)
print(str(args.taxon) + " read IDs found. Extracting contigs...")

os.system('usearch9.2_linux64 -fastx_getseqs ' + args.trin_fasta + ' -labels ' + seqIDs + ' -fastaout ' + fasta)
print("\n")

print("Indexing " + str(args.taxon) + " contigs...")
os.system('bwa index -p ' + prefix + ' ' + fasta)
print("\n")

print("Mapping reads to " + str(args.taxon) + " with bwa mem...")
os.system('bwa mem -t 8 -v 0 ' + prefix + ' ' + args.R1 + ' ' + args.R2 + ' | samtools sort > ' + bam)
print("Read mapping complete.")
print("\n")


# Get number of taxon contigs
fasta_file = open(args.trin_fasta)
n = 0
for line in fasta_file:
	if line.startswith(">"):
		n = n + 1
fasta_file.close()

# Get number of reads mapped to taxon contigs
os.system('~/samtools-1.10/samtools flagstat -@ 8 ' + bam + ' > ' + flagstat)
df = pd.read_csv(flagstat, sep = '\s+', header = None, names=["col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11"])
num_mapped = (df.iloc[0, 0] - df.iloc[1, 0] - df.iloc[2, 0]) - (df.iloc[0, 0] - df.iloc[4, 0]) 
#result = (num_mapped / num_reads) * 100
#rounded_result = round(result, 2)

# Get number of reads mapped to whole assembly
trin_flagstat = flagstat.tmp
os.system('~/samtools-1.10/samtools flagstat -@ 8 ' + args.trinity_bam + ' ' + '>' + ' ' + trin_flagstat)
df2 = pd.read_csv(flagstat.tmp, sep = '\s+', header = None, names=["col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11"])
num_mapped2 = (df2.iloc[0, 0] - df2.iloc[1, 0] - df2.iloc[2, 0]) - (df2.iloc[0, 0] - df2.iloc[4, 0]) 
os.system('rm -f ' + trin_flagstat)

# Calc pct
pct = (num_mapped / num_mapped2) * 100
rounded_pct = round(pct, 2)

# Print results
orig_stdout = sys.stdout  # Save reference to orignal stdout (print to screen)
with open(args.out_file, "w") as file:
	sys.stdout = file  # changes that standard output to the file created above
	print('#_reads_to_assembly' + "\t" + '#_' + str(args.taxon) + '_contigs' + "\t" + '#_reads_to_' + str(args.taxon) + "\t" +  'pct_reads_to_' + str(args.taxon))
	print(str(num_mapped2) + "\t" + str(n) + "\t" + str(num_mapped) + "\t" + str(rounded_pct) + '%')
	sys.stdout = orig_stdout  # revert std out back to orig
