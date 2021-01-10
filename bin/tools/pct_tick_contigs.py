#!/usr/bin/python3

# Alexander W. Gofton, CSIRO, Jan 2021 alexander.gofton@csiro.au; alexander.gofton@gmail.com

import pandas as pd
import os
import argparse

# Init argparse
parser = argparse.ArgumentParser(description="Calculates the number & pct of contigs in a given taxon grouping. Run this after running extract_and_map.sh")
parser.add_argument("-l", "--lca_sum", help = "lca summary file (output of contig_lca_sum.py)", type = str, required = True)
parser.add_argument("-t", "--trin_fasta", help = "trinity contigs fasta file", type = str, required = True)
parser.add_argument("-i", "--host_fasta", help = "host contigs fasta file", type = str, required = True)
parser.parse_args()
args = parser.parse_args()

# Count number of assigned contigs (from blastn or diamond)
file = open(args.lca_sum)
num_ass = len(file.readlines())
file.close()

# Count total number of transcripts in trinity.fasta
fasta = open(args.trin_fasta)
n = 0
for line in fasta:
	if line.startswith(">"):
		n = n + 1
fasta.close()

# Count total number of tick transcripts in host_contigs.fasta
fasta = open(args.host_fasta)
x = 0
for line in fasta:
	if line.startswith(">"):
		x = x + 1
fasta.close()

# Calc %
host_pct = (x / num_ass) * 100
rounded_host_pct = round(host_pct, 2)

ass_pct = (num_ass / n) * 100
rounded_ass_pct = round(ass_pct, 2)

# Print result
print(args.host_fasta)
print("Total" + "\t" + "num_assigned" + "\t" + "pct_assigned" + "\t" + "num_host" + "\t" + "pct_host")
print(str(n) + "\t" + str(num_ass) + "\t" + str(rounded_ass_pct) + "%" + "\t" + str(x) + "\t" + str(rounded_host_pct) + "%")

