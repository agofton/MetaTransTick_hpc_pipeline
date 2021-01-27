#!/usr/bin/python3

# Alexander W. Gofton, CSIRO, Jan 2021 alexander.gofton@csiro.au; alexander.gofton@gmail.com

import pandas as pd
import os
import sys
import argparse

####
# Init. argparse
parser = argparse.ArgumentParser(description="Calculates the number and percentage of reads mapped to a given taxon based on the flagstat.txt outputs of samtools flagstat.txt files.")
parser.add_argument("-i", "--trin_flagstat", help = "samtools flagstat output of trinity_bwa.bam", type = str, required = True)
parser.add_argument("-f", "--tax_flagstat", help = "samtools flagstat output of taxon_bwa.bam", type = str, required = True)
parser.add_argument("-t", "--threads", help = "Num threads", type = str, required = True)
parser.parse_args()
args = parser.parse_args()
####

# Calc num of reads mapped to whole assembly
df = pd.read_csv(args.trin_flagstat, sep = '\s+', header = None, names=["col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11"])
num_reads_all = df.iloc[0, 0] - df.iloc[1, 0] - df.iloc[2, 0]
num_unmapped_all = df.iloc[0, 0] - df.iloc[4, 0]
num_mapped_all = num_reads_all - num_unmapped_all

# Calc num of reads mapped to taxon contigs
df2 = pd.read_csv(args.tax_flagstat, sep = '\s+', header = None, names=["col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11"])
num_reads_tax = df2.iloc[0, 0] - df2.iloc[1, 0] - df2.iloc[2, 0]
num_unmapped_tax = df2.iloc[0, 0] - df2.iloc[4, 0]
num_mapped_tax = num_reads_tax - num_unmapped_tax

# Calc %
pct = (num_mapped_tax / num_mapped_all) * 100
rnd_pct = round(pct, 2)

# Print results
print(args.tax_flagstat)
print("num_reads" + "\t" + "pct")
print(str(num_mapped_tax) + "\t" + str(rnd_pct) + "%")
print("Number of reads mapped to whole assembly: " + str(num_mapped_all))

