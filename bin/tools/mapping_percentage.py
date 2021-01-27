#!/usr/bin/python3

# Alexander W. Gofton, CSIRO, Jan 2021 alexander.gofton@csiro.au; alexander.gofton@gmail.com

import pandas as pd
import os
import argparse

####
# Init argparse
parser = argparse.ArgumentParser(description="Calculates the number and percentage of reads mapped to a bam file, based on flagstat file")
parser.add_argument("-i", "--flagstat", help = "Output flagstat.txt file", type = str, required = True)
parser.parse_args()
args = parser.parse_args()
####

# Import flagstat.txt as pandas df
df = pd.read_csv(args.flagstat, sep = '\s+', header = None, names=["col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11"])

# Get relivant values
total = df.iloc[0, 0]
sec = df.iloc[1, 0]
sup = df.iloc[2, 0]
mapped = df.iloc[4, 0]

# Do math
num_reads = total - sec - sup
num_unmapped = total - mapped
num_mapped = num_reads - num_unmapped 
result = (num_mapped / num_reads) * 100
rounded_result = round(result, 2)

# Print results
print(str(args.flagstat) + ':')
print("num_mapped" + "\t" + "pct_mapped")
print(str(num_mapped) + "\t" + str(rounded_result) + "%")
print()







