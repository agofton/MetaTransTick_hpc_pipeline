#!/usr/bin/python3

import pandas as pd
import argparse

# Initialise user arguments
parser = argparse.ArgumentParser(description = "Calculateds mapping states from samtools flagstat and fasta files.")
parser.add_argument("-t", "--trin_fasta", help = "", type = str, required = True)
parser.add_argument("-l", "--lca", help = "", type = str, required = True)
parser.add_argument("-f", "--tax_fasta", help = "", type = str, required = True)
parser.add_argument("-x", "--tax_flagstat", help = "", type = str, required = True)
parser.add_argument("-y", "--trin_flagstat", help = "", type = str, required = True)
parser.add_argument("-z", "--tax", help = "", type = str, required = True)
parser.add_argument("-o", "--out", help = "", type = str, required = True)
parser.parse_args()
args = parser.parse_args()

# Total number of contigs
ntot = len([1 for line in open(args.trin_fasta) if line.startswith(">")])

# Number of assigned contigs
nass = len(open(args.lca).readlines(  ))

# Number of tax contigs
ntax = len([1 for line in open(args.tax_fasta) if line.startswith(">")])

# Tax contig as pct of assigned contigs
pct_tax = (ntax / nass) * 100
rounded_pct_tax = round(pct_tax, 2)

# nreads mapped to tax
tax_fs = pd.read_csv(args.tax_flagstat, sep = '\s+', header = None, names=["col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11"])
nmapped_tax = (tax_fs.iloc[0, 0] - tax_fs.iloc[1, 0] - tax_fs.iloc[2, 0]) - (tax_fs.iloc[0, 0] - tax_fs.iloc[4, 0])

# nreads mapped to whole assembly
trin_fs = pd.read_csv(args.trin_flagstat, sep = '\s+', header = None, names=["col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9", "col10", "col11"])
nmapped_tot = (trin_fs.iloc[0, 0] - trin_fs.iloc[1, 0] - trin_fs.iloc[2, 0]) - (trin_fs.iloc[0, 0] - trin_fs.iloc[4, 0])

# Calcs.
pct_mapped = (nmapped_tax / nmapped_tot) * 100
rounded_pct_mapped = round(pct_mapped, 2)

pct_ass_contigs = (nass / ntot) * 100
rounded_pct_ass_contigs = round(pct_ass_contigs, 2)

# Print results as df
results_list = [{"Trinity_contigs": ntot, "Assigned_contig": nass, "pct_assigned_contig": rounded_pct_ass_contigs, str(args.tax) + "_contigs": ntax, "pct_" + str(args.tax) + "_contigs": rounded_pct_tax, str(args.tax) + "_mapped_reads": nmapped_tax, "pct_" + str(args.tax) + "_mapped_reads": rounded_pct_mapped}]
results_df = pd.DataFrame(results_list)
print(results_df.to_string(index = False))
results_df.to_csv(args.out, sep = "\t", header = True, index = False)

