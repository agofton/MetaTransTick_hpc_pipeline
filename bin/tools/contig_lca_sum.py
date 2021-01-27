#!/usr/bin/python3

# Alexander W. Gofton, CSIRO, Jan 2021 alexander.gofton@csiro.au; alexander.gofton@gmail.com

import os
import sys
import pandas as pd
import argparse
from termcolor import colored, cprint

# Init argparse
parser = argparse.ArgumentParser(description="Combines lca taxonomy .csv, blastn top hits .txt, blastx (diamond) top hits .txt, and contigs statistics .txt into once .csv summary file.")

# Add arguments
parser.add_argument("-b", "--lca_bn", help = "blastn_lca.csv (complete path)", type = str, required = True)
parser.add_argument("-x", "--lca_bx", help = "blastx_lca.csv (complete path)", type = str, required = True)
parser.add_argument("-n", "--bn", help = "blastn_topHits.txt (complete path)", type = str, required = True)
parser.add_argument("-d", "--bx", help = "blastx_topHits.txt (complete path)", type = str, required = True)
parser.add_argument("-t", "--tpm", help = "contig_stats.txt (complete path)", type = str, required = True)
parser.add_argument("-o", "--out", help = "output.txt (complete path)", type = str, required = True)

# Finish args init
parser.parse_args()
args = parser.parse_args()

################################
# Import dfs

# Import blastn_lca.csv and give col headings
lca_bn = pd.read_csv(args.lca_bn, header = None)
lca_bn.columns = ["contig_id", "blastn_tax_path"]
cprint("Imported blastn lca file: " + str(args.lca_bn), "blue", attrs = ["bold"])

# Import blastx(diamond) lca.csv and give col headings
lca_bx = pd.read_csv(args.lca_bx, header = None)
lca_bx.columns = ["contig_id", "blastx_tax_path"]
cprint("Imported diamond lca file: " + str(args.lca_bx), "blue", attrs = ["bold"])

# Import blastn results and give col names
bn = pd.read_csv(args.bn, sep = "\t", header = None, dtype = str)
bn.columns = ["contig_id", "bn_sciname", "bn_accession", "bn_title", "bn_tax_id", "bn_perc_id", "bn_length", "bn_missmatch", "bn_gapopen", "bn_qstart", "bn_qend", "bn_sstart", "bn_send", "bn_evalue", "bn_bitscore"]
cprint("Imported blastn top hits file: " + str(args.bn), "blue", attrs = ["bold"])

# Import blastx (diamond) results and give col names
bx = pd.read_csv(args.bx, sep = "\t", header = None)
bx.columns = ["contig_id", "bx_sciname", "bx_accession", "bx_title", "bx_tax_id", "bx_perc_id", "bx_length", "bx_missmatch", "bx_gapopen", "bx_qstart", "bx_qend", "bx_sstart", "bx_send", "bx_evalue", "bx_bitscore"]
cprint("Imported diamond top hits file: " + str(args.bx), "blue", attrs = ["bold"])

# Import TPM.txt, tab separated, and rename 1st col header
tpm = pd.read_csv(args.tpm, sep = "\t")
tpm.rename(columns={"#rname" : "contig_id"}, inplace = True)
tpm.rename(columns={"endpos" : "contig_length"}, inplace = True)
cprint("Imported contig read mapping stats file: " + str(args.tpm), "blue", attrs = ["bold"])

# Filter out unassigned reads from lca dfs
unass_bn = ["NCBI;", "NCBI;cellular organisms;", "NCBI;No hits;", "NCBI;Not assigned;"]
lca_bn_ass = lca_bn[lca_bn["blastn_tax_path"].isin(unass_bn) == False]

unass_bx = ["NCBI;", "NCBI;cellular organisms;", "NCBI;Not assigned;"]
lca_bx_ass = lca_bx[lca_bx["blastx_tax_path"].isin(unass_bx) == False]

# All dfs now have "contig_id" col in common - will use that to merge dfs
# 1. Merge lca dfs, keeping contig_id keys from boths dfs to give all contigs that have an assigned taxon string
# 2. Add tpm data to lca_merged, using contig_id keys from lca_merged df only
# 3. Merge blastn and blastx results, keeping contig_id keys from both dfs
# 4. Merge blast_merged with lca_tpm, keeping contig_id keys from lca_tpm only - only contig_ids that have an assigned taxon string

lca_merged = pd.merge(lca_bn_ass, lca_bx_ass, on = "contig_id", sort = False, how = "outer")
lca_tpm = pd.merge(lca_merged, tpm, on = "contig_id", sort = False, how = "left")
blast_merged = pd.merge(bn, bx, on = "contig_id", sort = False, how = "outer")
lca_tpm_blast = pd.merge(lca_tpm, blast_merged, sort = False, how = "left")

# Filtering relivant cols 
cols_to_keep = ["contig_id", "blastn_tax_path", 
                "blastx_tax_path", "contig_length", 
                "numreads", "coverage", "meandepth", 
                "TPM", "bn_sciname", "bn_accession", 
                "bn_title", "bn_tax_id", "bn_perc_id", 
                "bn_length", "bn_evalue", "bn_bitscore", 
                "bx_sciname", "bx_accession", "bx_title", 
                "bx_tax_id", "bx_perc_id", "bx_length", 
                "bx_evalue", "bx_bitscore"]

final_df = lca_tpm_blast[cols_to_keep]

# Write file as tap separated .txt
final_df.to_csv(args.out, index = False, sep = "\t")
cprint("Writing contig summary file " + str(args.out) + " ...Complete.", "green", attrs = ["bold", "blink"])
