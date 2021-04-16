#!/usr/bin/python3

# Alexander W. Gofton, CSIRO, Jan 2021 alexander.gofton@csiro.au; alexander.gofton@gmail.com

import pandas as pd
import argparse


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

# Import blastn lca and filter unassigned transcripts
blastn_lca = pd.read_csv(args.lca_bn, header = None, sep = "\t", names = ["transcript_id", "blastn_tax_path"])
unass = ["NCBI;", "NCBI;cellular organisms;", "NCBI;No hits;", "NCBI;Not assigned;"]
blastn_lca_unass = blastn_lca[blastn_lca["blastn_tax_path"].isin(unass) == False]

# Import blastx lca and filter unassigned transcripts
blastx_lca = pd.read_csv(args.lca_bx, header = None, sep = "\t", names = ["transcript_id", "blastx_tax_path"])
blastx_lca_unass = blastx_lca[blastx_lca["blastx_tax_path"].isin(unass) == False]

# Import blastn topHits
blastn_topHits = pd.read_csv(args.bn, sep = "\t", header = None, dtype = str,
	names = ["transcript_id", "blastn_sciname", "blastn_accession", "blastn_title", "blastn_tax_id", 
	"blastn_perc_id", "blastn_length", "blastn_missmatch", "blastn_gapopen", "blastn_qstart", 
	"blastn_qend", "blastn_sstart", "blastn_send", "blastn_evalue", "blastn_bitscore"])

# Import blastx topHits
blastx_topHits = pd.read_csv(args.bx, sep = "\t", header = None, dtype = str, 
	names = ["transcript_id", "blastx_sciname", "blastx_accession", "blastx_title","blastx_tax_id", 
	"blastx_perc_id", "blastx_length", "blastx_missmatch", "blastx_gapopen", "blastx_qstart", 
	"blastx_qend", "blastx_sstart", "blastx_send", "blastx_evalue", "blastx_bitscore"])

# Import tpm
tpm = pd.read_csv(args.tpm, sep = "\t")
tpm.rename(columns={"#rname" : "transcript_id"}, inplace = True)
tpm.rename(columns={"endpos" : "transcript_length"}, inplace = True)

# merge lca files, keeping all records
lca_merged = pd.merge(blastn_lca_unass, blastx_lca_unass, on = "transcript_id", how = "outer")
# Add tpm to lca_merged, keeping only records in lca_merged
lca_tpm = pd.merge(lca_merged, tpm, on = "transcript_id", sort = False, how = "left")
# Add in blastn topHits
lca_tpm_blastn = pd.merge(lca_tpm, blastn_topHits, on = "transcript_id", how = "left")
# Add in blastx topHits
lca_tpm_blastn_blastx = pd.merge(lca_tpm_blastn, blastx_topHits, on = "transcript_id", how = "left")


# Filtering unwanted cols
cols_to_keep=["transcript_id", "transcript_length", "TPM", "blastn_tax_path", "blastx_tax_path", 
	"blastn_sciname", "blastn_title", "blastn_accession", "blastn_perc_id", "blastn_length", 
	"blastn_evalue", "blastn_bitscore","blastx_sciname", "blastx_title", "blastx_accession", 
	"blastx_perc_id", "blastx_length", "blastx_evalue", "blastx_bitscore"]

final_df = lca_tpm_blastn_blastx[cols_to_keep]
final_df.to_csv(args.out, index = False, sep = "\t", na_rep = 'N/A')











