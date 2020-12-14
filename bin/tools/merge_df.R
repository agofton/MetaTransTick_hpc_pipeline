#!/bin/env Rscript

# Usage: merge_df.R blastn.txt blastx.txt TPM.txt output.txt

# Init command line args
args = commandArgs(trailingOnly=TRUE)

# Import blastn & add headers
blastn <- read.delim(args[1], sep="\t", header=FALSE)
colnames(blastn) <- c("SeqID","ssciname_bn","saccver_bn","stitle_bn","staxid_bn","pident_bn","length_bn","mismatch_bn","gapopen_bn","qstart_bn","qend_bn","sstart_bn","send_bn","evalue_bn","bitscore_bn")

# Import blastx & add headers
blastx <- read.delim(args[2], sep="\t", header=FALSE)
colnames(blastx) <- c("SeqID","ssciname_bx","saccver_bx","stitle_bx","staxid_bx","pident_bx","length_bx","mismatch_bx","gapopen_bx","qstart_bx","qend_bx","sstart_bx","send_bx","evalue_bx","bitscore_bx")

# Import TPM.txt & replace headers
TPM <- read.delim(args[3], sep="\t", header=TRUE)
colnames(TPM) <- c("SeqID", "bwa_startpos", "bwa_endpos", "bwa_numreads", "bwa_covbases", "bwa_coverage", "bwa_meandepth", "bwa_meanbaseq", "bwa_meanmapq", "bwa_kb", "bwa_RPK", "bwa_TPM")

# Merge blastn & blastx by SeqID, keep all rows
blastnx <- merge(blastn, blastx, by="SeqID", all=TRUE)

# Merge blastnx & TPM by SeqID, keep only records found in blastnx
blastnxTPM <- merge(blastnx, TPM, by="SeqID", all.x=TRUE)

# Write df
write.table(blastnxTPMm file=args[4], sep="\t", rownames=FALSE)
