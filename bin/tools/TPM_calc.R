#!/usr/bin/env Rscript

# Usage: ./TPM_calc.R samtools_cov_table.txt
# TPM calcs will be appened to the end of samtools_cov_table.txt

library(dplyr)

# input arguments after the scrip.R handle to be parsed as args
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
	  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
	  # default output file
	  args[2] = "out.txt"
}

#Read in the bt2 mapping stats file and select cols
tpm <- read.csv(args[1], header=T, sep="\t") %>% select(X.rname:meanmapq)
# Get length of contig in Kb
tpm$kb <- tpm[,3]/1000
# Divide num_reads / contig length(kb) to give RPK
tpm$RPK <- tpm[,4]/tpm[,10]
# Sum all RPK and divide by 1 mil to give tpm scaling factor
tpm_scaling_factor <- sum(tpm$RPK)/1000000
# Divide RPK by scaling factor
tpm$TPM <- tpm[,11]/tpm_scaling_factor
# Get relivant cols and make new df
tpm <- select(tpm, c(X.rname:meanmapq, TPM))
# write file
write.csv(tpm, file=args[1], quote=FALSE, row.names=FALSE)
