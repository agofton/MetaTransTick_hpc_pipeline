#!/usr/bin/env Rscript

# arg1 = input.txt (samtools coverage output file)
# arg2 = output.txt (can be same as input)

# inpute arguments after the scrip.R handle to be parsed as args
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
	  stop("At least one argument must be supplied (input file).n", call.=FALSE)

} else if (length(args)==1) {
	  # default output file
	  args[2] = "out.txt"

}

# R steps to calc TMP from samtools coverage

# import samtools coverage data frame
TPM_calc <- read.delim(args[1], header=T)

# add column which given contig length in kb
TPM_calc$kb <- TPM_calc[,3]/1000

# add RPK column
TPM_calc$RPK <- TPM_calc[,4]/TPM_calc[,10]

# calc total RPK/1mil and store as varibale TotRPK
totRPK <- sum(TPM_calc$RPK)/1000000

# add final TMP col
TPM_calc$TPM <- TPM_calc[,11]/totRPK

# save as .txt
write.table(TPM_calc,args[1],sep="\t",row.names=FALSE)
