#!/usr/bin/python3

# Alexander W. Gofton, 2021, alexander.gofton@gmail.com; alexander.gofton@csiro.au

import pandas as pd
import datetime
import argparse
from termcolor import colored, cprint

# Call current date and time
def date_and_time():
    now = datetime.datetime.now()
    print(now.strftime("%d-%m-%Y %H:%M:%S"))

# Initialise user arguments
parser = argparse.ArgumentParser(description = "Takes output from samtools (v1.10) coverage and does RPK and TPM calculations for each contig, adding the result to the last two columns of the input file. Output will overwrite input unless --output [-o] is specified.")
parser.add_argument("-i", "--input", help = "samtools coverage output.txt", type = str, required = True)
parser.add_argument("-o", "--output", help = "output.txt", type = str, required = False)
parser.parse_args()
args = parser.parse_args()

# Check is args.output argument exists, if not, args.output = args.input

cprint("Input is: " + str(args.input), "cyan", attrs = ["bold"])
cprint("Output is: " + str(args.output), "cyan", attrs = ['bold'])

df = pd.read_csv(args.input, sep = "\t")

df["RPK"] = df["numreads"] / (df["endpos"] / 1000)
total_RPK = df["RPK"].sum()
df["TPM"] = df["RPK"] / total_RPK

df.to_csv(args.output, sep = "\t", index = False, header = True)

cprint("TPM calculations completed ", "green", attrs = ["bold"])
date_and_time()