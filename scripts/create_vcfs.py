#!/usr/bin/env python

import re
from collections import Counter
import random as rand
import string
import time
import os
import glob
import sys
import pandas as pd
import numpy as np

def main(input_dir_current,output_dir_current):
	
	os.chdir(output_dir_current)
	csv_files = glob.glob(os.path.join("*_processed.csv"))
	for file in csv_files:

		filename_in = file
		filename_out = file.split('.')[0]+".vcf"
		genome_type = "chrMT"
		heteroplasmies = pd.read_csv(filename_in, engine="python", sep=",", index_col=0)

		#=======OUTPUT=======#

		with open(filename_out, 'w') as output:

			output.write('##fileformat=VCFv4.0\n')
			output.write('##reference=chrRSRS\n')
			output.write('##FORMAT=<ID=GT,Number=.,Type=String,Description="Genotype">\n')
			output.write('##FORMAT=<ID=DP,Number=.,Type=Integer,Description="Reads covering the REF position">\n')
			output.write('##FORMAT=<ID=HF,Number=.,Type=Float,Description="Heteroplasmy Frequency of variant allele">\n')
			output.write('##FORMAT=<ID=CILOW,Number=.,Type=Float,Description="Value defining the lower limit of the confidence interval of the heteroplasmy fraction">\n')
			output.write('##FORMAT=<ID=CIUP,Number=.,Type=Float,Description="Value defining the upper limit of the confidence interval of the heteroplasmy fraction">\n')
			output.write('##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes">\n')
			output.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">)\n')
			output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTESTSAMPLENAME"+"\n")
			print(heteroplasmies)
			for row in heteroplasmies.index:
				position = str(row)
				array = heteroplasmies.loc[row]
				ref = str(array["Ref"]).capitalize()
				bases = ["A","C","G","T"]
				alt_bases = list()
				bases.remove(ref)
				total = str(array["Total"])
				total_reads = int()

				alts = str()
				alt_reads = str()

				for base in bases:
					if str(array[base]) != '0':
						alt_bases.append(base)

				for base in alt_bases:
					alts = alts + str(base) + ","
					alt_reads = alt_reads + str(array[base]) + ","
					total_reads = total_reads + int(str(array[base]))
				alts = alts.strip(',')
				alt_reads = alt_reads.strip(',')

				hf = round(int(total_reads)/int(total),4)

				if str(alts) != '':
					new_line = str(genome_type+"\t"+position+"\t.\t"+ref+"\t"+alts+"\t.\tPASS\tAC=1;AN=1\tGT:DP:HF:CILOW:CIUP:TYPE\t"+alt_reads+":"+total+":"+str(hf)+":0.0:1.0:misml\n")
					output.write(new_line)
