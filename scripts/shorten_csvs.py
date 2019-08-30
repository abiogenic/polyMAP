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
import numpy as nps

def main(input_dir_current,output_dir_current,min_location,max_location):

	os.chdir(input_dir_current)
	csv_files = glob.glob(os.path.join("*q20.csv"))
	for file in csv_files:

		os.chdir(input_dir_current)
		filename_in = file
		filename_out = str(file.split(".")[0] + "_processed.csv")

		heteroplasmies = pd.read_csv(filename_in, engine="python", sep=",", index_col=1)
		heteroplasmies = heteroplasmies.sort_values(by="Pos")
		heteroplasmies = heteroplasmies.loc[:, ~heteroplasmies.columns.str.match('Unnamed')]
		heteroplasmies = heteroplasmies.drop(columns=['Ea','Ec','Eg','Et','Ed','Ei','Score','GeneProduct'])
		heteroplasmies = heteroplasmies.loc[min_location:max_location]

		os.chdir(output_dir_current)
		print("\t" + filename_out + " processed.")
		heteroplasmies.to_csv(filename_out)