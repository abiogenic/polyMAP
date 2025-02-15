#!/usr/bin/env python3

import os
import copy
import glob
import shutil
import random
import time
from operator import *
import hap_functions as hf
import pandas as pd

def main(output_dir_current,repeat,organelle):

	""" GET output_dir_current AND FILES """

	os.chdir(output_dir_current)

	fileNameList = glob.glob(os.path.join("*.vcf"))

	print(fileNameList)

	""" GET SEEDS """

	if repeat:
		seedList = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
	else:
		seedList = [123]

	""" DECLARE MAJOR LISTS """

	rawData=[]
	listOfPositionsForLoci=[]

	""" OPEN THE OUTPUT FILE FOR WRITING """

	with open(str(organelle+'.txt'), 'w') as output:

		""" ITERATE FOR EACH SEED IN SEEDLIST"""

		for seedDraw in seedList:
			random.seed(seedDraw)
			random.shuffle(fileNameList)

			minSitePosition = 0
			sortedListOfLociForEachIndividual = []
			rawDataSorted = []
			rawData = []
			listOfPositionsForLociSorted = []

			""" ITERATE OVER FILES IN output_dir_current """

			# for each file, collect a list of lists of [[refAllele, refFreq, numReads], [varAllele, varFreq, numReads]]

			for fileIteration in fileNameList:
				with open(fileIteration) as currentFile:
					lociPositionList = hf.findLociPositions(currentFile)
					#print(lociPositionList)
					for i in lociPositionList:
						listOfPositionsForLoci.append(i)
			#print(listOfPositionsForLoci)

			# for all files, find the unique positions

			listOfPositionsForLoci = hf.unique(listOfPositionsForLoci)

			#print(listOfPositionsForLoci)

			# for each file, collect a list of lists of [[refAllele, newRefFreq, numReads], [varAllele, newVarFreq, numReads], [lociPosition, refAllele]]

			for fileIteration in fileNameList:
				with open(fileIteration) as currentFile:
					lociList = hf.findLoci(currentFile, listOfPositionsForLoci)
					rawData.append(lociList)

			#quit()

			#print(rawData)
			#print(rawData[0])

			# copy the list just made, this copy can be edited

			rawDataCopy = copy.deepcopy(rawData)

			# re-sort loci so all individuals have same order of loci

			rawDataSorted = []
			for i in rawData:
				sortedListOfLociForEachIndividual = sorted(i, key = itemgetter(-1))
				rawDataSorted.append(sortedListOfLociForEachIndividual)

			#print(rawDataSorted)
			#print(rawDataSorted[0])		#[[['A', 1.0, 1.0], ['A', 0.0, 1.0], ['244953', 'A']], [['A', 1.0, 1.0], ['A', 0.0, 1.0], ['244954', 'A']], [['C', 1.0, 1.0], ['C', 0.0, 1.0], ['244956', 'C']]]

			# as;ldg;alksdg;alskdjf ai dunno what dis is fuck it why do i need so muchsoorting/1/1/1/!?!?!?

			listOfPositionsForLociSorted = sorted(listOfPositionsForLoci, key = itemgetter(0))

			#print(listOfPositionsForLociSorted)			#[['244953', 'A'], ['244954', 'A'], ['244956', 'C']]
			#print(listOfPositionsForLociSorted[0])			#['244953', 'A']

			new_df_list = list()
			for individual in rawDataSorted:
				new_line = list()
				for site in individual:
					new_line.append(site[0:2])
				new_df_list.append(new_line)
			columns_labels = [str(i[0]+" "+i[1]) for i in listOfPositionsForLociSorted]
			df = pd.DataFrame(new_df_list, columns=columns_labels, index=fileNameList)
			print(columns_labels,df)
			df.to_csv(str(organelle+'.csv'), sep="\t")


			""" COLLECT HAPLOTYPES """

			haplotypeList = []

			runCycleCount = 0

			haplotype_not_added = False

			while not haplotype_not_added:
			
				for file_key in range(0,len(fileNameList)):
					het_count = 0
					potential_haplotype = {}

					individual = rawDataSorted[file_key]
					number_of_sites = len(individual)
					for locus in range(0,number_of_sites):
						site = individual[locus]
						alleles = [i[0] for i in site[:-1]]
						
						if len(alleles) == 1:
							pass
						else:
							het_count = het_count + 1
						
						if het_count <= 1:
							potential_haplotype[site[-1][0]] = alleles[0]
						else:
							break
					if het_count <= 1:
						if potential_haplotype not in haplotypeList:
							haplotypeList.append(potential_haplotype)

				# print(haplotypeList)
				# print(len(haplotypeList))

				#for individual in individualList:
				for file_key in range(0,len(fileNameList)):
					individual = rawDataSorted[file_key]
					number_of_sites = len(individual)

				# iterate over each individual
					# check each site, if they are heterozygous, add 1
					# until 1, add that site to potential haplotype
					# at 1, copy the list to add the second haplotype
					# at above 1, drop the potential haplotypes
				# all individuals with 0 have their potential haplotype confirmed in the haplotypeList
				# all individuals with 1 have both of their potential haplotypes confirmed in the haplotypeList

				# iterate over each haplotype in the haplotypeList
					# iterate over each individual
						# iterate over each site
							# if they match the haplotype, set the minimum number they match
						# if a full match was found, iterate over each site again
							# where they match, subtract the minimum number


