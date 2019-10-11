
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
	quit()
	""" GET SEEDS """
	if iterateOverSeedList:
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
			for fileIteration in fileNameList:
				with open(fileIteration) as currentFile:
					lociPositionList = hf.findLociPositions(currentFile)
					for i in lociPositionList:
						listOfPositionsForLoci.append(i)
			listOfPositionsForLoci = hf.unique(listOfPositionsForLoci)
			for fileIteration in fileNameList:
				with open(fileIteration) as currentFile:
					lociList = hf.findLoci(currentFile, listOfPositionsForLoci)
					rawData.append(lociList)
			rawDataCopy = copy.deepcopy(rawData)
			rawDataSorted = []
			for i in rawData:
				sortedListOfLociForEachIndividual = sorted(i, key = itemgetter(2))
				rawDataSorted.append(sortedListOfLociForEachIndividual)
			listOfPositionsForLociSorted = sorted(listOfPositionsForLoci, key = itemgetter(0))
			new_df_list = list()
			for individual in rawDataSorted:
				new_line = list()
				for site in individual:
					new_line.append(site[0:2])
				new_df_list.append(new_line)
			columns_labels = [str(i[0]+" "+i[1]) for i in listOfPositionsForLociSorted]
			df = pd.DataFrame(new_df_list, columns=columns_labels, index=fileNameList)
			print(df)
			df.to_csv(str(organelle+'.csv'), sep="\t")
			""" COLLECT HAPLOTYPES """
			haplotypeList = []
			addedHaplotype1 = True
			addedHaplotype2 = True
			runCycleCount = 0
			while addedHaplotype1 or addedHaplotype2:
				hetMin = 1000000
				runCycleCount = runCycleCount + 1
				""" FIND LEAST-SITE HETETROZYGOTES """
				for key in fileNames:
					hf.heterozygote_count(rawDataSorted[key])
					quit()
				quit()
