#!/usr/bin/env python3

import os
import glob
import shutil
import copy
import math

lgm_to_round_to = 0.05

def round_nearest(x, a):
    return round(round(x / a) * a, -int(math.floor(math.log10(a))))


def findLociPositions(f):
	lociList = []
	lociPositionList = []
	content = f.readlines()
	for line in content:
		ifComment = line.startswith("#")		
		if not ifComment:
			#print(line)		#chrMT	244956	.	C	A	.	PASS	AC=1;AN=1	GT:DP:HF:CILOW:CIUP:TYPE	2:66:0.0303:0.0:1.0:misml
			currentLine = line.split('\t')
			lociPosition = [currentLine[1], currentLine[3]]
			lociPositionList.append(lociPosition)

	return(lociPositionList)

def findLoci(currentFile, lociPositionList):
	lociList = []
	content = currentFile.readlines()
	for line in content:
		ifComment = line.startswith("#")	
		if not ifComment:
			currentLine = line.split('\t')
			lociPosition, refAllele, allVarAllele, freqData = currentLine[1], currentLine[3], currentLine[4], currentLine[9]
			
			varReads, allReads, allVarFreq, lowCI, highCI, varType = freqData.split(':')

			#to do: work with multiple variations of the sample
			numReads = float(allReads)
			varAllele = allVarAllele.split(',')[0]

			varFreq = round(float(allVarFreq.split(',')[0]),3)
			refFreq = round(float(1-varFreq), 3)

			sumFreq = varFreq + refFreq

			newVarFreq = float(round(varFreq / sumFreq, 3))
			newVarFreq = float(round(newVarFreq, 3))
			newRefFreq = 1 - newVarFreq

			newVarFreq = round_nearest(newVarFreq, lgm_to_round_to)
			#print(newVarFreq)
			newRefFreq = 1 - newVarFreq

			newLine = [refAllele, newRefFreq, numReads], [varAllele, newVarFreq, numReads], [lociPosition, refAllele]
			lociList.append(newLine)
			
	newLociList = copy.deepcopy(lociList)

	# IF NO VARIANTS AT THIS SITE, HOMOZYGOTE WITH RATIO 1:0 IS CREATED
	for i in lociPositionList:
		if i not in [locus[-1] for locus in newLociList]:
			blankList = [i[1],1.0,1.0],i
			newLociList.append(blankList)
	return(newLociList)

def unique(list1):
 
    # intilize a null list
    unique_list = []
     
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    # return list
    return(unique_list)
    # throw error if references do not match

def collectLoci(fileNameList):

	loci = []

	# iterate over each file in the list and collect unique loci

	for fileIteration in fileNameList:
		with open(fileIteration) as currentFile:
			lociPositionList = findLociPositions(currentFile)
			for i in lociPositionList:
				loci.append(i)

	# for all files, find the unique positions

	loci = unique(loci)
	return(loci)

def heterozygote_count(individual):
	number_of_sites = len(individual)
	for locus in range(0,number_of_sites):
		site = individual[locus]
		number_of_alleles = len(site)-1
		for allele in range(0,number_of_alleles):
			print(site[allele])
