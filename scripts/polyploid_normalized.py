#!/usr/bin/env python

import os
import copy
import glob
import shutil
import random
from operator import *

##### 
def findLociPositions(file, f, showStatus):
	if showStatus:
		print("FINDING SITE LOCATIONS:")
	#Fetch file name
	fileName = file
	#Create list of all loci for this file
	lociList = []
	lociPositionList = []
	#For each line in current .vcf file
	content = f.readlines()
	for line in content:
		#See if line starts with #
		ifComment = line.startswith("#")
		#Print the line if not a comment			
		if ifComment == 0:			
			# print(line)
			loci = line.split('\t')[1]
			#print(loci)
			#to do: work with multiple variations of the sample
			freqData = line.split('\t')[9]
			allVarFreq = freqData.split(':')[2]

			allVarAllele = line.split('\t')[4]
			varAllele = allVarAllele.split(',')[0]

			varFreq = float(allVarFreq.split(',')[0])
			refFreq = float(round(1-varFreq, 3))

			refAllele = line.split('\t')[3]
			numReads = float(freqData.split(':')[1])
			lociPosition = [line.split('\t')[1], refAllele]
			newLine = [[refAllele, refFreq, numReads], [varAllele, varFreq, numReads]]
			lociList.append(newLine)
			lociPositionList.append(lociPosition)
	#Print when finished with current file
	if showStatus:
		print("FINISHED WITH " + file)
		print(lociPositionList)
	return(lociPositionList)

def findLoci(fileIteration, currentFile, lociPositionList, showStatus):
	#Fetch file name
	fileName = fileIteration
	#Create list of all loci for this file
	lociList = []
	#lociPositionList = []
	#For each line in current .vcf file
	content = currentFile.readlines()
	for line in content:
		#See if line starts with #
		ifComment = line.startswith("#")
		#Evaluate the line if not a comment			
		if not ifComment:			
			#print(line)
			loci = line.split('\t')[1]
			#print(loci)
			#to do: work with multiple variations of the sample
			freqData = line.split('\t')[9]
			refAllele = line.split('\t')[3]
			numReads = float(freqData.split(':')[1])
			lociPosition = [line.split('\t')[1], refAllele]

			allVarAllele = line.split('\t')[4]
			varAllele = allVarAllele.split(',')[0]
			
			allVarFreq = freqData.split(':')[2]

			sumVarFreq = 0
			print(allVarFreq)
			for freq in allVarFreq.split(','):
				sumVarFreq = sumVarFreq + float(freq)

			varFreq = float(allVarFreq.split(',')[0])
			refFreq = float(round(1-sumVarFreq, 3))

			newVarFreq = varFreq / varFreq + refFreq
			newRefFreq = 1 - newVarFreq

			newLine = [[refAllele, newRefFreq, numReads], [varAllele, newVarFreq, numReads], lociPosition]
			lociList.append(newLine)
			#lociPositionList.append(lociPosition)
	#Print when finished with current file
	#print("Finished with " + file)
	#print(type(lociList))
	newLociList = copy.deepcopy(lociList)

	if showStatus:
		print("CREATING LOCI LIST FOR: " + str(fileIteration))
	#print(lociPositionList)
	if showStatus:
		print("FOUND INPUT:")	
	for i in lociPositionList:
		#print(i)
		#print(lociList)
		if i not in [locus[2] for locus in lociList]:
			if showStatus:
				print(str(i) + " NOT FOUND")
			blankList = [[i[1],1.0,1.0],[i[1],0.0,1.0],i]
			newLociList.append(blankList)
		else:
			if showStatus:
				print("IN")
	if showStatus:
		print("CREATED INPUT:")
		for i in newLociList:
			print(i)
	return(newLociList)

#####
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

def main(output_dir_current,runOnce,run,showStatus,extendedOutput,iterateOverSeedList):
	##### READ ME
	#Put in directories with .vcf files to compare

	# wd = os.path.abspath(__file__)
	# wd = wd.rsplit("/", 1)[0]+"/"
	# #print(wd)
	# wdBack = wd.rsplit("/", 2)[0]+"/"
	# #print(wdBack)

	##### EXAMPLE FORMAT
	#chrMT	460	.	G	T	.	PASS	AC=1;AN=1	GT:DP:HF:CILOW:CIUP	1:11:1.0:0.741:1.0

	listOfListsForLociForEachIndividual=[]
	listOfPositionsForLoci=[]

	##### PARSE FILES
	#Run through each .vcf file in wd (working directory)
	#Each file is "currentFile"
	#lociList is dict of each site and variant
	#{'4694': [('C', 0.0, '10'), ('A', 1.0, '10')], '6571': ... etc
	#masterLociList is list of lociLists

	if iterateOverSeedList:
		seedList = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
	else:
		seedList = [123]

	#seedDraw = int(random.sample(range(len(seedList)), 1)[0])

	#random.seed(int(seedList[seedDraw]))

	os.chdir(output_dir_current)

	fileNameList = glob.glob(os.path.join("*_processed.vcf"))

	with open('test.txt', 'w') as output:

		# for seedDraw in seedList:

		# 	output.write("SEED: " + str(seedDraw) + "\n")
		# 	random.seed(seedDraw)
		# 	random.shuffle(fileNameList)
		# 	output.write(str(fileNameList) + "\n")

		for seedDraw in seedList:

			random.seed(seedDraw)
			print("SEED: " + str(seedDraw))

			random.shuffle(fileNameList)
			print(fileNameList)

			minSitePosition = 0
			sortedListOfLociForEachIndividual = []
			listOfListsForLociForEachIndividualSorted = []
			listOfListsForLociForEachIndividual = []
			listOfPositionsForLociSorted = []

			output.write("\nSEED: " + str(seedDraw) + "\n")

			print("FILES FOUND:")
			print(fileNameList)

			for fileIteration in fileNameList:
				#For current file (currentFile)
				with open(fileIteration) as currentFile:
					lociPositionList = findLociPositions(fileIteration, currentFile, showStatus)
					for i in lociPositionList:
						listOfPositionsForLoci.append(i)

			listOfPositionsForLoci = unique(listOfPositionsForLoci)

			output.write("\nPOSITIONS FOUND: \n")
			output.write(str(listOfPositionsForLoci) + "\n")

			if showStatus:

				print("UNIQUE POSITIONS FOUND:")
				print(listOfPositionsForLoci)

			for fileIteration in fileNameList:
				#For current file (currentFile)
				with open(fileIteration) as currentFile:
					lociList = findLoci(fileIteration, currentFile, listOfPositionsForLoci, showStatus)
					#print(lociList)
					listOfListsForLociForEachIndividual.append(lociList)

			# output.write("\nINDIVIDUALS FOUND: \n")
			# for individual in listOfListsForLociForEachIndividual:
			# 	output.write(str(individual) + "\n")

			# for i in range(0,len(listOfListsForLociForEachIndividual)):
			# 	#print(i)
			# 	if len(listOfListsForLociForEachIndividual[i]) == 0:
			# 		print(fileNameList[i])
			# 		print("0")
			# 	else:
			# 		print(fileNameList[i])
			# 		print(listOfListsForLociForEachIndividual[i])
			# 		#for j in i:
			# 		#	print(j)

			# print(listOfPositionsForLoci)

			if run:
				if showStatus:
					print("LIST OF EACH GENOTYPE FOR THESE LOCI:")
					for i in range(0,len(listOfPositionsForLoci)):
						print(i)
						#print(listOfPositionsForLoci[i])
						print(listOfListsForLociForEachIndividual)
						for listOfLociForEachIndividual in listOfListsForLociForEachIndividual:
							print(listOfLociForEachIndividual[i])

					print("EACH INDIVIDUAL:")
					for i in range(0,len(listOfListsForLociForEachIndividual)):
						print(fileNameList[i])
						for j in listOfListsForLociForEachIndividual[i]:
							if showStatus:
								print(j)

			listOfListsForLociForEachIndividualCopy = copy.deepcopy(listOfListsForLociForEachIndividual)

			listOfListsForLociForEachIndividualSorted = []
			for i in listOfListsForLociForEachIndividual:
				sortedListOfLociForEachIndividual = sorted(i, key = itemgetter(2))
				listOfListsForLociForEachIndividualSorted.append(sortedListOfLociForEachIndividual)

			listOfPositionsForLociSorted = sorted(listOfPositionsForLoci, key = itemgetter(0))

			# output.write("\nSORTED INDIVIDUALS FOUND: \n")
			# for individual in listOfListsForLociForEachIndividualSorted:
			# 	output.write(str(individual) + "\n")

			if showStatus:
				print(listOfPositionsForLociSorted)

			haplotypeList = []

			addedHaplotype1 = True
			addedHaplotype2 = True

			runCycleCount = 0

			while addedHaplotype1 or addedHaplotype2:
				hetMin = 1000000
				runCycleCount = runCycleCount + 1
				if showStatus:
					print("Run Cycle: " + str(runCycleCount))
				addedHaplotype1 = False
				addedHaplotype2 = False

				if not run:
					break 

				if showStatus:
					print("FINDING LEAST-SITE HETETROZYGOTES:")
				hetCountList = []
				for i in range(0,len(listOfListsForLociForEachIndividualSorted)):
					if showStatus:
						print(fileNameList[i])
					hetCountForSite = -1
					for j in listOfListsForLociForEachIndividualSorted[i]:
						if showStatus:
							print(j[0],j[1])
							print(j[0][0], j[1][0])
						# check for frequencies of 0, if present site is homozygote
						if j[0][1] > 0 and j[1][1] > 0:
								if j[0][0] != j[1][0]:
									if hetCountForSite < 0:
										hetCountForSite = 0
									hetCountForSite = hetCountForSite + 1
						if xor(j[0][1] == 0, j[1][1] == 0):
							if hetCountForSite < 0:
								hetCountForSite = 0
					if hetCountForSite >= 0:
						if hetCountForSite < hetMin:
							hetMin = hetCountForSite
						if hetCountForSite == hetMin:
							minSite = listOfListsForLociForEachIndividualSorted[i]
							minSitePosition = i
						if showStatus:
							print("het count = " + str(hetCountForSite))
						hetCountList.append(hetCountForSite)

				if showStatus:
					print("hetMin: ")
					print(hetMin)
				#print(listOfListsForLociForEachIndividualSorted)
				#print("----------")
				#print(listOfListsForLociForEachIndividualSorted)

				for j in minSite:
					if showStatus:
						print(j)

				if showStatus:
					print("FINDING UNAMBIGIOUS HAPLOTYPES:")
				
				if hetMin > 1:
					print("ERROR: no remaining single-site (or no-site) HETETROZYGOTES")

				elif hetMin == 0 or hetMin == 1:
					haplotype1 = []
					haplotype2 = []
					if showStatus:
						print("Min Site:")
						print(minSite)

					#thisLocus = thisIndividual[j]
					# minFraction = 1
					# 	for thisAllele in thisLocus[0:1]:
					# 		if showStatus:
					# 			print("thisAllele: " + str(thisAllele))
					# 		if haplotype[j] == thisAllele[0] and thisAllele[1] > 0:
					# 		 	doesMatch = True
					# 		 	if minFraction > thisAllele[1]:
					# 		 		minFraction = thisAllele[1]
					# 		 	break

					if showStatus:
						print(minSite)

					#listOfPossibleHaplotypes = []
					
					for j in minSite:
						print(j[0])
						if j[0][1] > 0:
							haplotype1.append(j[0][0])
							if not runOnce:
								addedHaplotype1 = True
							#j[0][1] = 0.0
						elif j[1][1] > 0:
							haplotype1.append(j[1][0])
							if not runOnce:
								addedHaplotype1 = True
							#j[1][1] = 0.0


						if j[1][1] > 0:
							haplotype2.append(j[1][0])
							if not runOnce:
								addedHaplotype2 = True
							#j[1][1] = 0.0
						elif j[0][1] > 0:
							haplotype2.append(j[0][0])
							if not runOnce:
								addedHaplotype2 = True
							#j[1][1] = 0.0

					if addedHaplotype1:
						if len(haplotype1) == len(listOfPositionsForLociSorted):
							if haplotype1 not in haplotypeList:
								haplotypeList.append(haplotype1)
								if showStatus:
									print("ADDING: " + str(haplotype1))
						else:
							if showStatus:
								print(haplotype1)
							print("ERROR: haplotype1 for " + str(fileNameList[minSitePosition]) + " is not as long as expected")

					if addedHaplotype2:
						if len(haplotype2) == len(listOfPositionsForLociSorted):
							if haplotype2 not in haplotypeList:
								haplotypeList.append(haplotype2)
								if showStatus:
									print("ADDING: " + str(haplotype2))
						else:
							if showStatus:
								print(haplotype2)
							print(listOfListsForLociForEachIndividualSorted)
							print("ERROR: haplotype2 for " + str(fileNameList[minSitePosition]) + " is not as long as expected")

					#del listOfListsForLociForEachIndividualSorted[minSitePosition]

				if showStatus:
					print(haplotypeList)

				if showStatus:
					print("FIND EXISTING HAPLOTYPES IN OTHER SAMPLES:")

				for haplotype in haplotypeList:
					if showStatus:
						print("---------------------------new haplotype")
						print(haplotype)
					for i in range(0,len(listOfListsForLociForEachIndividualSorted)):
						thisIndividual = listOfListsForLociForEachIndividualSorted[i]
						#print(thisIndividual)
						if showStatus:
							print("-------------new individual:" + str(thisIndividual))
						
						resolvable = True
						if showStatus:
							print(range(0,len(thisIndividual)))
						minFraction = 1
						for j in range(0,len(thisIndividual)):
							if showStatus:
								print("allele position: " + str(j))
							thisLocus = thisIndividual[j]
							#if thisIndividual[2] == [['C', 0.25, '4'], ['T', 0.75, '4'], ['4694', 'C']]:
								#pdb.set_trace()
							#print("---new locus")
							doesMatch = False
							for thisAllele in thisLocus[0:2]:
								if showStatus:
									print("thisAllele: " + str(thisAllele))
									print("haplotype: " + str(haplotype))
								if haplotype[j] == thisAllele[0] and thisAllele[1] > 0:
								 	doesMatch = True
								 	if thisAllele[1] < minFraction:
								 		minFraction = thisAllele[1]
								 		print(minFraction)
								 	break
							if not doesMatch:
								resolvable = False
								break
						if resolvable:
							#print("TRUE")
							#print("minFraction: " + str(minFraction))
							for j in range(0,len(thisIndividual)):
								thisLocus = thisIndividual[j]
								if showStatus:
									print("---new locus")
								for thisAllele in thisLocus:
									if haplotype[j] == thisAllele[0] and thisAllele[1] > 0:
										#print("HERE")
										if showStatus:
											print(thisAllele[0])
											print(thisAllele[1])
										thisAllele[1] = thisAllele[1] - minFraction
										if showStatus:
											print("newAllele: " + str(thisAllele))
										break

						# if haplotype == ['C', 'T', 'T']:
						# 	if listOfListsForLociForEachIndividualSorted[-1] != [[['A', 0.5, '2'], ['C', 0.5, '2'], ['12545', 'A']], [['G', 0.0, '2'], ['T', 1.0, '2'], ['460', 'G']], [['C', 0.5, '2'], ['T', 0.5, '2'], ['4694', 'C']]]:
						# 		pdb.set_trace()

						resolvable = True
						if showStatus:
							print(range(0,len(thisIndividual)))
						for j in range(0,len(thisIndividual)):
							if showStatus:
								print(j)
							thisLocus = thisIndividual[j]
							if showStatus:
								print("---new locus")
							doesMatch = False
							for thisAllele in thisLocus[0:2]:
								if showStatus:
									print("thisAllele: " + str(thisAllele))
								if haplotype[j] == thisAllele[0] and thisAllele[1] > 0:
								 	doesMatch = True
								 	break
							if not doesMatch:
								resolvable = False
								break

						if resolvable:
							if showStatus:
								print("TRUE")
								print("minFraction: " + str(minFraction))
							for j in range(0,len(thisIndividual)):
								thisLocus = thisIndividual[j]
								if showStatus:
									print("---new locus")
								for thisAllele in thisLocus:
									if haplotype[j] == thisAllele[0] and thisAllele[1] > 0:
										if showStatus:
											print("HERE")
											print(thisAllele[0])
											print(thisAllele[1])
										thisAllele[1] = thisAllele[1] - minFraction
										if showStatus:
											print("newAllele: " + str(thisAllele))
										break

			#print("HAPLOTYPE LIST:")
			#print(haplotypeList)

			#output.write(listOfPositionsForLociSorted)

			if run:

				output.write('RESOLVED TOTAL OF ' + str(len(haplotypeList)) + ' HAPLOTYPES: \n')
				output.write('[')
				lengthOfHaplotype = len(haplotype)
				counter = 0
				for position in listOfPositionsForLociSorted:
					output.write("'" + str(position[0]) + "'")
					counter = counter + 1
					if counter != lengthOfHaplotype:
						output.write(',')
						
				output.write(']')
				output.write('\n')

				for haplotype in haplotypeList:
					output.write('[')
					counter = 0
					for i in haplotype:
						lengthOfHaplotype = len(haplotype)
						for j in i:
							if type(j) == float:
								j = round(j,3)
							output.write("'" + str(j) + "'")
							counter = counter + 1
							if counter != lengthOfHaplotype:
								output.write(',')
					output.write(']\n')

			if extendedOutput:		
				output.write('\n')
				output.write('ADJUSTED INPUT: \n')
				for fileIteration in range(0,len(fileNameList)):
					output.write("\n")
					output.write(fileNameList[fileIteration])
					for listOfLociForThisIndividual in listOfListsForLociForEachIndividualSorted[fileIteration]:
						output.write('\n\t')
						for i in listOfLociForThisIndividual:
							output.write('[')
							lengthOfLoci = len(i)
							counter = 0
							for j in i:
								if type(j) == float:
									j = round(j,3)
								output.write("'" + str(j) + "'")
								counter = counter + 1
								if counter != lengthOfLoci:
									output.write(',')
							output.write(']\t')

					output.write('\n')

				output.write('\n')
				output.write('INPUT FOUND: \n')
				for fileIteration in range(0,len(fileNameList)):
					output.write("\n")
					output.write(fileNameList[fileIteration])
					for listOfLociForThisIndividual in listOfListsForLociForEachIndividualCopy[fileIteration]:
						output.write('\n\t')
						for i in listOfLociForThisIndividual:
							output.write('[')
							lengthOfLoci = len(i)
							counter = 0
							for j in i:
								if type(j) == float:
									j = round(j,3)
								output.write("'" + str(j) + "'")
								counter = counter + 1
								if counter != lengthOfLoci:
									output.write(',')
							output.write(']\t')

					output.write('\n')
				output.write('\n')


	#print(len(fileNameList))
	#print(fileNameList)
	#print(len(hetCountList))
	#print(hetCountList)
