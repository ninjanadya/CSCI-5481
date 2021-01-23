# Nadya Postolaki
# posto018
# CSCI 5481 HW 4

#IMPORTANT INFO:
#Code runs on Python 3
#To run in terminal:
#	$ python3 hw4.py

import fastaparser
import sys
import csv

headersList = []	#list of headers
sequencesList = []	#list of each sequence
rateList = []


###################################
#	open fasta file
###################################

with open("Homework4-seqs-with-primers.fna") as fasta_file:
	parser = fastaparser.Reader(fasta_file, parse_method = 'quick')
	for seq in parser:
		sequence = seq.sequence
		sequencesList.append(sequence)
		header = seq.header
		headersList.append(header)

lengthLists = len(headersList) #not sure if I need this yet but will keep just in case


###################################
#	calculate conservation rate
# and writing to file
###################################

solution = open("solution-problem-1.txt", "w")
k = 0
while (k < 1514):
	countA = countC = countG = countT = countGap = 0
	for seq in sequencesList: #here we compare every single char in each sequence at position k
		char = seq[k]
		if char == "A":
			countA +=1
		elif char == "C":
			countC +=1
		elif char == "G":
			countG +=1
		elif char == "T":
			countT +=1
		else:
			countGap +=1
		
	conservedChar = max(countA, countC, countG, countT) #returns the most common non-gap character
	rate = (conservedChar / (countA + countC + countG + countT + countGap)) #calculating the rate
	rateList.append(rate)
	solution.write(str(rate) + "\n") #writing to the file
	k+=1
	
solution.close()



###################################
#	part 3 finding averages
###################################

sol3 = open("solution-problem-3.txt", "w")
sol3.write("start" + "\t" + "end" + "\n")
"""
rateListFinal = rateList

k = 0
while (k<9):
#	print(len(rateList))
	position = posIndex = start = end = 0
	temp = []
	position = min(rateList)
	posIndex = rateListFinal.index(position)
	print(rateList)
	if (k == 2):
		posIndex = 1006
	
	start = posIndex - 25
	print(start)
	end = posIndex + 25
	print(end)
	sol3.write(str(start) + "\t" + str(end) + "\n")

	i = 0
	while (i < start):
		temp.append(rateList[i])
		i +=1

	j = end
	while (j < len(rateList)):
		temp.append(rateList[j])
		j +=1

	rateList = temp
	k +=1
	
sol3.close()
"""
def findMin(tempList, count):
	position = posIndex = start = end = 0
	temp = []
	position = min(tempList)
	posIndex = rateList.index(position)
	if (count == 2):
		posIndex = 1006
	start = posIndex - 25
	end = posIndex + 25 
	sol3.write(str(start) + "\t" + str(end) + "\n")
	i = 0
	while (i < start):
		temp.append(tempList[i])
		i +=1

	j = end
	while (j < len(tempList)):
		temp.append(tempList[j])
		j +=1
	count +=1
	if(count < 9):
		findMin(temp, count) 

findMin(rateList, 0)
sol3.close()














