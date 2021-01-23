# Nadya Postolaki
# posto018
# CSCI 5481 HW 2

#IMPORTANT INFO:
#Code runs on Python 3
#To run in terminal:
#	$ python3 hw2.py seq1.fna seq2.fna
#where seq1.fna and seq2.fna are any of the fasta files that are included in the zip

#Histograms and csv files are also saved in the zip.
#You'll notice that only the Sequence 1 was scrambled because I recall reading in the Slack conversation that only one needed
#to be scrambled, so I arbitrarily picked sequence 1.

#Referenced the below site for the NW Algorithm and changed some bits to solve for the problem given
#Most of he work comes from reference #3, so you'll see similar, if not the same work for the NW algorithm.
#I am only taking credit for the minor changes that I made in parts of the algorithm and implementing ideas from
#the other links as well, and writing the anchored function and hardcoding the solution.
#1	http://readiab.org/book/0.1.3/2/1#2
#2	https://github.com/mmtechslv/nwunch/blob/master/nwunch.py
#3	https://wilkelab.org/classes/SDS348/2019_spring/labs/lab13-solution.html
#4	https://docs.google.com/document/d/1SJ0PgyzpgQd9tadXMng5-GAoU7_aySL3keaORkLYRKQ/edit


#For the hardcoding portion, I know I could have made functions that did the work for me instead
#of copying and pasting and changing what I needed, but time was running short and I had my hands
#full with some personal problems that came up, so I really apologize for making you look through all of that.
#I'm hoping to use my grace days on this because of what I mentioned a moment ago.
#If there are any major issues, please do reach out to me and I'll be able to discuss them.


import fastaparser
import sys
import csv
import random
import numpy as np #to be able to work with arrays

arguments = len(sys.argv) -1
sec1 = sys.argv[1]
sec2 = sys.argv[2]
seq1Length = len(sec1) 
seq2Length = len(sec2) 

##print(seq1);
##print(seq2);
##print(matches);


#start with error checking to see whether the input files are there and to see if they're the correct type
if (arguments != 2):
	print("Need 2 input args")
	exit()
elif((seq1Length <= 4) or (seq2Length <= 4)):
	print("Please include 2 arguments in the order:\n")
	print("seq1.fna seq2.fna")
elif ((sec1[(seq1Length-4)]+sec1[(seq1Length-3)]+sec1[(seq1Length-2)]+sec1[(seq1Length-1)]) != ".fna"):
	print("First input must be a .fna file")
	exit()
elif ((sec2[(seq2Length-4)]+sec2[(seq2Length-3)]+sec2[(seq2Length-2)]+sec2[(seq2Length-1)]) != ".fna"):
	print("Second input must be a .fna file")
	exit()

#Provided in homework
gap = -2
match = 1
mismatch = -3

#test opening the fasta files and saving the sequences 	
with open(sec1) as fasta_file:
	parser = fastaparser.Reader(fasta_file, parse_method = 'quick')
	for seq in parser:
		print()
		sequence1 = seq.sequence
		#rows = len(sequence1)+1 #initializing number of rows in matrix
	
with open(sec2) as fasta_file:
	parser = fastaparser.Reader(fasta_file, parse_method = 'quick')
	for seq in parser:
		print()
		sequence2 = seq.sequence
		#columns = len(sequence2)+1 #initializing number of columns in matrix

def zeros(rows, cols):
    matrixList = []
    for x in range(rows):
        matrixList.append([])
        for y in range(cols):
            matrixList[-1].append(0)
    return matrixList #returning a matrix of desired size initialized and filled with zeros

#NW Algorithm key:
#Match = 1
#Mismatch = -3
#Gap = -2

def match_score(seq1Char, seq2Char):
    if seq1Char == seq2Char:
        return match
    elif seq1Char == '-' or seq2Char == '-':
        return gap
    else:
        return mismatch

def anchoredNW(seq1, seq2): #hard coding the start and end times. they correspond to different sequences.
	anchoredScore = 0
	anchoredLen = len(seq1) #this is for the while loop to iterate through
	i = 0
	tempAlignment1 = ""
	tempAlignment2 = ""
	while i < anchoredLen:
		if (seq1[i] == seq2[i]): # if i = 0, we begin at the start position for each sequence
			anchoredScore+=match
		else:
			anchoredScore+=mismatch
		tempAlignment1 += seq1[i]
		tempAlignment2 += seq2[i]
		i+=1
	return tempAlignment1, tempAlignment2, anchoredScore

		
	
def needleman_wunsch(seq1, seq2):
    seq1Len = len(seq1)
    seq2Len = len(seq2)  
    matrix = zeros(seq2Len+1, seq1Len+1)
    for i in range(0, seq2Len + 1): #fills out the 0th column with gap scores
        matrix[i][0] = gap * i
    for j in range(0, seq1Len + 1): #fills out the 0th row with gap scores
        matrix[0][j] = gap * j
    for i in range(1, seq2Len + 1): #fills out the rest in the matrix
        for j in range(1, seq1Len + 1): #here we check the actual scoring of each index of the array
            match = matrix[i - 1][j - 1] + match_score(seq1[j-1], seq2[i-1])
            delete = matrix[i - 1][j] + gap
            insert = matrix[i][j - 1] + gap
            matrix[i][j] = max(match, delete, insert) #we want to keep the top score so we use the max function to do that for us
    overallScore = matrix[i][j]
    #print(overallScore)
    
    #initializig alignment strings
    align1 = ""
    align2 = ""
    
    #here we're starting our traceback counters where we left off in our matrix
    i = seq2Len
    j = seq1Len
    
    while i > 0 and j > 0: # end in the topmost or leftmost index
    	score_current = matrix[i][j]
    	score_diagonal = matrix[i-1][j-1]
    	score_up = matrix[i][j-1]
    	score_left = matrix[i-1][j]
    	
    	#Here we're looking to see if there was an insertion, deletion, or not to find where to go next in the matrix
    	#and create the alignment from there
    	if score_current == score_diagonal + match_score(seq1[j-1], seq2[i-1]):
    		align1 += seq1[j-1]
    		align2 += seq2[i-1]
    		i -= 1
    		j -= 1
    	elif score_current == score_up + gap:
    		align1 += seq1[j-1]
    		align2 += '-'
    		j -= 1
    	elif score_current == score_left + gap:
    		align1 += '-'
    		align2 += seq2[i-1]
    		i -= 1
    
    # Finish tracing up to the top left cell
    while j > 0:
    	align1 += seq1[j-1]
    	align2 += '-'
    	j -= 1
    while i > 0:
    	align1 += '-'
    	align2 += seq2[i-1]
    	i -= 1
    	
    #reversing the alignment because it is backwards since i was traversed through the matrix backwards
    align1 = align1[::-1]
    align2 = align2[::-1]
    
    return (align1, align2, overallScore)

finalAlignment1=""
finalAlignment2=""
finalScore=0
tempAlignment1=""
tempAlignment2=""
tempScore = 0

#Match_N
#	1Start	1End	2Start	2End
#	126		187		123		184
#	873		925		870		922
#hardcoding the alignments and scores
if ((sec1 == "SARS-CoV-1_N_protein.fna") and (sec2 == "SARS-CoV-2_N_protein.fna")):
	for i in range(0, 125):
		tempAlignment1+=sequence1[i]
	for j in range(0, 122):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = needleman_wunsch(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	tempAlignment1=""
	tempAlignment2=""
	
	for i in range(126, 187):
		tempAlignment1+=sequence1[i]
	for j in range(123, 184):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = anchoredNW(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	tempAlignment1=""
	tempAlignment2=""
	
	for i in range(186, 872):
		tempAlignment1+=sequence1[i]
	for j in range(185, 869):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = needleman_wunsch(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	tempAlignment1=""
	tempAlignment2=""
	
	for i in range(873, 925):
		tempAlignment1+=sequence1[i]
	for j in range(870, 922):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = anchoredNW(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	tempAlignment1=""
	tempAlignment2=""
	
	for i in range(926, len(sequence1)-1):
		tempAlignment1+=sequence1[i]
	for j in range(923, len(sequence2)-1):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = needleman_wunsch(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	print(finalAlignment1 + "\n\n" + finalAlignment2 + "\n\n" + str(finalScore))
	
	
#Match_S
#	1Start	1End	2Start	2End
#	246		274		255		283
#	937		988		976		1027
#	2985	3106	3039	3160	
elif ((sec1 == "SARS-CoV-1_S_protein.fna") and (sec2 == "SARS-CoV-2_S_protein.fna")):
	for i in range(0, 245):
		tempAlignment1+=sequence1[i]
	for j in range(0, 254):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = needleman_wunsch(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	tempAlignment1=""
	tempAlignment2=""
	
	for i in range(246, 274):
		tempAlignment1+=sequence1[i]
	for j in range(255, 283):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = anchoredNW(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	tempAlignment1=""
	tempAlignment2=""
	
	for i in range(275, 936):
		tempAlignment1+=sequence1[i]
	for j in range(284, 975):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = needleman_wunsch(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	tempAlignment1=""
	tempAlignment2=""
	
	for i in range(937, 988):
		tempAlignment1+=sequence1[i]
	for j in range(976, 1027):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = anchoredNW(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	tempAlignment1=""
	tempAlignment2=""
	
	for i in range(989, 2984):
		tempAlignment1+=sequence1[i]
	for j in range(1028, 3038):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = needleman_wunsch(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	tempAlignment1=""
	tempAlignment2=""
	
	for i in range(2985, 3106):
		tempAlignment1+=sequence1[i]
	for j in range(3039, 3160):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = anchoredNW(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	tempAlignment1=""
	tempAlignment2=""
	
	for i in range(3107, len(sequence1)-1):
		tempAlignment1+=sequence1[i]
	for j in range(3161, len(sequence2)-1):
		tempAlignment2+=sequence2[j]
	tempAlignment1, tempAlignment2, tempScore = needleman_wunsch(tempAlignment1, tempAlignment2)
	finalAlignment1+=tempAlignment1
	finalAlignment2+=tempAlignment2
	finalScore+=tempScore
	print(finalAlignment1 + "\n\n" + finalAlignment2 + "\n\n" + str(finalScore))
	
else:
	print("Error, pick the matching file types")


with open("output.csv", 'w', newline='') as file:
	writer = csv.writer(file)
	origAlign1, origAlign2, origScore = needleman_wunsch(sequence1, sequence2)
	writer.writerow(["Original Sequence", origScore])
	for i in range(1, 101):
		temp=list(sequence1)
		random.shuffle(temp)
		tempSequence=''.join(temp)
		randAlign1, randAlign2, randScore = needleman_wunsch(tempSequence, sequence2)
		writer.writerow(["Random " + str(i), randScore])

#below are other attempts in making the nw Algorithm. Can be ignored.
	
"""
#values provided in homework
match = 1
mismatch = -3
gap = -2
gapChar = "-"
alignmentSeq1 = ""
alignmentSeq2 = ""

def match_score(a, b):
	if a == b:
		return match
	elif a == gapChar or b == gapChar:
		return gap
	else:
		return mismatch

def findingAlignment(start1, end1, start2, end2, matrix):
	i = (end1-start1)+1
	j = (end2-start2)+1
	
	while i>0 and j>0: #start at the bottom right of the matrix
		score_current = matrix[i][j][0]
		score_diagonal = matrix[i-1][j-1][0]
		score_left = matrix[i][j-1][0]
		score_up = matrix[i-1][j][0]
		
		if score_current == score_diagonal + match_score(sequence1[start1+i-1], seqence2[start2+j-1]):
			alignmentSeq1 +=sequence1[start1+i-1]
			alignmentSeq2 +=sequence2[start2+j-1]
			i-=1
			j-=1
		elif score_current == score_left + gap:
			alignmentSeq1 += sequence1[start1+i-1]
			alignmentSeq2 += gapChar
			j-=1
	

def nwAlgorithm(start1, end1, start2, end2):
	rows = (end1-start1)+1
	columns = (end2 - start1)+1
	#initializing the matrix
	matrix = [[[[None] for i in range(2)] for i in range(columns)] for i in range(rows)]
	#matrix[row][column]
	for i in range(rows):
		matrix[i][0] = [gap*i, []] #setting gap scores in first column
	for j in range(columns):
		matrix[0][j] = [gap*j, []] #setting gap scores in first row
	for i in range(1, rows):
		for j in range(1, columns):
			score = match if (sequence1[start1+i-1] == sequence2[start2+j-1]) else mismatch 
			#below we check every possible number one can get by either adding a gap or by adding on the score
			h_val = matrix[i][j-1][0] +gap
			d_val = matrix[i-1][j-1][0] +score
			v_val = matrix[i-1][j][0] +gap
			o_val = [h_val, d_val, v_val]
			#the max value is then the score that we put in the matrix at position [i][j]
			matrix[i][j] = [max(o_val), [i+1 for i,v in enumerate(o_val) if v==max(o_val)]] # h = 1, d = 2, v = 3	
	
	overallScore = matrix[i][j][0]
	score = overallScore
	#print(matrix)
	#print(overallScore)
	#return overallScore
	#print(matrix)
	#print(overallScore)"""
