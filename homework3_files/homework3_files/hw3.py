# Nadya Postolaki
# posto018
# CSCI 5481 HW 3

#IMPORTANT INFO:
#Code runs on Python 3
#To run in terminal:
#	$ python3 hw3.py sequence.fna

#so my code only partially works. 
#I have no idea how many grace days I have left, 
#but if there's any way to get partial credit, I'd really appreciate it.
#My brain is fried and I have too much on my plate with personal stuff going on
#so I'm happy to take whatever I can get from this assignment.
#sorry to make you go through a couple hundred lines of code
#I've literally been blasting coconut mall from mario kart for days on end to say focused
#and i think I'm about to go insane, not even gonna lie lol
#anyways, happy grading??? idk if you actually enjoy this :(

import fastaparser
import sys
import csv
import math
from node import Node

seq = sys.argv[1]

sequenceList = [] #List of all the sequences
headersList = [] #List of the names of each sequence
listPDistances = [] #List of lists :-) building a matrix sort of thing
#Reminder to myself:
#Length of each sequence is 1486! Start at 0 :)

##############################################################
#Gets sequence
##############################################################
with open(seq) as fasta_file:
	parser = fastaparser.Reader(fasta_file, parse_method = 'quick')
	for seq in parser:
		sequence = seq.sequence
		header = seq.header
		sequenceList.append(sequence)
		headersList.append(header)
length = len(sequenceList)


##############################################################
#method to check whether 2 chars match
#wrote this to clean up my work below
##############################################################
def checkDifference(sec1Char, sec2Char): 	
	if (sec1Char == sec2Char):		
		return 0
	else:
		return 1

##############################################################
#Acquires p-distance between each sequence including same sequences
#stores values in a list, which is then stored in another list, thus creating a distance matrix
##############################################################
pDistance = 0
differences = 0	
pDisList = []
tempDisVal = []
pDisValues = []
for i in range(length):
	pDisList = [] 			#resets temp list
	tempDisVal = []
	differences = 0 		#resets number of differences
	pDistance = 0 			#this reset isnt needed but I put it in as precaution
	pDisList.append(headersList[i])	#this is the right axis
	for j in range(length):
		differences = 0
		pDistance = 0
		for k in range (1486):
			differences += checkDifference(sequenceList[i][k], sequenceList[j][k])
		pDistance = (differences/1486)
		tempDisVal.append(pDistance)
		pDisList.append(pDistance)
	listPDistances.append(pDisList) #appends each list to distance matrix
	pDisValues.append(tempDisVal)

##############################################################
#creating our genetic distance file
#each index of the matrix is it's own row in the file
#tab separated
#we begin with the top axis
##############################################################	
geneticDistances = open("genetic-distances.txt", "w") #this is where my p-distances will be stored
geneticDistances.write("\t")
for i in range(len(headersList)):
	geneticDistances.write(str(headersList[i]) + "\t") #writing out the top axis

##############################################################
#The below for loop writes each list of p-distances as a separate row in the txt file
##############################################################
for i in range(len(listPDistances)):
	geneticDistances.write("\n")
	for j in range(len(listPDistances[i])):
		geneticDistances.write(str(listPDistances[i][j]) + "\t")
geneticDistances.close()

##############################################################
#update matrix to make a new smaller one to create that tree
#not really updating, more so making a new matrix
##############################################################
def calcNewMatrix(matrix, mini, minj, length):
	newdistances = [[0] * (length + 1) for _ in range(length + 1)]
	for i in xrange(length):
        	for j in xrange(length):
            		newdistances[i][j] = matrix[i][j]
    
	for k in range(length):
        	newdistances[length][k] = (0.5) * (matrix[mini][k] + matrix[minj][k] - matrix[mini][minj])
        	newdistances[k][length] = newdistances[length][k]

	newmatrix = [[0] * (length - 1) for _ in range(length - 1)]
	newi = newj = 0
##############################################################
#replacing nodes
##############################################################
	for i in range(length + 1):
        	if i == mini or i == minj:
            		continue
        	newj = 0
        	for j in range(length + 1):
            		if j == mini or j == minj:
                		continue
            		newmatrix[newi][newj] = newdistances[i][j]
            		newj += 1
        	newi += 1

	return newmatrix

root = None
#pDisValues is matrix with all the values
def joinNeighbors(count):
	global root #making sure value of root is available outside the function
	relationships = {} #child/parent dictionary
	orderedSeqIds = {}
	for i, ele in enumerate(headersList):
        	orderedSeqIds[ele] = str(i + 1) #creating a new map for the ids

##############################################################
#Gets the next neighbor by taking in a new matrix
#and new set of ids
#changes count (starts at 120)
##############################################################
	def nextN(ids, matrix, count):
        	global root
        	lenMatrix  = len(matrix)
		
##############################################################
#this if statement checks the last nodes and sets it as the root
##############################################################
        	if lenMatrix <= 2:
            		if ids[0] in orderedSeqIds:
                		tempFirstNode = orderedSeqIds[ids[0]]
            		else:
                		tempFirstNode = ids[0]
            		if ids[1] in orderedSeqIds:
                		tempSecNode = orderedSeqIds[ids[1]]
            		else:
                		tempSecNode = ids[1]
           		
            		relationships[tempFirstNode] = (tempSecNode, matrix[0][1])
            		root = tempSecNode
            		return

##############################################################
#else: we know that we aren't at our last nodes, here we make the new q matrix
#could clean this up and make a new function, but Im already late :(
##############################################################
        	qmatrix = [[0]*lenMatrix for _ in range(lenMatrix)]
        	mini = minj = 0 #initialize minimums
        	minVal = float('inf') #pos infinity
        	for i in range(lenMatrix):
        		for j in range(lenMatrix):
        			if i == j:
        				continue
        			qmatrix[i][j] = ((lenMatrix-2) * matrix - sum(matrix[i]) - sum(matrix[j]))
        			if qmatrix[i][j]<minVal:
        				mini = i
        				minj = j
        				minVal = qmatrix[i][j]
        	
        	edgei = (1 / 2.0 * matrix[mini][minj] + 1 / (2.0 * (lenMatrix - 2)) * (sum(matrix[mini]) - sum(matrix[minj])))
        	edgej = (matrix[mini][minj] - edgei)
        	
        	tempFirstNode = tempSecNode = None

##############################################################
#setting up tips
##############################################################
        	if ids[mini] in orderedSeqIds:
            		tempFirstNode = orderedSeqIds[ids[mini]]
        	else:
            		tempFirstNode = ids[mini]
        	if ids[minj] in orderedSeqIds:
            		tempSecNode = orderedSeqIds[ids[minj]]
        	else:
            		tempSecNode = ids[minj]
        
        	relationships[tempFirstNode] = (str(count), edgei) #add to tree
        	relationships[tempSecNode] = (str(count), edgej) #add to tree

        	#reset id
        	ids = [ele for ele in ids if ele not in (ids[mini], ids[minj])] + [str(count)]

        	count -=1
        	newMatrix = calcNewMatrix(matrix, mini, minj, lenMatrix)

        	nextN(ids, newMatrix, count) #do it again bruv

	nextN(headersList, pDisValues, count)

	root_node = Node(root)
	queue = [root_node]
	while len(queue) > 0:
		node_list = []
		for node in queue:
            		for child, (parent, distance) in relationships.items():
                		if parent == node.id:
                    			childNode = Node(child)
                    			node.add_child(childNode, distance)
                    			node_list.append(childNode)
		for node in node_list:
			del relationships[node.id]
		queue = node_list
	return root_node


root = joinNeighbors(120)

##############################################################
#preOrder traversal of a binary tree
##############################################################
def preOrder(root, check):
	if check is None:
		check = []
	if root is None or root.children is None:
		return
	for child, distance in root.children.items():
		check.append((root.id, child.id, distance))
		preOrder(child, check)
	return check

##############################################################
#printing the tree using preorder traversal into the file
##############################################################
check = preOrder(root, None)
edgesFile = open("nj-edges.txt", "w") #this is where my edges will be stored
for (parent, child, distance) in check:
	edgesFile.write(parent + "\t" + child + "\t" + str(distance) + "\n") #printing to file
edgesFile.close()

"""
##############################################################
#uses postOrder traversal of a newick tree to write newick file
##############################################################
def postOrder(node):
	if not node.children:
		if 1 <= int(node.id) <= 61: #there are 61 total sequences, counted earlier
			return headersList
"""
"""
orderedIds = dfs(root, None) #depth first search
origPart = bfs(root) #breadth first search
matches = [0] * 59 #59 nodes

##############################################################
#running 100 times
##############################################################
for _ in range(100):
	bSeq = {}
	index = []
	numIndexes = len
"""











