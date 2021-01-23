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
	#print(overallScore)
