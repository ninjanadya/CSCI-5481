# Nadya Postolaki
# posto018
# CSCI 5481 HW 1

#IMPORTANT INFO:
#Code runs on Python 3
#To run in terminal:
#	$ python3 count_codons.py input.fna output.csv
# where input.fna can be any name file of type .fna
# and output.csv can be any file name of type .csv (the program will make a .csv file with that name)

# References:
# https://pypi.org/project/fastaparser/

import fastaparser
import sys
import csv

#sys.argv[0] == hw1.py

arguments = len(sys.argv) -1
fnaInput = sys.argv[1]
csvOutput = sys.argv[2]
fnaLength  = len(fnaInput)
csvLength  = len(csvOutput)


#If statements below are error checks on the commandline args
if (arguments != 2):
	print("Chief, we've got a problem. You need two arguments, and here you are giving me NOT THAT AMOUNT >:(")
	exit()
elif (fnaLength <= 4):
	print("You gotta enter a .fna file, my dude :(")
	exit()
elif (csvLength <= 4):
	print("Man, you're hurting me. You gotta enter a .csv file, bro :'(")
	exit()
elif ((fnaInput[(fnaLength-4)]+fnaInput[(fnaLength-3)]+fnaInput[(fnaLength-2)]+fnaInput[(fnaLength-1)]) != ".fna"):
	print("E N T E R  A '.fna' F I L E  F I R S T(OuO)")
	exit()
elif ((csvOutput[(csvLength-4)]+csvOutput[(csvLength-3)]+csvOutput[(csvLength-2)]+csvOutput[(csvLength-1)]) != ".csv"):
	print("pls. ppleeeasseeee. just. enter a '.csv' file second. ur making me cry")
	exit()

codonList = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"] 	#this is a surprise tool that will help us later
		#https://i.kym-cdn.com/photos/images/original/001/264/842/220.png
print()
print("Below is what has been saved to a .csv file:")

with open(csvOutput, 'w', newline='') as file: 	#Makes a new .csv file and writes to it
	writer = csv.writer(file)

	with open(sys.argv[1]) as fasta_file:		#.fna file parser
		parser = fastaparser.Reader(fasta_file, parse_method='quick')
		for seq in parser:
			print()
			tSeq = seq.sequence
			length = len(tSeq)-1 		 #removing terminator for the while loop, consequently making my whole life more difficult than it honestly needs to be but I'm just gonna leave it here because my code works and I personally just do not want to go back and change it, so pls just accept :(
			length = length - (length%3) #helps ignore last 1 or 2 extra characters. 
										 #new length is perfectly divisible by 3 now
			i = 0
			while i < (len(codonList)):  #:O remember the list earlier?? we're using it now!
				count = 0
				j = 0
				while j < length:
					if (codonList[i] == tSeq[j]+tSeq[j+1]+tSeq[j+2]): #comparing codons from the codon list to the codons in the genome! :D
						count+=1 
					j+=3 #we skip by 3 because we're comparing strings of 3 characters
				print(codonList[i]+","+str(count)) #making the terminal all pretty with info UwU
				writer.writerow([codonList[i], count])
				i+=1
