Nadya Postolaki

2)
-Code runs on Python 3
-To run in terminal:
-	$ python3 count_codons.py input.fna output.csv
- where input.fna can be any name file of type .fna
- and output.csv can be any file name of type .csv (the program will make a .csv file with that name)

3)
-The fake genome file named "fakeGenome.fna" is expected to have exactly 1 of each codon, as described in the header, with an extra "T" to test whether it will count the extra character. Output is found in "fakeGenomeOutput.csv" and shows the correct amount of each codon, as expected.

4)
-The actual covid genome files were run and saved to "wholeGenomeOutput.csv" and "separateGenomeOutput.csv"

5)
The figure is called "SequencesVsWholeGenome.png"

6)
The figure is called "AminoAcidSepVsWhole.png"

7)
The largest discrepency is definitely Valine with a difference of 604 counts where the separated genome has 1152 count and the whole genome only 548 count. This is due to the shifting of the frame which allows us to find more distinguishable codons and easily identify the amino acids vs a random frame shift in the genome as a whole.
