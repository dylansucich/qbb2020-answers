#!/Users/dylansucich/miniconda3/bin/python3

import numpy as np


# HoxD70 matrix of Chiaromonte, Yap, Miller 2002,

def sigma(s,t): # s and t are the bases; want A to be zero, C to be 1, G is 2, T is 3
#use the bases as the index and then return the sigma value as part of the function
    sigma = [ [   91, -114,  -31, -123 ],
              [ -114,  100, -125,  -31 ],
              [  -31, -125,  100, -114 ],
              [ -123,  -31, -114,   91 ] ]    
    index = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
    score = sigma[index[s]][index[t]]
    return score
    
gap = 300
def needlemanwunch(seq1, seq2, gap = -300):
    n = len(seq1) + 1 # number of columns in our matrix, plus header
    m = len(seq2) + 1 # number of rows in our matrix, plus header 

    # initialize alightnment matrix of zeroes
    alignment = np.zeros((m, n), dtype = float) 
    walkback = np.zeros((m, n), dtype= int) # where the pointer will go 

    for i in range(m): 
        alignment[i][0] = (-gap)*i
    for j in range(n): 
        alignment[0][j] = (-gap)*j

    for i in range(1, m):
        for j in range(1, n):
            v = alignment[i - 1][j] + gap # define as 300
            h = alignment[i][j - 1] + gap
            d = alignment[i - 1][j - 1] - sigma(seq1[j-1],seq2[i-1])
            alignment[i][j] = max(v, h, d)

            if alignment[i][j] == d: # diagonal --> match or mismatch
                walkback[i][j] = 1
            elif alignment[i][j] == h: #
                walkback[i][j] = 2
            elif alignment[i][j] == v: #
                walkback[i][j] = 3
	
    return alignment, walkback

def print_alignment(alignment, walkback, seq1, seq2): # Now tracing back
	align1 = ""
	align2 = ""
	i = len(seq2)
	j = len(seq1)

	align_score = alignment[i][j]
	while i > 0 and j > 0:
		if walkback[i][j] == 1: # diagonal --> match or mismatch
			align1 += seq2[i-1]
			align2 += seq1[j-1]
			i -= 1
			j -= 1

		elif walkback[i][j] == 2: # horizontal
			align1 += seq2[i-1]
			align2 += "-"
			i -= 1

		elif walkback[i][j] == 3: # vertical
			align1 += "-"
			align2 += seq1[j-1]
			j -= 1
	
	align1 = align1[::-1]
	align2 = align2[::-1]

	print ("Gap:" + "\n" + str(gap) + "\n")
	print ("Alignment 1:" + "\n" + align1 + "\n") 
	print ("Alignment 2:" + "\n" + align2 + "\n") 
	print ("Score: " + "\n" + str(align_score))

def main():
    seq2 = 'GGGGCTGCCAACACAGAGGTGCAACCATGGTGCTGTCCGCTGCTGACAAGAACAACGTCAAGGGCATCTTCACCAAAATCGCCGGCCATGCTGAGGAGTATGGCGCCGAGACCCTGGAAAGGATGTTCACCACCTACCCCCCAACCAAGACCTACTTCCCCCACTTCGATCTGTCACACGGCTCCGCTCAGATCAAGGGGCACGGCAAGAAGGTAGTGGCTGCCTTGATCGAGGCTGCCAACCACATTGATGACATCGCCGGCACCCTCTCCAAGCTCAGCGACCTCCATGCCCACAAGCTCCGCGTGGACCCTGTCAACTTCAAACTCCTGGGCCAATGCTTCCTGGTGGTGGTGGCCATCCACCACCCTGCTGCCCTGACCCCGGAGGTCCATGCTTCCCTGGACAAGTTCTTGTGCGCCGTGGGCACTGTGCTGACCGCCAAGTACCGTTAAGACGGCACGGTGGCTAGAGCTGGGGCCAACCCATCGCCAGCCCTCCGACAGCGAGCAGCCAAATGAGATGAAATAAAATCTGTTGCATTTGTGCTCCAG'
    seq1 = 'CATAAACCCTGGCGCGCTCGCGGCCCGGCACTCTTCTGGTCCCCACAGACTCAGAGAGAACCCACCATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCCAGGTTAAGGGCCACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCTGCGGTGCACGCCTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAAGCTGGAGCCTCGGTGGCCATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCCCTTCCTGCACCCGTACCCCCGTGGTCTTTGAATAAAGTCTGAGTGGGCGGCAAAAAAAAAAAAAAAAAAAAAA'
    alignment, walkback = needlemanwunch( seq2, seq1 )
    print_alignment( alignment, walkback, seq2, seq1 )

if __name__ == "__main__":
    main()
# note initialize first row and column
# Fill in the reamining rows and columns







