############################

# Project 1
# Cæcilia 201806070
# Lea 201808501
# Peter 201806679

############################

import numpy as np

# Consider the following substitution matrix for DNA sequences:

#    A  C  G  T    
# A 10  2  5  2    
# C  2 10  2  5    
# G  5  2 10  2    
# T  2  5  2 10  

# Question 1: What is the optimal (here maximal) cost of an alignment of AATAAT and AAGG using the above substitution matrix and gap cost -5?

# import timeit

# start = timeit.default_timer()

seq1 = 'AATAAT'
seq2 = 'AAGG'
   
match_scores = {'A': {'A': 10, 'T': 2, 'G': 5, 'C': 2}, 
				'T': {'A': 2, 'T': 10, 'G': 2, 'C': 5}, 
				'G': {'A': 5, 'T': 2, 'G': 10, 'C': 2}, 
				'C': {'A': 2, 'T': 5, 'G': 2, 'C': 10}}

def C(s1,s2, gapscore = -5):
    t_mat = np.zeros(((len(str(s1))+1),len(str(s2))+1))
    for i in range(len(t_mat)):
        t_mat[i][0] = gapscore * i
    for j in range(len(t_mat[0])):
        t_mat[0][j] = gapscore * j
    for i in range(1,len(s1)+1):
        for j in range(1,len(s2)+1):
           v1 = t_mat[i][j-1] + gapscore
           v2 = t_mat[i-1][j] + gapscore
           v3 = t_mat[i-1][j-1] + match_scores[s1[i-1]][s2[j-1]] 
           t_mat[i][j] = max(v1,v2,v3)
    return t_mat

# # print(C(seq1,seq1)) 

# stop = timeit.default_timer()

# print('Time: ', stop - start)

# # The optimal cost is therefore 20. 

# # Question 2: What is the optimal (here maximal) cost of an alignment of seq1.fasta and seq2.fasta using the same substitution matrix and gap cost? (You probably want to implement the algorithm for computing the cost of an optimal alignment.)

def read_fasta_file(filename):
    lines = []
    for l in open(filename).readlines():
        if l[0] != ">" and l[0] != ';':
            lines.append(l.strip())
    return "".join(lines).replace(' ','')

fasta1 = read_fasta_file('seq1.fasta')
fasta2 = read_fasta_file('seq2.fasta')

print(C(fasta1,fasta2)) 

# The optimal cost is therefore 1346. 

# Question 3 (optional): How does an optimal alignment look like for the above two pairs of sequences using the given substitution matrix and gap cost -5? (you probably want to implement the algorithm for finding an optimal alignment by backtracking through the dynamic programming table.)

def optimal_alignment(s1,s2, gapscore = -5):
    optimal_alignment_1 = ''
    optimal_alignment_2 = ''
    current_entry = C(s1,s2, gapscore = -5)[-1][-1]
    i = 1
    j = 1
    while i < len(s1) + 1 and j < len(s2) +1:
        if current_entry ==  C(s1,s2, gapscore = -5)[-i-1][-j] + gapscore:
            optimal_alignment_1 = s1[-i] + optimal_alignment_1 
            optimal_alignment_2 = '-' + optimal_alignment_2
            current_entry = C(s1,s2, gapscore = -5)[-i-1][-j]
            i+=1
        elif current_entry == C(s1,s2, gapscore = -5)[-i][-j-1] + gapscore:
            optimal_alignment_1 = '-' + optimal_alignment_1
            optimal_alignment_2 = s2[-j] + optimal_alignment_2
            current_entry = C(s1,s2, gapscore = -5)[-i][-j-1]
            j+=1
        elif current_entry == C(s1,s2, gapscore = -5)[-i-1][-j-1] + match_scores[s1[-i]][s2[-j]]:
            optimal_alignment_1 = s1[-i] + optimal_alignment_1
            optimal_alignment_2 = s2[-j] + optimal_alignment_2
            current_entry = C(s1,s2, gapscore = -5)[-i-1][-j-1]
            i+=1
            j+=1
            print(i)
    return optimal_alignment_1 + '\n' + optimal_alignment_2 

print(optimal_alignment(seq1, seq2))
print(optimal_alignment(fasta1, fasta2))

# Question 4 (optional): How many optimal alignments are for the above two pairs of sequences using the given substitution matrix and gap cost -5? Explain how you can compute the number of optimal alignments.

#######################################################
from cmath import inf
import numpy
'''
Dependant on packages: numpy & math
Also needs CHR's fasta parser in same dir as this is run from 
That and the input files ofcourse... 
'''


def cost(i:str, j:str)->int:
    ''' 
    Score matrix for match/mismatch
    With translation from Nucleotide to idx
    Takes two nucleotides as input (strings of length 1)
    outputs the score for the given match/mismatch as int
    '''

    score = numpy.array([[10, 2, 5, 2],
                         [2, 10, 2, 5],
                         [5, 2, 10, 2],
                         [2, 5, 2, 10]])

    t = {"A":0,
         "C":1,
         "G":2,
         "T":3}
    return score[t[i],t[j]]

def setup_matrix(seq1:str, seq2:str):
    '''
    Takes two strings as input
    Strings should be two nucleotide sequences we wish to align
    Returns a matrix of size seq1 X seq2

    '''

    matrix = numpy.full([len(seq1)+1,len(seq2)+1], None)
    return matrix

def c(i, j, matrix, gap_cost):
    ''' 
    Input: 
    Two integers i and j, each equivalent to the length of the corresponding nucleotide sequence
    A matrix to work on (numpy array)  as output by setup_matrix function
    A linear gap modifier (int)

    Returns the optimal alignment cost of the two sequenses
    '''

    if matrix[i,j] is None:
        v1, v2, v3, v4 = -inf, -inf, -inf, -inf
        if i > 0 and j > 0:
            v1 = c(i-1,j-1, matrix, gap_cost) + cost(seq1[i-1],seq2[j-1])
        if i > 0 and j >= 0:
            v2 = c(i-1,j, matrix, gap_cost) + gap_cost
        if i >= 0 and j > 0:
            v3 = c(i, j-1, matrix, gap_cost) + gap_cost
        if i == 0 and j == 0:
            v4 = 0
        matrix[i,j] = max(v1, v2, v3, v4)
    return matrix[i,j]

def build_alignment(seq1:str, seq2:str, gap_score:int):
    '''
    Wrapper for the previous functions, call them in order to pipe inputs/outputs
    Takes 3 inputs, 2 sequences, and a gap score
    returns a tuple with the filled out matrix, and the an int the optimal alignment cost
    '''

    matrix = setup_matrix(seq1, seq2)
    i, j = len(seq1), len(seq2)
    out = c(i, j, matrix, gap_score)
    return  (matrix, out)

