
import numpy as np

def global_linear(s1, s2, gap_score = 5):
    # Define substitution matrix:
    def cost(i:str, j:str)->int:
        score = np.array([[10, 2, 5, 2],
                            [2, 10, 2, 5],
                            [5, 2, 10, 2],
                            [2, 5, 2, 10]])
        t = {"A":0,
            "C":1,
            "G":2,
            "T":3}
        return score[t[i],t[j]]

    # Fill out score matrix:
    def C(s1,s2, gapscore = 5):
        t_mat = np.zeros(((len(str(s1))+1),len(str(s2))+1))
        for i in range(len(t_mat)):
            t_mat[i][0] = gapscore * i
        for j in range(len(t_mat[0])):
            t_mat[0][j] = gapscore * j
        for i in range(1,len(s1)+1):
            for j in range(1,len(s2)+1):
                v1 = t_mat[i][j-1] + gapscore
                v2 = t_mat[i-1][j] + gapscore
                v3 = t_mat[i-1][j-1] + cost(s1[i-1],s2[j-1])
                t_mat[i][j] = min(v1,v2,v3)
        return t_mat

    # Make alignment:
    def get_alignment(s1, s2, gap_score):
        align1, align2 = '', ''
        score = C(s1, s2, gapscore = 5)
        i,j = len(s1),len(s2) # start from the bottom right cell
        while i > 0 and j > 0: # end toching the top or the left edge
            score_current = score[i][j]
            score_diagonal = score[i-1][j-1]
            score_up = score[i][j-1]
            score_left = score[i-1][j]
            if score_current == score_diagonal + cost(s1[i-1], s2[j-1]):
                align1 = s1[i-1] + align1
                align2 = s2[j-1] + align2
                i -= 1
                j -= 1
            elif score_current == score_left + gap_score:
                align1 = s1[i-1] + align1
                align2 = '-' + align2
                i -= 1
            elif score_current == score_up + gap_score:
                align1 = '-' + align1
                align2 = s2[j-1] + align2
                j -= 1
        # Finish tracing up to the top left cell
        while i > 0:
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        while j > 0:
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            j -= 1

        return align1 + "\n" + align2 + "\n" + 'Optimal cost: {}'.format(score[-1][-1])
    return get_alignment(s1, s2, gap_score)

def global_affine(s1, s2):
    # Define score matrix:
    def cost(i:str, j:str)->int:
        score = np.array([[0, 5, 2, 5],
                            [5, 0, 5, 2],     
                            [2, 5, 0, 5],
                            [5, 2, 5, 0]])
        t = {"A":0,
            "C":1,
            "G":2,
            "T":3}
        return score[t[i],t[j]]

    # Fill out score matrix:
    def affine_align(s1, s2, gap_start = 5, gap_extend = 5):
        S = np.zeros(((len(str(s1))+1),len(str(s2))+1))
        I = np.zeros(((len(str(s1))+1),len(str(s2))+1))
        D = np.zeros(((len(str(s1))+1),len(str(s2))+1))
        for i in range(1, len(s1)+1):
            S[i][0] = float('inf')
            I[i][0] = float('inf')
            D[i][0] = gap_start + (i * gap_extend)
        for i in range(1, len(s2)+1):
            S[0][i] = float('inf')
            I[0][i] = gap_start + (i * gap_extend)
            D[0][i] = float('inf')
        for i in range(1, len(s1)+1):
            for j in range(1, len(s2)+1):
                S[i][j] = cost(s1[i-1], s2[j-1]) + min(
                    S[i-1][j-1], 
                    I[i-1][j-1], 
                    D[i-1][j-1])
                I[i][j] = min(
                    gap_start + gap_extend + S[i-1][j], 
                    gap_extend + I[i-1][j], 
                    gap_start + gap_extend + D[i-1][j])
                D[i][j] = min(
                    gap_start + gap_extend + S[i][j-1], 
                    gap_start + gap_extend + I[i][j-1], 
                    gap_extend + D[i][j-1])
        return S

    # Make alignment:
    def backtrace(s1, s2, gap_start = 5, gap_extend = 5):
        S = affine_align(s1, s2, gap_start=5, gap_extend=5)
        align1, align2 = [], []
        i, j = len(s1), len(s2)
        while i > 0 or j > 0:
            if i > 0 and j > 0 and S[i][j] == S[i-1][j-1] + cost(s1[i-1], s2[j-1]):
                align1.append(s1[i-1]) 
                align2.append(s2[j-1])
                i -= 1
                j -= 1
            else:
                k = 1
                while True:
                    if i >= k and S[i][j] == S[i-k][j] + (gap_start + gap_extend * k):
                        align1.append(s1[i-k:i])
                        align2.append('-' * k)
                        i = i-k
                        break
                    elif j >= k and S[i][j] == S[i][j-k] + (gap_start + gap_extend * k):
                        align1.append('-' * k)
                        align2.append(s2[j-k:j])
                        j = j-k
                        break
                    else:
                        k = k + 1
        return "".join(reversed(align1)) + "\n" + "".join(reversed(align2))
    return backtrace(s1, s2)    

#######################################

# Test code

#######################################

t_seq1 = "tatggagagaataaaagaactgagagatctaatgtcgcagtcccgcactcgcgagatactcactaagaccactgtggaccatatggccataatcaaaaag".upper()
t_seq2 = "atggatgtcaatccgactctacttttcctaaaaattccagcgcaaaatgccataagcaccacattcccttatactggagatcctccatacagccatggaa".upper()

test_seq1 = "TCCAGAGA"
test_seq2 = "TCGAT"

seq1 = 'AATAAT'
seq2 = 'AAGG'


print(global_linear(seq1,seq2, 5)) 
print(global_affine(seq1, seq2))
