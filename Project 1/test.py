import numpy as np

# def prepare_matrix(length_seq1, length_seq2, gap_score):
# 	matrix = []
# 	for i in range(length_seq1):
# 		row = []
# 		for j in range(length_seq2):
# 			row.append(None)
# 		matrix.append(row)
# 	for i in range(len(matrix)): 
# 		matrix[i][0] = gap_score * i 
# 	for j in range(len(matrix[0])): 
# 		matrix[0][j] = gap_score * j 
# 	return matrix 

def prepare_matrix(length_seq1, length_seq2):
	matrix = np.zeros((length_seq1+1,length_seq2+1))
	for i in range(len(matrix)): 
		matrix[i][0] = -5 * i 
	for j in range(len(matrix[0])): 
		matrix[0][j] = -5 * j 
	return matrix 

# print(prepare_matrix(6, 4))

def fill_matrix(s1, s2, match_scores, gap_score):
	matrix = np.zeros((len(s1)+1,len(s2)+1))
	for i in range(len(matrix)): 
		matrix[i][0] = gap_score * i 
	for j in range(len(matrix[0])): 
		matrix[0][j] = gap_score * j
	for i in range(1,len(s1)+1):
		for j in range(1, len(s2)+1):
			base1 = s1[i-1]
			base2 = s2[j-1]
			cell_left = matrix[i][j-1] + gap_score
			cell_above = matrix[i-1][j] + gap_score
			cell_diagonal = matrix[i-1][j-1] + match_scores[base1][base2]
			matrix[i][j] = max(cell_left, cell_above, cell_diagonal)
	return matrix

def get_traceback_arrow(matrix, row, col, match_scores, gap_score):
	score_diagonal = matrix[row-1][col-1]
	score_left = matrix[row][col-1]
	score_up = matrix[row-1][col]
	score_current = matrix[row][col]
	if score_current == score_diagonal + match_scores:
		return "diagonal"
	elif score_current == score_left + gap_score:
		return "left"
	elif score_current == score_up + gap_score:
		return "up"

def trace_back(seq1, seq2, matrix, score_matrix, gap_score):
	aligned1 = ""
	aligned2 = ""
	row = len(seq1)
	col = len(seq2)
	while row > 0 and col > 0:
		base1 = seq1[row-1]
		base2 = seq2[col-1]
		match_score = score_matrix[base1][base2]
		traceback_arrow = get_traceback_arrow(matrix, row, col, match_score, gap_score)
		if traceback_arrow == "diagonal":
			aligned1 = base1 + aligned1
			aligned2 = base2 + aligned2
			row -= 1
			col -= 1
		elif traceback_arrow == "up":
			aligned1 = base1 + aligned1
			aligned2 = "-" + aligned2
			row -= 1
		elif traceback_arrow == "left":
			aligned1 = "-" + aligned1
			aligned2 = base2 + aligned2
			col -= 1
	while row > 0:
		base1 = seq1[row-1]
		aligned1 = base1 + aligned1
		aligned2 = "-" + aligned2
		row -=1
	while col > 0:
		base2 = seq2[col-1]
		aligned1 = "-" + aligned1
		aligned2 = base2 + aligned2
		col -= 1
	return [aligned1, aligned2]

def align(seq1, seq2, match_score, gap_score):
	matrix = fill_matrix(seq1, seq2, match_score, gap_score)
	return trace_back(seq1, seq2, matrix, match_score, gap_score)

#############################################################
# Code for calling and testing your functions should be below
# here. If you separate function definitions from the rest of
# your script in this way, you are less likely to make mistakes.
#############################################################

s1 = "TCCAGAGA"
s2 = "TCGAT"

match_scores = {'A': {'A': 10, 'T': 2, 'G': 5, 'C': 2}, 
				'T': {'A': 2, 'T': 10, 'G': 2, 'C': 5}, 
				'G': {'A': 5, 'T': 2, 'G': 10, 'C': 2}, 
				'C': {'A': 2, 'T': 5, 'G': 2, 'C': 10}}


# print(fill_matrix(s1, s2, match_scores, -5))

# Reading in fasta file:
def read_fasta(filename):
    seqs = list()
    info = list()
    seq = ''
    inf = ''
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq != '':
                    seqs.append(seq)
                    info.append(inf)
                    seq = ''
                inf = line[1:]
            else:
                seq += line
        seqs.append(seq)
        info.append(inf)
    return "".join(seqs).replace(" ", "")

seq1 = read_fasta("seq1.fasta")
seq2 = read_fasta("seq2.fasta")

# optimal = fill_matrix(seq1, seq2, match_scores, -5)
# print(optimal[-1][-1])

alignment = align(seq1, seq2, match_scores, -5)
for s in alignment:
    print(s)
