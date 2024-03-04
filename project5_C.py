import numpy as np
import time

# Recursively build a Newick-format string from an adjacency list
def make_newick_string(node_i, tree, taxon_labels):
	if len(tree[node_i]) == 0: # no outgoing edges, so must be a leaf node
		return taxon_labels[node_i]
	else: # an internal node
		newick_substrings = []
		for child_i in sorted(tree[node_i]):
			branch_length = tree[node_i][child_i]
			substring = make_newick_string(child_i, tree, taxon_labels)
			newick_substrings.append("%s:%f" % (substring, branch_length))

		return "(" + ",".join(newick_substrings) + ")"

# A function to write a tree as a Newick-format string 
def write_tree(newick_path, tree, taxon_labels):
    newick_string = make_newick_string(len(tree) - 1, tree, taxon_labels) + ";"
    newick_file = open(newick_path, "w")
    newick_file.write(newick_string)
    newick_file.close()

# Find the coordinates of whatever off-diagonal element of a matrix has the lowest value
def get_lowest_off_diagonal_value_coordinate(matrix):
	lowest_value = None
	lowest_value_coordinate = None

	for i, row in enumerate(matrix):
		for j, value in enumerate(row):
			if i != j:
				if lowest_value == None or value < lowest_value:
					lowest_value = value
					lowest_value_coordinate = (i, j)

	return lowest_value_coordinate

# Compute the neighbor joining Q-matrix from a distance matrix
def compute_q_matrix(distance_matrix):
    matrix_length = len(distance_matrix)
    q_matrix = (matrix_length - 2.0) * distance_matrix
    for i in range(matrix_length):
        for j in range(matrix_length):
            if i != j:
                q_matrix[i][j] -= np.sum(distance_matrix[i]) + np.sum(distance_matrix[j])
    return q_matrix

# Produce tree from distance matrix
def nj(distance_matrix, taxon_labels):
    matrix_length = n_taxa = len(distance_matrix)
    n_nodes = n_taxa + n_taxa - 2
    matrix_map = [n for n in range(n_taxa)] # mapping matrix rows and columns to node indices
    # build the tree in newick format
    tree = []   
    for i in range(n_nodes):
        tree.append({})

    for u in range(n_taxa, n_nodes): # we call internal nodes "u"
        if u == n_nodes - 1:
            f, g = 0, 1 # when this is the seed node, don't have to find the next nodes to branch off
        else:
            q_matrix = compute_q_matrix(distance_matrix)
            f, g = get_lowest_off_diagonal_value_coordinate(q_matrix) # these are the next nodes to branch off

        fg_distance = distance_matrix[f][g]
        f_length = 0.5 * fg_distance + (np.sum(distance_matrix[f]) - np.sum(distance_matrix[g])) / (2.0 * matrix_length - 2)
        g_length = fg_distance - f_length

        # add the edges and branch lengths
        tree[u][matrix_map[f]] = f_length
        tree[u][matrix_map[g]] = g_length

        # if this is the seed node, fill in the last root branch length and stop calculating
        if u == n_nodes - 1:
            tree[u][matrix_map[2]] = distance_matrix[0][2] - f_length
            break
        new_distance_matrix = np.zeros((matrix_length - 1, matrix_length - 1))

        # a and b are the old indices, i and j are the new indices
        i = 0
        new_matrix_map = [u]
        for a in range(matrix_length):
            if (a != f) and (a != g): # skip the rows to be merged
                i += 1
                j = 0

                new_matrix_map.append(matrix_map[a])

                ua_distance = 0.5 * (distance_matrix[f][a] + distance_matrix[g][a] - fg_distance)
                new_distance_matrix[0][i] = ua_distance
                new_distance_matrix[i][0] = ua_distance

                for b in range(matrix_length): # skip the columns to be merged
                    if (b != f) and (b != g):
                        j += 1
                        new_distance_matrix[i][j] = distance_matrix[a][b]

        distance_matrix = new_distance_matrix
        matrix_map = new_matrix_map
        matrix_length = matrix_length - 1

        output_path = "test_C.newick"
    return write_tree(output_path, tree, taxon_labels)


############################################

# Test code

############################################

test_mat = np.array([
                    [0.00, 0.23, 0.16, 0.20, 0.17],
                    [0.23, 0.00, 0.23, 0.17, 0.24],
                    [0.16, 0.23, 0.00, 0.20, 0.11],
                    [0.20, 0.17, 0.20, 0.00, 0.21],
                    [0.17, 0.24, 0.11, 0.21, 0.00]
                    ])


def get_dist_mat(filepath):
    phylip_file = open(filepath, 'r')
    l = []
    l = [line.split() for line in phylip_file]
    dist_mat = []
    for i in range(1,len(l)):
        dist_mat.append(l[i][1:])
    dist_mat = np.array(dist_mat, dtype=float)
    return dist_mat

def get_taxa_labels(filepath):
    phylip_file = open(filepath, 'r')
    l = []
    l = [line.split() for line in phylip_file]
    taxon_labels = []
    for j in range(1,len(l)):
        taxon_labels.append(l[j][0])
    taxon_labels = np.array(taxon_labels)
    return taxon_labels



dist_mat = get_dist_mat("89_Adeno_E3_CR1.phy")
tax_lab = get_taxa_labels("89_Adeno_E3_CR1.phy")

test_tax = ["A", "B", "C", "D", "E"]

# start = time.time()

# nj(test_mat, test_tax)

# end = time.time()

# print("The time of execution of above program is :", round(end-start, 3))

############################################

# RF-distance

############################################

