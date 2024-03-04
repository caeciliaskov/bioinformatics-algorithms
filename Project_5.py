from ete3 import Tree
import time
import numpy as np

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

# A function to write a tree to a file as a Newick-format string 
def write_tree(newick_path, tree, taxon_labels):
    newick_string = make_newick_string(len(tree) - 1, tree, taxon_labels) + ";"
    newick_file = open(newick_path, "w")
    newick_file.write(newick_string)
    newick_file.close()

# Find the coordinates of the lowest value in the matrix
def min_value_coordinate(matrix):
	lowest_value = None
	lowest_value_coordinate = None

	for i, row in enumerate(matrix):
		for j, value in enumerate(row):
			if i != j:
				if lowest_value == None or value < lowest_value:
					lowest_value = value
					lowest_value_coordinate = (i, j)

	return lowest_value_coordinate

# Produce tree from distance matrix
def nj(D, taxon_labels):
    matrix_length = n_taxa = len(D)
    n_nodes = n_taxa + n_taxa - 2
    mat_map = [i for i in range(n_taxa)] # mapping matrix rows and columns to node indices
    r = np.zeros(len(mat_map))
 
    # Compute the neighbor joining N-matrix from a distance matrix
    def compute_N(D):
        N = (1 / (matrix_length - 2.0)) * D 
        for i in range(matrix_length):
            r[i] = np.sum(D[i])
            for j in range(matrix_length):
                r[j] = np.sum(D[j])
                if i != j:
                    N[i][j] -= r[i] + r[j] 
        return N

    # build the tree in newick format   
    T = []   
    for i in range(n_nodes):
        T.append({})

    # Internal nodes
    for k in range(n_taxa, n_nodes): # we call internal nodes k
        if k == n_nodes - 1: 
            i, j = 0, 1 # when this is the seed node, don't have to find the next nodes to branch off
        else: 
            N = compute_N(D)
            i, j = min_value_coordinate(N) # these are the next nodes to branch off
        
        ij_dist = D[i][j]
        i_len = 0.5 * (ij_dist + (1/(len(mat_map)-2) * np.sum(D[i]) - 1/(len(mat_map)-2) * np.sum(D[j]))) 
        j_len = ij_dist - i_len

        # add the edges and branch lengths
        T[k][mat_map[i]] = i_len
        T[k][mat_map[j]] = j_len

        # if this is the seed node, fill in the last root branch length and stop calculating
        if k == n_nodes - 1:
            T[k][mat_map[2]] = D[0][2] - i_len
            break

        new_D = np.zeros((matrix_length-1, matrix_length-1))

        l = 0
        new_mat_map = [k]

        # a and b are the old indices, l and m are the new indices
        for a in range(matrix_length):
            if (a != i) and (a != j): # skip the rows to be merged
                l += 1
                m = 0
                new_mat_map.append(mat_map[a])

                ka_dist = 0.5 * (D[i][a] + D[j][a] - ij_dist)
                new_D[0][l] = ka_dist
                new_D[l][0] = ka_dist

                for b in range(matrix_length): # skip the columns to be merged
                    if (b != i) and (b != j):
                        m += 1
                        new_D[l][m] = D[a][b]
        
        D = new_D
        mat_map = new_mat_map
        matrix_length = matrix_length -1
        
        output_path = "401_our.newick"
    return write_tree(output_path, T, taxon_labels)


################################################

# Test code

################################################

test_mat = np.array([
                    [0.00, 0.23, 0.16, 0.20, 0.17],
                    [0.23, 0.00, 0.23, 0.17, 0.24],
                    [0.16, 0.23, 0.00, 0.20, 0.11],
                    [0.20, 0.17, 0.20, 0.00, 0.21],
                    [0.17, 0.24, 0.11, 0.21, 0.00]
                    ])

test_taxon = ["A", "B", "C", "D", "E"] 

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

dist_mat = get_dist_mat("401_DDE.phy")
tax_lab = get_taxa_labels("401_DDE.phy")

start = time.time()

nj(dist_mat, tax_lab)

end = time.time()

print(end-start)

############################################

# RF-distance

############################################

from ete3 import Tree 

def rfdist(t1, t2):
  result = t1.compare(t2, unrooted=True)
  return result["rf"]

tree1 = Tree(newick = "401_our.newick")
tree2 = Tree(newick = "401_RNJ.newick")

#print(rfdist(tree1, tree2))