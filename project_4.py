
#####################################
#           RF Distance             #
#####################################

from ete3 import Tree 

def rfdist(t1, t2):
  result = t1.compare(t2, unrooted=True)
  return result["rf"]

#print(rfdist(test_t1, test_t2))

#####################################
#           Experiment 1            #
#####################################

clustal_rapidnj = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 1/clustal_rapidnj.newick")
kalign_rapidnj = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 1/kalign_rapidnj.newick")
muscle_rapidnj = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 1/muscle_rapidnj.newick")
clustal_quicktree = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 1/clustal_quicktree.newick")
kalign_quicktree = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 1/kalign_quicktree.newick")
muscle_quicktree = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 1/muscle_quicktree.newick")

clustal_rapidnj_permuted = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 2/clustalw_rapidnj_tree_permuted.newick")
kalign_rapidnj_permuted = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 2/kalign_rapidnj_tree_permuted.newick")
muscle_rapidnj_permuted = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 2/muscle_rapidnj_tree_permuted.newick")
clustal_quicktree_permuted = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 2/clustalw_quicktree_permuted.newick")
kalign_quicktree_permuted = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 2/kalign_quicktree_permuted.newick")
muscle_quicktree_permuted = Tree(newick="/Users/caeciliaskov-jensen/Documents/Uni/8. semester/Algorithms in Bioinformatics/15. Projects/Project 4/Experiment 2/muscle_quicktree_permuted.newick")

print(rfdist(muscle_quicktree_permuted, muscle_quicktree))
