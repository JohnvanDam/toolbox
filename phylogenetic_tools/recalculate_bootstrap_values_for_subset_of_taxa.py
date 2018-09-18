'''
    Author: John van Dam
    Created: February 12th, 2018
    Version: 1.0

    Purpose of script: The purpose of this script is to select a subset of taxa
    in a phylogenetic tree and recalculate bootstrap values for those taxa only
    from the bootstrap trees as provided by RAxML, and possibly other programs.
    Currently accepts newick format files only...
'''


import sys
import os
import argparse
import logging
import re
from ete3 import Tree


# Some general stuff
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)-15s %(levelname)-8s \
                    %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)  # Set up the logger
logger.setLevel("INFO")


# Setup the argument parser
parser = argparse.ArgumentParser(description="Prunes the tree and recalculates the bootstrap values.")
parser.add_argument("tree", type=str,
                    help="The phylogenetic tree from which to do the sub taxa \
                    analysis.")
parser.add_argument("bootstraptrees", type=argparse.FileType('r'),
                    help="The file containing all corresponding bootstrap trees.")
parser.add_argument("taxa", type=str,
                    help="Identifying string or regex for requested taxa.")
args = parser.parse_args()


# Compile a regular expression to determine the wanted taxa
try:
    leaf_regex = re.compile(args.taxa)  # Yes, I know. No sanitation of user input at all.
except Exception as e:
    logger.error("Error in regular expression.")
    logger.error(e)
    exit(1)

# Loading and preparing the primary phylogenetic tree
logger.info("Loading and pruning the primary tree.")
try:
    main_tree = Tree(args.tree)
except Exception as e:
    logger.error("Error loading the phylogenetic tree from file {}.".format(args.tree))
    logger.error(e)
    exit(1)
main_tree_leaves_all = set(main_tree.get_leaf_names())

# Get the leaves we want to preserve
selected_main_leave_names = [leaf for leaf in main_tree_leaves_all if leaf_regex.match(leaf)]
logger.info("The following leaf nodes matching your query were found:")
for leaf in selected_main_leave_names:
    logger.info(leaf)

# Prune the main tree
logger.info("Pruning the tree.")
main_tree.prune(selected_main_leave_names, preserve_branch_length=True)  # we want to preserve branch lengths!
logger.info(main_tree.write())

# Loading and preparing the bootstrap trees
logger.info("Loading and pruning the bootstrap trees.")
bootstrap_trees = []  # Empty list
for index, bstree in enumerate(args.bootstraptrees):
    try:
        bs_tree = Tree(bstree.rstrip())
    except Exception as e:
        logger.error("Error loading bootstrap tree no. {} from file {}.".format(index,args.bootstraptrees))
        logger.error(e)
        exit(1)

    # Get the leaves and test against the main tree if they are the same
    if main_tree_leaves_all != set(bs_tree.get_leaf_names()):
        # They are not! Throw the error!
        logger.error("The leaves in the bootstrap tree no. {} are not the same as in the primary tree!".format(index))
        exit(1)

    # Prune the bootstrap tree
    bs_tree.prune(selected_main_leave_names, preserve_branch_length=False)  # we don't need branchlengths
    # And add to the collection of bootstrap trees
    bootstrap_trees.append(bs_tree)
logger.info("Loaded {} bootstrap trees.".format(len(bootstrap_trees)))


# Calculate the new bootstrap scores
# For each node in main_tree, that is not a leaf, count how often you find the same clade in the bs trees
for main_node in main_tree.traverse(strategy="levelorder"):
    if main_node.is_leaf():
        continue
    new_support = 0
    # Get all leaf names from the main tree
    clade_leaf_names = main_tree.get_leaf_names()
    # Now check for each bs_tree if the common ancestor of these same leaves have more leaves
    for bs_tree in bootstrap_trees:

        # Get all node objects for all the leaves by name
        clade_leaf_nodes_in_bs_tree = [bs_tree.get_leaves_by_name(leaf_name)[0] for leaf_name in clade_leaf_names]
        # Get common ancestor in bs_tree
        common_ancestor_node = bs_tree.get_common_ancestor(clade_leaf_nodes_in_bs_tree)
        # Get leafnames of the common ancestor node and check if they are the same
        bs_tree_clade_leaf_names = common_ancestor_node.get_leaf_names()
        if set(clade_leaf_names) == set(bs_tree_clade_leaf_names):
            # Clades match!
            new_support = new_support + 1
    new_support = new_support / len(bootstrap_trees) * 100
    logger.debug("Support for internal node was {}, now is {}".format(main_node.support, new_support))
    main_node.support = new_support

# Output the new tree
logger.info("Writing new bootstrapped tree to STDOUT.")
print(main_tree.write())
