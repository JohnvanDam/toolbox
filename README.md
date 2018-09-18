# toolbox
A collection of small scripts for bioinformatics tasks used in publications.

- Phylogenetic tools:
  - recalculate_bootstrap_values_for_subset_of_taxa.py

## Phylogenetic tools
A collection of scripts dealing with phylogenetic data and analyses.

### recalculate_bootstrap_values_for_subset_of_taxa.py
usage: recalculate_bootstrap_values_for_subset_of_taxa.py [-h] tree bootstraptrees taxa

Prunes the tree for specific taxa and recalculates the bootstrap values.

positional arguments:
  tree            The phylogenetic tree on which to do the sub taxa analysis.
  bootstraptrees  The file containing all corresponding bootstrap trees.
  taxa            Identifying string or regex for requested taxa.

optional arguments:
  -h, --help      show this help message and exit
