#!/usr/bin/env python

import phylo
import numpy as np
import pandas as pd

datafile = 'NCBI2.csv'
notefile = 'annotation.txt'
go_file = 'humanGO.txt'
phyloTree = phylo.phylo_score(datafile, 'Homo_sapiens', notefile)
go_list = pd.read_csv(go_file, sep='\t', names=['source', 'bioentity_label', 'bioentity_name', 'evidence_type'])
go_gene_set = set(go_list['bioentity_label'])
go_gene_set = list(go_gene_set)
print go_gene_set
phyloTree.wrapper_get_mrs(go_gene_set, 'go_out.csv')

