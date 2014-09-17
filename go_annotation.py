#!/bin/env python

import pandas as pd
import numpy as np

go_file = 'humanGO.txt'
phylo_file = 'phylo_extended.csv'

if __name__ == '__main__':
    go_list = pd.read_csv(go_file, sep='\t', names=['source', 'bioentity_label', 'bioentity_name', 'evidence_type'])
    phylo_list = pd.read_csv(phylo_file)
    go_list_set = set(go_list['bioentity_label'])
    phylo_list['in_goterm'] = [(_ in go_list_set) for _ in phylo_list['gene']]
    phylo_list.to_csv('phylo_go.csv', index=False) 
    