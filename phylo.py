#!/usr/bin/env python
import os
import os.path, time
from os import listdir 
from os.path import join
import pandas as pd
import numpy as np
import math
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import scipy.stats as stats
import fastcluster

class phylo_score:
    def __init__(self, datafile, host_organ, notefile):
        self.datafile = datafile
        self.blast_data = pd.read_csv(datafile, index_col=0)
        self.host_organ = host_organ
        self.annotation = gene_note(notefile)
        self.preprocess()
        self.hclust()
        
    def preprocess(self, homology_cutoff=50, conserve_cutoff=5):
        self.blast_data[self.blast_data < homology_cutoff] = 1
        conserved = (self.blast_data>=50).sum(1)
        self.normalized = self.blast_data[conserved > conserve_cutoff]
        self.normalized = self.normalized.div(self.normalized[self.host_organ], axis=0)
        self.normalized = np.log(self.normalized)/np.log(2)
        self.normalized = self.normalized.drop(self.host_organ, 1)
        self.normalized = self.normalized.add(-self.normalized.mean(axis=0), axis=1)
        self.normalized = self.normalized.div(self.normalized.std(axis=0), axis=1)


    def hclust(self):
        link_file = self.datafile + '.link.npy'
        if os.path.isfile(link_file) and os.path.getmtime(link_file) >= os.path.getmtime(self.datafile):
            self.link_matrix = np.load(link_file)
        else:
            blast_score = self.normalized.as_matrix()
            self.link_matrix = fastcluster.linkage(blast_score, method='average', 
                                                   metric='correlation', 
                                                   preserve_input=False)
            del blast_score
            np.save(link_file, self.link_matrix)
            
        self.gene_num = self.normalized.shape[0]
        self.node_num = self.gene_num + self.link_matrix.shape[0]
        self.parent_tree = np.array(np.arange(self.node_num))
        self.leaf_num = np.array([1] * self.gene_num + 
                                 [0] * (self.node_num - self.gene_num))
        for i in range(self.link_matrix.shape[0]):
            assert(self.parent_tree[self.link_matrix[i, 0]] == int(self.link_matrix[i, 0]))
            assert(self.parent_tree[self.link_matrix[i, 1]] == int(self.link_matrix[i, 1]))
            assert(self.leaf_num[self.gene_num + i] == 0)
            self.parent_tree[self.link_matrix[i, 0]] = self.gene_num + i
            self.parent_tree[self.link_matrix[i, 1]] = self.gene_num + i
            self.leaf_num[i + self.gene_num] = self.leaf_num[self.link_matrix[i, 0]] + \
                                            self.leaf_num[self.link_matrix[i, 1]]
            
    def cluster_analysis(self, genes, outfile, cluster_cutoff=2):
        gene_idx = self.genes_to_idx(genes)
        dist_matrix_file = self.datafile + '.npdist.npy'
        if os.path.isfile(dist_matrix_file) and os.path.getmtime(dist_matrix_file) >= os.path.getmtime(self.datafile):
            dist_matrix = np.load(dist_matrix_file)
        else:
            if not hasattr(self, 'corr_matrix'):
                corr_file = self.datafile + '.corr.npy'
                if os.path.isfile(corr_file) and os.path.getmtime(corr_file) >= os.path.getmtime(self.datafile):
                    self.corr_matrix = np.load(corr_file)
                else:
                    self.corr_matrix = np.corrcoef(self.normalized.as_matrix())
                    np.save(corr_file, self.corr_matrix)
            top_item_num = 50
            temp = self.corr_matrix.argsort(axis=1)
            rank = temp.argsort(axis=1)
            dist_matrix = top_item_num + rank - self.gene_num + 1 
            dist_matrix[dist_matrix < 0 ] = 0
            temp = np.transpose(dist_matrix) * dist_matrix
            self.dist_matrix = np.sqrt(temp)
            np.save(self.datafile + '.npdist.npy', self.dist_matrix)
        print "distance calculated..."
        linkage_corr = fastcluster.linkage(dist_matrix, 'weighted', 'euclidean')
        clusters = fcluster(linkage_corr, cluster_cutoff, criterion='distance')
        np.save(self.datafile + '.npcluster.npy', clusters) 
        print "clusters generated..."
        pvalues = np.array([1] * len(clusters))
        goi_clusters = set(clusters[gene_idx])
        significant_clusters = []
        for goi_cluster in goi_clusters:
            cluster_size = sum(clusters == goi_cluster)
            intersect_size = sum(clusters[gene_idx] == goi_cluster)
            pvalue = stats.hypergeom.sf(intersect_size, intersect_size + self.gene_num - cluster_size,
                                        len(gene_idx), cluster_size)
            pvalues[np.where(clusters == goi_cluster)] = pvalue
            
        self.cluster_result = self.normalized.iloc[:, :1]
        self.cluster_result.iloc[:, 0] = clusters
        idx_in_cluster = [ i for i in range(len(clusters)) if pvalues[i] < 0.05 ]
        self.cluster_result = self.cluster_result.iloc[idx_in_cluster, :]
        self.cluster_result['gene'] = self.annotation.get_gene(self.cluster_result.index)
        self.cluster_result['description'] = self.annotation.get_description(self.cluster_result.index)
        self.cluster_result.columns = ['cluster', 'gene', 'description']
        self.cluster_result.to_csv(outfile)
            
    def common_ancestor(self, genes):
        path = self.path_top(genes[0])
#        print "Depth of tree {}".format(len(path))
                             
        for gene in genes:
            path_new = self.path_top(gene)
            path = [ x for x in path if x in path_new ]
#        print "Depth: {}".format(len(path))
        return path[0]
    
    def is_child(self, child, parent):
        pos = child
        while (pos != self.parent_tree[pos] and pos != parent):
            pos = self.parent_tree[pos]
        if pos == parent:
            return True
        return False
    
    def get_children(self, node):
#        Cannot use recursion. Manage own queue
        node_list = [node]
        children_list = []
        while len(node_list) > 0:
            current_node = node_list.pop(0)
            if current_node < self.gene_num:
                children_list.append(current_node)
            else:
                node_list += [int(self.link_matrix[current_node - self.gene_num, 0]),
                              int(self.link_matrix[current_node - self.gene_num, 1])]
        return children_list
    
    def genes_to_idx(self, genes):
        gene = genes[0]
        if isinstance(gene, (int, long ,float, complex)):
            genes_idx = genes
        elif gene.startswith('ENS'):
            idx_list = np.arange(self.gene_num)
            genes_idx = [idx_list[np.where(self.normalized.index == x)] for x in genes]
            genes_idx = [x[0] for x in genes_idx]
        else:
            genes_list = self.annotation.get_id(genes)
            idx_list = np.arange(self.gene_num)
            genes_idx = [idx_list[np.where(self.normalized.index == x)] for x in genes_list]
            genes_idx = [x[0] for x in genes_idx if len(x) > 0 ]
        return genes_idx
    
    def wrapper_top_correlated_genes(self, gene, outfile, num_hit=50):
        gene_idx = self.genes_to_idx([gene])
        if not hasattr(self, 'corr_matrix'):
            corr_file = self.datafile + '.corr.npy'
            if os.path.isfile(corr_file) and os.path.getmtime(corr_file) >= os.path.getmtime(self.datafie):
                self.corr_matrix = np.load(corr_file)
            else:
                corr_matrix = pdist(self.normalized, metric='correlation')
                self.corr_matrix = squareform(corr_matrix)
                del corr_matrix
                np.save(corr_file, self.corr_matrix)
        gene_corr = self.corr_matrix[gene_idx][0]
        deco_gene = [ (x, i) for i, x in enumerate(gene_corr) ]
        deco_gene.sort()
        gene_idx = [ i for (x, i) in deco_gene ]
        gene_list = pd.DataFrame(self.normalized.index[gene_idx[:num_hit]])
        gene_list.columns = ['stable_id']
        gene_list = gene_list.set_index(['stable_id'])
        gene_list['gene'] = self.annotation.get_gene(gene_list.index)
        gene_list['description'] = self.annotation.get_description(gene_list.index)
        gene_list.to_csv(outfile)
        
    
    def wrapper_cluster_gene_names(self, genes):
        genes_idx = self.genes_to_idx(genes)
        ancestor_node = self.common_ancestor(genes_idx)
        gene_nodes = self.get_children(ancestor_node)
        stable_ids = self.normalized.index[gene_nodes]
        gene_names = self.annotation.get_gene(stable_ids)
        return gene_names
    
    def path_top(self, gene):
        path_list = [gene]
        pos = gene
        while (pos != self.parent_tree[pos]):            
            pos = self.parent_tree[pos]
            path_list.append(pos)
        return path_list
    
    def mrs_score(self, gene, ref):
        max_mrs = 0.0
        max_leafgroup = []
        for i in range(len(ref)):
            common_ancestor = self.common_ancestor([gene, ref[i]])
            leaf_nodes = self.get_children(common_ancestor)
            gene_of_interest = 0
            for ref_gene in ref:
#                if self.is_child(ref_gene, common_ancestor):
                if ref_gene in leaf_nodes:
                    gene_of_interest += 1
            mrs = 1.0 * gene_of_interest / self.leaf_num[common_ancestor]
            if mrs > max_mrs:
                max_mrs = mrs
                max_leafgroup = leaf_nodes
        return (max_mrs, max_leafgroup)
    
    def wrapper_get_mrs(self, ref_genes, out_file):
        self.mrs_score_array = np.array([0.0] * self.gene_num)
        self.mrs_group_array = [''] * self.gene_num
        ref_idx = self.genes_to_idx(ref_genes)
        for i in range(self.gene_num):
            (self.mrs_score_array[i], mrs_group) = self.mrs_score(i, ref_idx)
            self.mrs_group_array[i] = str(mrs_group)
        self.mrs_score_frame = self.normalized.iloc[:, :1]
        self.mrs_score_frame.iloc[:, 0] = self.mrs_score_array
        self.mrs_score_frame['gene'] = self.annotation.get_gene(self.mrs_score_frame.index)
        self.mrs_score_frame['description'] = self.annotation.get_description(self.mrs_score_frame.index)
        self.mrs_score_frame.columns = ['mrs','gene','description']
        self.mrs_score_frame['group'] = self.mrs_group_array
        self.mrs_score_frame.to_csv(out_file)
        
        
class gene_note:
    def __init__(self, datafile):
        self.gene_notes = pd.read_csv(datafile, header=None, index_col=0, sep='\t',
                                      names=['gene', 'xref', 'description'])
    
    def get_id(self, gene):
        if isinstance(gene, (list, tuple)):
            return [self.gene_notes.index[np.where(self.gene_notes['gene']==x)] for x in gene]
        else:
            return self.gene_notes.index[np.where(self.gene_notes['gene']==gene)]
            
    def get_gene(self, stable_id):
        return self.gene_notes.loc[stable_id, 'gene']
    
    def get_description(self, stable_id):
        return self.gene_notes.loc[stable_id, 'description']
    
    def get_external_ref(self, stable_id):
        return self.gene_notes.loc[stable_id, 'xref']
        
        
        
