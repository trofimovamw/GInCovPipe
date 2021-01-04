#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:53:44 2020

@author: mariatrofimova
"""
from Bio import Phylo

class treeLookup:
    
    def __init__(self,tree_path):
        self.tree_path = tree_path
        self.tree = Phylo.read(self.tree_path, "newick")
        
        
    def _get_parent(self,tree, child):
        node_path = tree.get_path(child)
        return node_path[-2]
        
    def _tabulate_names(self,tree):
        names = {}
        for idx, clade in enumerate(tree.find_clades(terminal=False)):
            if clade.name:
                clade.name = "ID:%d_%s" % (idx,clade.name)
            elif clade.confidence:
                clade.name = "ID:%d_%s" % (idx,clade.confidence)
            else:
                clade.name = "ID:%d" % (idx)
            names[clade.name] = clade
        return names

    # Tabulate names and label homoplasies if they were already encountered
    # Label: ID:idx_H1(sign)-H2(sign)...
    def _tabulate_names_wtrack(self):
        names = {}
        # Key pos: number-of-times-seen
        seen_labels = dict()
        # Go through non-terminal nodes
        for idx, clade in enumerate(self.tree.find_clades(terminal=False)):
            # if clade confidence present - one homoplasy
            if clade.confidence:
                if str(clade.confidence) in seen_labels:
                    seen_labels[str(clade.confidence)] += 1
                    clade.name = "ID:%d_%d(%d)" % (idx,clade.confidence,seen_labels[str(clade.confidence)])
                else:
                    seen_labels[str(clade.confidence)] = 1
                    clade.name = "ID:%d_%d(%d)" % (idx,clade.confidence,seen_labels[str(clade.confidence)])
            # If clade name is present - multiple homoplasies
            elif clade.name:
                cname = clade.name
                print("Clade name avail")
                print(cname)
                split_name = cname.split("-")
                name = "ID:%d" % (idx)
                
                for pos in split_name:
                    if str(pos) in seen_labels:
                        seen_labels[str(pos)] += 1
                        name = "%s_%s(%d)" % (name,pos,seen_labels[str(clade.confidence)])
                    else:
                        seen_labels[str(clade.confidence)] = 1
                        name = "%s_%s(%d)" % (name,pos,seen_labels[str(clade.confidence)])
                print("Transformed clade name")
                print(name)
                clade.name = name
            else:
                clade.name = "ID:%d" % (idx)
            names[clade.name] = clade
        return names
    
    def _makeClustersFromTree(self):
        paths = []
        clusters_dict = dict()
        
        for terminal in (self.tree.get_terminals()):
            name = terminal.name
            print(name)
            path_to_root = self.tree.get_path(name)[:-1]
            #print(path_to_root)
            path_to_root.reverse()
            paths.append(path_to_root)
            found = False
            for l in range(len(path_to_root)):
                if found==False:
                    # Find closest marker
                    clade_name = path_to_root[l].name
                    split_clade_name = clade_name.split("_")
                    #print(clade_name)
                    if (len(split_clade_name)>1):
                        clusters_dict[name] = clade_name
#                        if clade_name in clusters_dict.keys():
#                            arr = clusters_dict[clade_name]
#                            arr.append(name)
#                            clusters_dict[clade_name] = arr
#                        else:
#                            arr = [name]
#                            clusters_dict[clade_name] = arr
                        found = True
        return clusters_dict


    def run(self):
        names = self._tabulate_names_wtrack()
        clusters = self._makeClustersFromTree()
        return self.tree, names, clusters
    

    
    
    
    
    
    
    
    