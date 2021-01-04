#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 13:03:59 2020

@author: mariatrofimova
"""
from bitarray import bitarray
import numpy as np
import math 

class MSN:
    
    def __init__(self, seqSets, reference):
        self.reference = reference
        self.seqSets = seqSets
        self.V = []
        self.names_dict = dict()
        self.graph = self.buildGraph()
        self.adj = dict()
        
    def _longestCommonInterval(self,d1,d2):
        pos_set1 = d1.split("-")
        pos_set2 = d2.split("-")
        
        a = bitarray(29903)
        a.setall(True)
        b = bitarray(29903)
        b.setall(True)
        
        positions = []
        
        if pos_set1[0]!='':
            for pos in pos_set1:
                posint = int(pos.split(">")[0])
                positions.append(posint)
                a[posint] = False
        if pos_set2[0]!='':
            for pos in pos_set2:
                posint = int(pos.split(">")[0])
                b[posint] = False
                positions.append(posint)
                if pos in pos_set1:
                    a[posint] = True
                    b[posint] = True
        match = a & b
        posset = sorted(list(set(positions)))
        mismatches = [0]
        for i in range(len(posset)):
            val = match[posset[i]]
            if val == False:
                mismatches.append(posset[i])
        mismatches.append(29903)
        distances = np.diff(np.array(mismatches))
        
        return max(distances)
    
    def _kimura(self,base1,base2):
        # 0 - transition
        # 1 - transversion
        # 2 - ambiguous
        pair = base1+base2
        transitions = ["CT", "TC", "AG", "GA"]
        transversions = [ "AC", "CA", "AT", "GT", "TG", "TA", "GC", "CG" ]
        if pair in transitions:
            return 0
        elif pair in transversions:
            return 1
        else: 
            return 2
        
    def _metricH(self,d1,d2):
        dist = 0
        if len(d1)!=0 and len(d2)!=0:
            # Get positions from FPs and reference
            pos_set1 = d1.split("-")
            pos_set2 = d2.split("-")
            # In pos_set1 and not in pos_set2
            difference1 = set(pos_set1).difference(set(pos_set2))
            # In pos_set2 and not in pos_set1
            difference2 = set(pos_set2).difference(set(pos_set1))
            # If this mutant only in pos_set1 - pos_set2 has reference
            # Transition count 
            ts = 0
            # Transverion count
            tv = 0
            for pos in difference1:
                #print("Difference 1 ",difference1)
                position = int(pos.split(">")[0])
                #print(pos.split(">")[1])
                otherside = self.reference[position-1]
                #print(otherside)
                kimura = self._kimura(pos.split(">")[1],otherside)
                if kimura==0:
                    ts += 1
                elif kimura==1:
                    tv += 1
            for pos in difference2:
                #print("Difference 2 ",difference2)
                position = int(pos.split(">")[0])
                #print(pos.split(">")[1])
                otherside = self.reference[position-1]
                #print(otherside)
                kimura = self._kimura(pos.split(">")[1],otherside)
                if kimura==0:
                    ts += 1
                elif kimura==1:
                    tv += 1
            p = ts/29903
            q = tv/29903
            dist = -0.5*math.log((1-2*p-q)*math.sqrt(1-2*q))
            #print(d1,d2)
            #print("Kimura distance: ",dist)
        elif len(d1)!=0 and len(d2)==0:
            # d2 is reference
             # Transition count 
            ts = 0
            # Transverion count
            tv = 0
            for pos in d1.split("-"):
                position = int(pos.split(">")[0])
                otherside = self.reference[position-1]
                kimura = self._kimura(pos.split(">")[1],otherside)
                if kimura==0:
                    ts += 1
                elif kimura==1:
                    tv += 1
            p = ts/29903
            q = tv/29903
            dist = -0.5*math.log((1-2*p-q)*math.sqrt(1-2*q))
            #print(d1,d2)
            #print("Kimura distance: ",dist)
        elif len(d1)==0 and len(d2)!=0:
            # d2 is reference
             # Transition count 
            ts = 0
            # Transverion count
            tv = 0
            for pos in d2.split("-"):
                position = int(pos.split(">")[0])
                otherside = self.reference[position-1]
                kimura = self._kimura(pos.split(">")[1],otherside)
                if kimura==0:
                    ts += 1
                elif kimura==1:
                    tv += 1
            p = ts/29903
            q = tv/29903
            dist = -0.5*math.log((1-2*p-q)*math.sqrt(1-2*q))
            #print(d1)
            #print(d2)
            #print("Kimura distance: ",dist)
        else:
            dist = 0
        return dist
            
        
    def _metric(self,d1,d2):
        ref = self.reference
        longest_interval = self._longestCommonInterval(d1,d2)
        
        pos_set1 = d1.split("-")
        pos_set2 = d2.split("-")
        leftmm = set(pos_set1).difference(set(pos_set2))
        rightmm = set(pos_set2).difference(set(pos_set1))
        mismatches = len(leftmm.union(rightmm))#/(len(pos_set1)+len(pos_set2))
        # Jaccard coefficient
        leftm = set(pos_set1).intersection(set(pos_set2))
        rightm = set(pos_set2).intersection(set(pos_set1))
        matchesJ = len(leftm.union(rightm))/(len(set(pos_set1).union(set(pos_set2))))
        
        s1 = 29903-mismatches
        
        s2 = longest_interval
        #if s2<0:
        #    s2 = 0
        s12 = s1+s2
        
        sii = 2*29903
       
        d = (sii-s12)/sii

        return d
    
    def find(self, parent, i):
        if parent[i] == i:
            return i
        return self.find(parent, parent[i])
 
    def connect(self, parent, r, u, v):
        uroot = self.find(parent, u)
        vroot = self.find(parent, v)
        
        if r[uroot] > r[vroot]:
            parent[uroot] = uroot
        elif r[uroot] < r[vroot]:
            parent[uroot] = vroot
        else:
            parent[vroot] = uroot
            r[uroot] += 1
        
    def buildGraph(self):
        graph = []
        vertices = []
        names_dict = dict()
        ident = 0
        for t, seqSet in enumerate(self.seqSets):
            for i, seq in enumerate(seqSet):
                vertices.append(seq)
        for i in range(len(vertices)-1):
            for j in range(i+1,len(vertices)):
                identi = 0
                identj = 0
                if not vertices[i][0] in names_dict:
                    names_dict[vertices[i][0]] = ident
                    ident = ident+1
                    identi = names_dict[vertices[i][0]]
                else: 
                    identi = names_dict[vertices[i][0]]
                    
                if not vertices[j][0] in names_dict:
                    names_dict[vertices[j][0]] = ident
                    ident = ident+1
                    identj = names_dict[vertices[j][0]]
                else:
                    identj = names_dict[vertices[j][0]]
                d = self._metric(vertices[i][1],vertices[j][1])
                graph.append([identi,identj,d])
                
        self.V = len(vertices)
        self.graph = graph
        self.names_dict = names_dict
            
    def KruskalMST(self):
        print("Calculating minimal spanning network")
        result = []
        i = 0
        e = 0
        self.graph =  sorted(self.graph, key=lambda item: item[2])
        parent, r = [], []
        for node in range(self.V):
            parent.append(node)
            r.append(0)
        while e < self.V-1:	
            u, v, w =  self.graph[i]
            i += 1
            x = self.find(parent, u)
            y = self.find(parent, v)
            if x != y:
                e += 1  
                result.append([u, v, w])
                self.connect(parent, r, x, y)                 
        return result, self.names_dict
    
    def dfs(self,vertices, v, visited):
        visited[v] = True
        vertices.append(v)
        for c in self.adj[v]:
            if not visited[c]:
                vertices = self.dfs(vertices,c,visited)
        return vertices
    
    def KruskalClustering(self,edges):
        print("Old edges: ",len(edges))
        # Create clusters by splitting at longest edges 
        new_edges = []
        new_adj = dict()
        rem_adj = dict()
        cutoff = 0.05
        for e in edges:
            if e[2] < cutoff:
                new_edges.append(e)
                if e[0] in new_adj.keys():
                    new_adj[e[0]].append(e[1])
                else:
                    new_adj[e[0]] = []
                    new_adj[e[0]].append(e[1])
                if e[1] in new_adj.keys():
                    new_adj[e[1]].append(e[0])
                else:
                    new_adj[e[1]] = []
                    new_adj[e[1]].append(e[0])
            else:
                if not e[0] in new_adj.keys():
                    new_adj[e[0]] = []
                if not e[1] in new_adj.keys():
                    new_adj[e[1]] = []
                if not e[0] in rem_adj.keys():
                    rem_adj[e[0]] = []
                    rem_adj[e[0]].append(e[1])
                if not e[1] in rem_adj.keys():
                    rem_adj[e[1]] = []
                    rem_adj[e[1]].append(e[0])
        print("New edges: ",len(new_edges))
        # Find connected components
        visited = [False for v in range(self.V)]
        components = []
        self.adj = new_adj
        for v in range(self.V):
            if visited[v] == False:
                vertices = []
                components.append(self.dfs(vertices,v,visited))
        print("Connected components with ",cutoff)
        components_dict = dict()
        
        for i in range(len(components)):
            for j in range(len(components[i])):
                components_dict[components[i][j]] = i
        
        # Clustered network
        clusters_edges = []
        for key,value in self.adj.items():
            fr = components_dict[key]
            new_value = []
            for v in value:
                nv = components_dict[v]
                new_value.append(nv)
            new_value_u = np.unique(np.array(new_value))
            for nvu in new_value_u:
                clusters_edges.append([fr,nvu])
        
        return components, components_dict, clusters_edges
        
  
    def getOriginsPerMutant(self,cl_edges,cl_nodes,inv_ids_dict,components_fps):
        adj_list = dict()
        for e in cl_edges:
            # Get sequence fps belonging to cluster
            fr = e[0]
            to = e[1]
            
            # Make adjacency list
            if fr in adj_list.keys():
                adj_list[fr].append(to)
            else: 
                adj_list[fr] = []
                adj_list[fr].append(to)
            if to in adj_list.keys():
                adj_list[to].append(fr)
            else: 
                adj_list[to] = []
                adj_list[to].append(fr)
        # Get mutants that are in cluster and not in its neighbours
        # key: cluster, value: isolated positions
        isolated_mutants = dict()
        for key,value in adj_list.items():
            # Fps in this cluster
            key_fps = components_fps[key]
            key_positions = []
            for k in key_fps:
                pos = k.split("-")
                for p in pos:
                    key_positions.append(p)
            u_key_positions = np.unique(np.array(key_positions))
            # Fps in disconnected clusters
            n_positions = []
            for v in value:
                nfps = components_fps[v]
                for fp in nfps:
                    pos = fp.split("-")
                    for p in pos:
                        n_positions.append(p)
            u_n_positions = np.unique(np.array(n_positions))
            # Positions in key and not in neighbours
            difference = set(u_key_positions).difference(set(u_n_positions))
            isolated_mutants[key] = []
            if len(difference)!=0:
                for pos in list(difference):
                    isolated_mutants[key].append(pos)
                    
        mutants_origins = dict()
        mutants_w_multiple_origins = dict()
        for key,value in isolated_mutants.items():
            for v in value:
                if not v in mutants_origins.keys():
                    mutants_origins[v] = []
                    mutants_origins[v].append(key)
                else:
                    mutants_origins[v].append(key)
                    if not v in mutants_w_multiple_origins.keys():
                        mutants_w_multiple_origins[v] = mutants_origins[v]
                    else:
                        mutants_w_multiple_origins[v].append(key)
        return mutants_w_multiple_origins
    
    def dfsc(self,vertices, v, visited, adj):
        visited[v] = True
        vertices.append(v)
        for c in adj[v]:
            print(c)
            if visited[c]==False:
                vertices = self.dfsc(vertices,c,visited,adj)
        return vertices
    
    def getConnectedMutClusters(self,cl_edges,cl_nodes,cl_fps):
        adj_list = dict()
        fps = []
        for e in cl_edges:
            # Get sequence fps belonging to cluster
            fr = e[0]
            fr_fps = cl_fps[e[0]]
            split_fr = fr_fps.split("-")
            to = e[1]
            to_fps = cl_fps[e[1]]
            split_to = to_fps.split("-")
            fps.append(split_fr)
            fps.append(split_to)
            
            # Make adjacency list
            if fr in adj_list.keys():
                adj_list[fr].append(to)
            else: 
                adj_list[fr] = []
                adj_list[fr].append(to)
            if to in adj_list.keys():
                adj_list[to].append(fr)
            else: 
                adj_list[to] = []
                adj_list[to].append(fr)
        fps_merged = [item for subset in fps for item in subset]
        fps_unique = np.unique(np.array(fps_merged))
        # For each unique FPs find connected components that contain the mutation
        fps_components_clusters = []
        for fp in fps_unique:
            # New nodes
            nodes = []
            # New adj list
            indiv_adj = dict()
            for i in range(len(cl_fps)):
                # If node FP contains the mutation - keep the node
                if cl_fps[i].find(fp) != -1:
                    nodes.append(i)
            # For edges in adj list
            for key,value in adj_list.items():
                # If node ID in nodes with FP - keep it and edges to other nodes with FP
                if key in nodes:
                    new_value = list(set(value) & set(nodes))
                    indiv_adj[key] = new_value
            visitedv = [False for v in range(len(nodes))]
            visited = dict(zip(nodes, visitedv))
            print(visited)
            components = []
            for v in nodes:
                if visited[v] == False:
                    vertices = []
                    # Component dfs util
                    components.append(self.dfsc(vertices,v,visited,indiv_adj))
            fps_components_clusters.append(components)
        return fps_unique, fps_components_clusters
        
        
            
        
        
        
 

            
            
            
            
        