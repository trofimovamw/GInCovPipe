#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 12:50:41 2020

@author: mariatrofimova
"""

# Infer global parameters and filter sites for origins counting
import math
import numpy as np
import random
import pandas as pd
from Bio import Phylo
from Bio import SeqIO
from scipy.stats import binom
import csv
from bitarray import bitarray


class parameterEstimation:
    
    def __init__(self, seqSets, digitSeqSets, reffile):
        self.seqSets = seqSets
        # Not needed anymore!
        self.digitSeqSets = digitSeqSets
        # Cutoff at 1/10 of all sequences in set
        self.freqCutoff = 2
        self.reffile = reffile

    
    # Count mutants
    def _countMutants(self,seqSets):
        mut_count = 0
        num_seqs = 0
        # Enumerate bins
        for t, seqSet in enumerate(seqSets):
            if len(seqSet)!=0:
                # Enumerate sequences
                for i, seq in enumerate(seqSet):
                    #print(seq)
                    if seq[1]!='':
                        mut_count += 1
                        num_seqs += 1
                    else:
                        num_seqs += 1
        return mut_count, num_seqs
    
    def _getAltBaseCounts(self,seqSets):
        # Dictionary with (key,value): "(Position>MutantBase,count)"
        positions = dict()
        for t, seqSet in enumerate(seqSets):
            if len(seqSet)!=0:
                for i, seq in enumerate(seqSet):
                    if seq[1]!='':
                        split_fp = seq[1].split("-")
                        for pos in split_fp:
                            split_pos = pos.split(">")
                            if int(split_pos[0]) in positions:
                                #positions[split_pos[0]] += 1
                                bases = positions[int(split_pos[0])]
                                if split_pos[1]=='A':
                                    bases[0] += 1
                                elif split_pos[1]=='C':
                                    bases[1] += 1
                                elif split_pos[1]=='G':
                                    bases[2] += 1
                                elif split_pos[1]=='T':
                                    bases[3] += 1
                                # If ambiguous base: TODO: check
                                else:
                                    bases[4] += 1
                                positions[int(split_pos[0])] = bases
                            else:
                                # A, C, G, T
                                bases = np.zeros(5)
                                if split_pos[1]=='A':
                                    bases[0] += 1
                                elif split_pos[1]=='C':
                                    bases[1] += 1
                                elif split_pos[1]=='G':
                                    bases[2] += 1
                                elif split_pos[1]=='T':
                                    bases[3] += 1
                                # If ambiguous base: TODO check
                                else:
                                    bases[4] += 1
                                positions[int(split_pos[0])] = bases
        return positions
    
    def _getAltPosCountsBinom(self,seq_binom):
        # Aggregate all mutant positions in a dictionary
        # (key,value):(position,mutant count)
        positions = dict()
        
        for i, seq in enumerate(seq_binom):
            for j, pos in enumerate(seq):
                if pos==1:
                    if j in positions:
                        positions[j] += 1
                    else:
                        positions[j] = 1
        return positions
    
    # Measure mutant position entropy
    def _calculateMutantEntropy(self,seq_binom):
        # Get mutant frequency and number of sequences
        #mut_count, num_seqs = self._countMutants(self.seqSets)
        # Get counts of mutant positions
        positions = self._getAltPosCounts(seq_binom)
        # (key,value): "(Position,mutant count)
        
        # Entropy dictionary - on mutant sites
        H = dict()
        for (key,value) in positions.items():
            # Counts of bases on position
            alt_count = value
            # Count of reference bases
            ref_base_count = num_seqs-sum(alt_count)
            # Get the entropy of position
            sum_entr = []
            # Variants array 
            alt_base_on_pos = [value/num_seqs,ref_base_count/num_seqs]
            for b in alt_base_on_pos:
                if b==0:
                    sum_entr.append(0)
                else:
                    lp = (b)*math.log(b)
                    sum_entr.append(lp)
            #print(sum_entr)
            # Sum of p*log(p) for mutation and no mutation
            hh = -sum(sum_entr)
            H[int(key)] = hh
        #print(H) 1- (1- 4.4e-7)^num_seq
        #         1- (1- 4.4e-7)^NrSEq
        return H
    
    def _proposeP(self):
        return random.uniform(math.exp(-10),math.exp(-5))
    
    # nseq: list of bin sizes
    # nbins: depr - amount of bins
    # lenseq: length of sequence
    # entropies: true measured entropies on positions
    def _runErrorRate(self,nseq,nbins,lenseq,entropies_r):
        # Initial probability fit
        p = 4.4*math.exp(-7)
        
        # Bin real entropies
        Harr = []
        for (key,value) in entropies_r.items():
            Harr.append(value)
        # TODO: Define edges instead of binning size, looking at mutations count
        bin_entropies = np.histogram(Harr, 100)
        bin_edges = bin_entropies[1]
        bin_probabilities = bin_entropies[0]/sum(bin_entropies[0])
        KL = 10
        eps = 0.1
        while KL>eps:
            for i in range(len(nbins)):
                p = self._proposeP()
                # Sequence set represented as 0-1 sequence
                seq_binom = []
                # Number of sequences in bin: nseq[i]
                for s in range(nseq[i]):
                    seq_binom.append(np.zeros(lenseq))
                # Draw number of mutations in bin
                binomial_draw = binom.pmf(p,nseq[i],1)
                fmut = binomial_draw/nseq[i]
                for j in range(len(seq_binom)):
                    for s in range(len(seq_binom[j])):
                        u = random.uniform(0,1)
                        if u<fmut:
                            seq_binom[j][i] = 1
                # Calculate entropy on positions
                # Get mutant entropies
                entropies = self._calculateMutantEntropy(seq_binom)
                # Bin proposed entropies based on bin edges of real data
                Hprop = []
                for (key,value) in entropies.items():
                    Harr.append(value)
                # Binned proposed entropies
                bin_sizes = []
                for e in range(1,len(bin_edges)):
                    bin_arr = []
                    for h in Hprop:
                        if e[i-1]<h<=e[i]:
                            bin_arr.append(h)
                    bin_sizes.append(len(bin_arr))
                binprop_probabilities = np.array(bin_sizes)/sum(bin_sizes)
                        
                # Calculate KL-divergence
                KLarr = []
                for b in range(len(bin_probabilities)):
                    if binprop_probabilities[b]!=0 and bin_probabilities[b]!=0:
                        div = binprop_probabilities[b]*np.log((binprop_probabilities[b]/bin_probabilities[b]))
                        KLarr.append(div)
                    else:
                        KLarr.append(0)
                KL = sum(KLarr)
        return p
    
    # Measure mutant base entropy
    def _calculateBaseEntropy(self,num_seqs):
        # Load reference
        ref_fasta = SeqIO.parse(open(self.reffile),'fasta')
        ref_seq = ''
        for fasta in ref_fasta:
            ref_seq = str(fasta.seq)
        # Get mutant frequency and number of sequences
        #mut_count, num_seqs = self._countMutants(self.seqSets)
        # Get counts of mutant positions
        positions = self._getAltBaseCounts(self.seqSets)
        # (key,value): "(Position>MutantBase,count)
        
        # Entropy dictionary - on mutant sites
        H = dict()
        for (key,value) in positions.items():
            # Counts of bases on position
            alt_base_on_pos = value
            # Count of reference bases
            ref_base_count = num_seqs-sum(alt_base_on_pos)
            # Base on ref
            ref_base = ref_seq[int(key)]
            # Put ref base count into the count array
            if ref_base=='A':
                alt_base_on_pos[0] += ref_base_count
            elif ref_base=='C':
                alt_base_on_pos[1] += ref_base_count
            elif ref_base=='G':
                alt_base_on_pos[2] += ref_base_count
            elif ref_base=='T':
                alt_base_on_pos[3] += ref_base_count
            # If it is an ambiguous base
            else:
                alt_base_on_pos[4] += ref_base_count
            # Get the entropy of position
            sum_entr = []
            for b in alt_base_on_pos:
                if b==0:
                    sum_entr.append(0)
                else:
                    lp = (b/sum(alt_base_on_pos))*math.log(b/sum(alt_base_on_pos))
                    sum_entr.append(lp)
            #print(sum_entr)
            # Sum of p*log(p) for mutation and no mutation
            hh = -sum(sum_entr)
            H[int(key)] = hh
        #print(H) 1- (1- 4.4e-7)^num_seq
        #         1- (1- 4.4e-7)^NrSEq
        return H
    
    def _getAltPosCounts(self):
        # Aggregate all mutant positions in a dictionary
        positions = dict()
        
        for t, seqSet in enumerate(self.seqSets):
            if len(seqSet)!=0:
                for i, seq in enumerate(seqSet):
                    if seq[1]!='':
                        split_fp = seq[1].split("-")
                        for pos in split_fp:
                            spl = pos.split(">")
                            pos = int(spl[0])
                            if pos in positions:
                                positions[pos] += 1
                            else:
                                positions[pos] = 1

        return positions
    
    # Simply removing positions with low prevalence
    def _filterSingletons(self):
        filteredSeqSets = []
        # First identify positions that pop up below a threshold

        positions = self._getAltPosCounts()
        #print(positions)
        #print(positions)
        
        # Identify not so frequent mutant positions
        minor_positions = []
        for (key,value) in positions.items():
            if value<self.freqCutoff:
                #print("Singleton position: ",key)
                #print("Count: ",value)
                minor_positions.append(key)
        #print("Positions below cutoff:")
        print(minor_positions)
        
        # Filter sequence fingerprint based on that
        for t, seqSet in enumerate(self.seqSets):
            filteredSets = []
            for i, seq in enumerate(seqSet):
                if seq[1]!='':
                    split_fp = seq[1].split("-")
                    to_keep = []
                    for pos in split_fp:
                        split_pos = pos.split(">")
                        posit = int(split_pos[0])
                        if not posit in minor_positions:
                            to_keep.append(pos)
                    merger = ''
                    if len(to_keep)!=0:
                        merger = "-".join(to_keep)
                    filteredSets.append((seq[0],merger))
                else:
                    filteredSets.append((seq[0],seq[1]))
            filteredSeqSets.append(filteredSets)
        
        return filteredSeqSets
    
    def _filterVariantsList(self):
        # Positions that are likely lab artifacts
        masking = ['187', '241', '335', '1059', '2094', '3037', '3130', 
                   '3145', '4050', '6255', '6990', '8022', '8782', '9223', 
                   '10323', '10741', '11074', '11083', '11704', '13402', 
                   '13408', '14408', '14724', '14786', '14805', '15324', 
                   '16887', '17247', '19684', '20148', '21137', '21575', 
                   '23403', '24034', '24378', '25563', '26144', '26461', 
                   '26681', '27384', '28077', '28826', '28854', '29353', '29700', '29736']
        filteredSeqSets = []
        # Filter sequence fingerprint based on that
        for t, seqSet in enumerate(self.seqSets):
            filteredSets = []
            for i, seq in enumerate(seqSet):
                if seq[1]!='':
                    split_fp = seq[1].split("-")
                    to_keep = []
                    for pos in split_fp:
                        split_pos = pos.split(">")
                        posit = int(split_pos[0])
                        if not posit in masking:
                            to_keep.append(pos)
                    merger = ''
                    if len(to_keep)!=0:
                        merger = "-".join(to_keep)
                    filteredSets.append((seq[0],merger))
                else:
                    filteredSets.append((seq[0],seq[1]))
            filteredSeqSets.append(filteredSets)
        return filteredSeqSets
    
    def _calculateMutantTrajectories(self,seqSets):
        trajectories = dict()
        
        for t, seqSet in enumerate(seqSets):
            for s, seq in enumerate(seqSet):
                if seq[1]!='':
                    split_fp = seq[1].split("-")
                    for pos in split_fp:
                        if pos in trajectories:
                            arr = trajectories[pos]
                            arr[t] += 1
                            trajectories[pos] = arr
                        else:
                            arr = np.zeros(len(seqSets))
                            arr[t] += 1
                            trajectories[pos] = arr
        for (key, value) in trajectories.items():
            traj = trajectories[key]
            for i in range(len(traj)):
                traj[i] = traj[i]/len(seqSets[i])
            trajectories[key] = traj
        return trajectories
    
    def _calculateAlleleTrajectories(self,seqSets):
        trajectories = dict()
        
        for t, seqSet in enumerate(seqSets):
            for s, seq in enumerate(seqSet):
                if seq[1]!='':
                    if seq[1] in trajectories:
                        arr = trajectories[seq[1]]
                        arr[t] += 1
                        trajectories[seq[1]] = arr
                    else:
                        arr = np.zeros(len(seqSets))
                        arr[t] += 1
                        trajectories[seq[1]] = arr
#        for (key, value) in trajectories.items():
#            traj = trajectories[key]
#            for i in range(len(traj)):
#                traj[i] = traj[i]/len(seqSets[i])
#            trajectories[key] = traj
        return trajectories

    
    # Removing mutated sites with low entropy
    def _filterMinorVariantsE(self,num_seqs):
        filteredSeqSets = []
        # First identify positions that pop up below a threshold
        entropies = self._calculateBaseEntropy(num_seqs)
        # Get the entropis dict
        rank_entropies = []
        for (key,value) in entropies.items():
            rank_entropies.append(value)
        # Sort entropies and make lower 5% cutoff
        rank_entropies = np.unique(sorted(rank_entropies))
        #print("Unique entropies: ",rank_entropies)
        le = len(rank_entropies)
        cutoff_ind = math.floor(le*0.05)
        cutoff = rank_entropies[cutoff_ind]
        #print("Enttropy cutoff: ",cutoff)

        # Filter sequence fingerprint based on that
        for t, seqSet in enumerate(self.seqSets):
            filteredSets = []
            for i, seq in enumerate(seqSet):
                if seq[1]!='':
                    split_fp = seq[1].split("-")
                    to_keep = []
                    for pos in split_fp:
                        split_pos = pos.split(">")
                        posit = int(split_pos[0])
                        if posit in entropies:
                            # If entropy of the position is above 5% cutoff - keep it
                            e = entropies[posit]
                            if e>cutoff:
                                to_keep.append(pos)
                    merger = ''
                    if len(to_keep)!=0:
                        merger = "-".join(to_keep)
                    filteredSets.append((seq[0],merger))
                else:
                    filteredSets.append((seq[0],seq[1]))
            filteredSeqSets.append(filteredSets)
        
        return filteredSeqSets, entropies
    
    # Removing mutated sites with low entropy
    def _filterMinorAlleles(self,num_seqs):
        filteredSeqSets = []
        # First identify positions that pop up below a threshold
        allele_trajectories = self._calculateAlleleTrajectories(self.seqSets)
        
        for t, seqSet in enumerate(self.seqSets):
            filteredSets = []
            for i, seq in enumerate(seqSet):
                if seq[1]!='':
                    if seq[1] in allele_trajectories:
                        arr = allele_trajectories[seq[1]]
                        freq = sum(arr)/num_seqs
                        if freq>0.001:
                            filteredSets.append(seq[1])
            filteredSeqSets.append(filteredSets)
        return filteredSeqSets
               
    # Proximity metric between two sequence fingerprints: size of intersection of the sets
  
    """
    Haplotype clustering
    """    
    # Proximity metric between two sequence fingerprints: size of intersection of the sets
    def _longestCommonInterval(self,d1,d2):
        pos_set1 = d1.split("-")
        pos_set2 = d2.split("-")
        
        a = bitarray(30000)
        a.setall(True)
        b = bitarray(30000)
        b.setall(True)
        
        for pos in pos_set1:
            posint = int(pos.split(">")[0])
            a[posint] = False
        for pos in pos_set2:
            posint = int(pos.split(">")[0])
            b[posint] = False
            if pos in pos_set1:
                a[posint] = True
                b[posint] = True
        match = a & b
        matches = []
        span = 0
        for i in range(len(match)):
            if match[i]==True:
                span += 1
            else:
                matches.append(span)
                span = 0
        matches.append(span)
        return span
    
    def _metric(self,d1,d2):
        longest_interval = self._longestCommonInterval(d1,d2)
        #first_pos = int(longest_interval[0].split(">")[0])
        #last_pos = int(longest_interval[-1].split(">")[0])
        #print("Longest matching interval: ",longest_interval)
        pos_set1 = d1.split("-")
        pos_set2 = d2.split("-")
        mismatches = set(pos_set1).difference(pos_set2)
        
        #print("Intersection: ",intersect)
        #difference = set(pos_set1).difference(set(pos_set2))
        #print("Difference: ",difference)
        #mismatches = [0 for k in difference]
        s1 = 30000-len(mismatches)
        print("Sequence pair ", d1, " \ ", d2)
        print(s1)
        #for mmatch in mismatches:
        #    s1 += mmatch
        #print("Component 1: ",s1)
        s2 = longest_interval
        if s2<0:
            s2 = 0
        s12 = s1+s2
        print(s2)
        #print("Component 2: ",s2)
        # The distance to itsef
        sii = 2*30000
        print("To itself:")
        print(sii)
        d = (sii-s12)/sii
        print(d)

        return d
        
    def _findNeighbours(self,data,d1,eps):
        nset = []
        for d in data:
            if d!=d1:
                #print("Calculating distance")
                dist = self._metric(d1,d)
                if dist<=eps:
                    nset.append(d)
        return set(nset)
    
    def _dbscan(self,data,minsize,eps):
        count = 0
        labels = dict()
        for d in data:
            # 0: no label
            labels[d] = 0
        for d in data:
            #print("Processing data point ",d)
            # Label is not yet assigned, i.e. it wasn't looked at before or was noise
            if labels[d] <= 0:
                # find neighbouring points
                neighbours = self._findNeighbours(data,d,eps)
                #print("Found neighbourhood:")
                #print(neighbours)
                if len(neighbours)<minsize:
                    # -1: noise
                    labels[d] = -1
                else:
                    count += 1
                    labels[d] = count
                    for n in neighbours:
                        if labels[n] == -1:
                            labels[n] = count
                        elif labels[n]>0:
                            print("Keeping it...")
                        else:
                            labels[n] = count
                            neighbours_ = self._findNeighbours(data,n,eps)
                            if len(neighbours)>=minsize:
                                neighbours = neighbours.union(neighbours_)
        return labels
    
    def _haplotypeClustering(self):
        allele_trajectories = self._calculateAlleleTrajectories(self.seqSets)
        #Identify a set of high frequency variants
        frequencies = dict()
        for (key,value) in allele_trajectories.items():
            sumcounts = sum(list(value))
            frequencies[key]= sumcounts
        #freq_distr = frequencies.values()
        
        # Keys of the dictionary are the allele types
        # Make dictionary with split fingerprint
        allele_trajectories_keys = allele_trajectories.keys()
        # Assign low frequency allele set to haplotype cluster
        labels = self._dbscan(allele_trajectories_keys,1,0.3)
        return labels
    
    def _filterHomoplasies(self):
        htable = pd.read_table("")
        positions = np.array(list(htable["Position"]))-1
        filteredSeqSets = []
        # Filter sequence fingerprint based on that
        for t, seqSet in enumerate(self.seqSets):
            filteredSets = []
            for i, seq in enumerate(seqSet):
                if seq[1]!='':
                    split_fp = seq[1].split("-")
                    to_keep = []
                    for pos in split_fp:
                        split_pos = pos.split(">")
                        posit = int(split_pos[0])
                        if  posit in positions:
                            to_keep.append(pos)
                    merger = ''
                    if len(to_keep)!=0:
                        merger = "-".join(to_keep)
                    filteredSets.append((seq[0],merger))
                else:
                    filteredSets.append((seq[0],seq[1]))
            filteredSeqSets.append(filteredSets)
        
        return filteredSeqSets

    def run(self):
        # In one of the trajectories - count mutants and number of sequences
        mut_count_pre, num_seqs_pre = self._countMutants(self.seqSets)
        print("Before filtering: ")
        print("   Number of sequences: ",num_seqs_pre)
        print("   Number of mutants: ",mut_count_pre)
        self.freqCutoff = math.floor(num_seqs_pre*0.02)
      
        filteredSets1 = self._filterSingletons()
        # Replace sequence set
        self.seqSets = filteredSets1
    
        mut_count_post, num_seqs_post = self._countMutants(filteredSets1)
        print("After filtering: ")
        print("   Number of sequences: ",num_seqs_post)
        print("   Number of mutants: ",mut_count_post)
      
        return mut_count_post/num_seqs_post, filteredSets1
    
    def runTree(self):
        mut_count_pre, num_seqs_pre = self._countMutants(self.seqSets)
        print("Before filtering: ")
        print("   Number of sequences: ",num_seqs_pre)
        print("   Number of mutants: ",mut_count_pre)
        self.freqCutoff = math.floor(num_seqs_pre*0.02)
        filteredSets1 = self._filterSingletons()
        self.seqSets = filteredSets1
        filteredSets2 = self._filterVariantsList()
        self.seqSets = filteredSets2
        #filteredSets2 = self._filterHomoplasies()
        mut_count_post, num_seqs_post = self._countMutants(filteredSets1)
        #labels = self._haplotypeClustering()
        return mut_count_post/num_seqs_post, self.seqSets
    
    def runSingles(self):
        mut_count_pre, num_seqs_pre = self._countMutants(self.seqSets)
        print("Before filtering: ")
        print("   Number of sequences: ",num_seqs_pre)
        print("   Number of mutants: ",mut_count_pre)
        self.freqCutoff = math.floor(2)
        filteredSets1 = self._filterSingletons()
        self.seqSets = filteredSets1
       
        mut_count_post, num_seqs_post = self._countMutants(filteredSets1)
        #labels = self._haplotypeClustering()
        return mut_count_post/num_seqs_post, self.seqSets
    
    def runMut(self):
        trajectories = self._calculateMutantTrajectories(self.seqSets)
        allele_trajectories = self._calculateAlleleTrajectories(self.seqSets)
        return trajectories, allele_trajectories
    
    def runHapl(self):
        labels = self._haplotypeClustering()
        return labels
    
    #def runError(self):
        #mut_count_pre, num_seqs_pre = self._countMutants(self.seqSets)
        #H = self._calculateBaseEntropy(self,num_seqs)
        # Nseq = array of binsizes
        # Nbins = number of bins
        # Lenseq = length of sequence
        # Entropies_r = H
        #p = self._runErrorRate(nseq,nbins,lenseq,entropies_r)                      
        
