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
from Bio import SeqIO
from scipy.stats import binom
from bitarray import bitarray


class parameterEstimation:
    
    def __init__(self, seqSets, digitSeqSets, reffile, freqCutoff):
        self.seqSets = seqSets
        self.digitSeqSets = digitSeqSets
        self.freqCutoff = freqCutoff
        self.reffile = reffile

    
    
    def _countMutants(self,seqSets):
        """
        Count mutant sequences in sample
        :param seqSets: set of bins of sequence fingerprints
        :return mut_count: list of mutant sequences counts per bin
        :return num_seqs: list of number of sequences per bin
        """
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
        """
        Count mutant positions in sample
        :param seqSets: set of bins of sequence fingerprints
        :return mut_count: list of mutant sequences counts per bin
        :return positions: dict of mutant positions and number of occurences 
        """
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
    
    def _calculateMutantEntropy(self,seq_binom):
        """
        Get entropies of mutant positions for error rate simulation
        :param seq_binom: set of bins of sequence fingerprints
        :return H: dictionary of mutant positions and entropies
        """
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
            # Sum of p*log(p) for mutation and no mutation
            hh = -sum(sum_entr)
            H[int(key)] = hh
       
        return H
    
    def _proposeP(self):
        """
        Proposal distribution for error rate
        """
        return random.uniform(math.exp(-10),math.exp(-5))
   
    def _runErrorRate(self,nseq,nbins,lenseq,entropies_r):
        """
        Get entropies of mutant positions for error rate simulation (in works)
        :param nseq: list of bin sizes
        :param nbins: number of bins
        :param lenseq: length of sequence
        :param entropies: true measured entropies on positions
        :return p: dictionary of mutant positions and entropies
        """
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
    
    def _calculateBaseEntropy(self,num_seqs):
        """
        Get entropies of mutant positions 
        :param num_seqs: number of sequences
        :return H: dictionary of mutant positions and entropies
        """
        # Load reference
        ref_fasta = SeqIO.parse(open(self.reffile),'fasta')
        ref_seq = ''
        for fasta in ref_fasta:
            ref_seq = str(fasta.seq)
        # Get counts of mutant positions
        positions = self._getAltBaseCounts(self.seqSets)        
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
            hh = -sum(sum_entr)
            H[int(key)] = hh
       
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
    
    def _filterSingletons(self):
        """
        Filter our mutant positions that only occur few times in entire sample
        frequency cutoff defined by user via config
        :return filteredSeqSets: sequence fingerprints without singletons
        """
        filteredSeqSets = []
        # First identify positions that pop up below a threshold

        positions = self._getAltPosCounts()
        
        # Identify not so frequent mutant positions
        minor_positions = []
        for (key,value) in positions.items():
            if value<=self.freqCutoff:
                #print("Singleton position: ",key)
                #print("Count: ",value)
                minor_positions.append(key)
        #print("Positions below cutoff:")
        
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
  
    def _calculateMutantTrajectories(self,seqSets):
        """
        Calculate trajectories of frequencied of mutant bases 
        :return trajectories: dict with key - mutant position and value - frequencies
            in individual bins
        """
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
        """
        Calculate trajectories of frequencied of mutant sequences 
        :return trajectories: dict with key - mutant sequence and value - frequencies
            in individual bins
        """
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

        return trajectories

    
    # Removing mutated sites with low entropy
    def _filterMinorVariantsE(self,num_seqs):
        """
        Removing mutated sites with low entropy
        :return filteredSeqSets: sequence fingerprints without low entropy 
            positions
        :return entropies: entropies of variants
        """
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
  
    """
    Haplotype clustering
    """    
    # Proximity metric between two sequence fingerprints: size of intersection of the sets
    def _longestCommonInterval(self,d1,d2):
        """
        Haplotype clustering metric - longest common subsequence
        :param d1: sequence fingerprint 1
        :param d2: sequence fingerprint 2
        :return span: length of longest common subsequence
        """  
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
        """
        Haplotype clustering metric
        :param d1: sequence fingerprint 1
        :param d2: sequence fingerprint 2
        :return d: float distance 
        """  
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
        """
        Clustering - find neighbours of current point in dbscan
        :param data: list of sequence fingerprints
        :param d1: center of neighbourhood
        :param eps: float - neighbourhood tolerance value
        """  
        nset = []
        for d in data:
            if d!=d1:
                #print("Calculating distance")
                dist = self._metric(d1,d)
                if dist<=eps:
                    nset.append(d)
        return set(nset)
    
    def _dbscan(self,data,minsize,eps):
        """
        Clustering with dbscan
        :param data: list of sequence fingerprints
        :param minsize: minimal size of neighbourhood to use
        :param eps: float - neighbourhood tolerance value
        :return labels: list of labels of each sequence - which cluster it belongs to
        """  
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
        """
        Haplotype clustering - find neighbours of current point in dbscan
        :return labels: list of labels of each sequence - which cluster it belongs to
        """  
        allele_trajectories = self._calculateAlleleTrajectories(self.seqSets)
        #Identify a set of high frequency variants
        frequencies = dict()
        for (key,value) in allele_trajectories.items():
            sumcounts = sum(list(value))
            frequencies[key]= sumcounts
        
        # Keys of the dictionary are the allele types
        # Make dictionary with split fingerprint
        allele_trajectories_keys = allele_trajectories.keys()
        # Assign low frequency allele set to haplotype cluster
        labels = self._dbscan(allele_trajectories_keys,1,0.3)
        return labels

    def run(self):
        """
        Run filter - remove variants that are below predefined cutoff 
        :return mut_count_post/num_seqs_post: mutant sequences proportion
        :return filteredSets1: sequences set filtered based on mutant positions
            frequency cutoff
        """ 
        # In one of the trajectories - count mutants and number of sequences
        mut_count_pre, num_seqs_pre = self._countMutants(self.seqSets)
        print("           Before filtering: ")
        print("              Number of sequences: ",num_seqs_pre)
        print("              Number of mutant sequences: ",mut_count_pre)
      
        filteredSets1 = self._filterSingletons()
        # Replace sequence set
        self.seqSets = filteredSets1
    
        mut_count_post, num_seqs_post = self._countMutants(filteredSets1)
        print("           After filtering: ")
        print("              Number of sequences: ",num_seqs_post)
        print("              Number of mutant sequences: ",mut_count_post)
      
        return mut_count_post/num_seqs_post, filteredSets1
    
    def runTree(self):
        """
        Run filter for tree-based analysis- remove variants that are below predefined cutoff
        :return mut_count_post/num_seqs_post: mutant sequences proportion
        :return filteredSets1: sequences set filtered based on mutant positions
            frequency cutoff
        """
        mut_count_pre, num_seqs_pre = self._countMutants(self.seqSets)
        print("           Before filtering: ")
        print("              Number of sequences: ",num_seqs_pre)
        print("              Number of mutant sequences: ",mut_count_pre)
        filteredSets1 = self._filterSingletons()
        self.seqSets = filteredSets1
        filteredSets2 = self._filterVariantsList()
        self.seqSets = filteredSets2
        mut_count_post, num_seqs_post = self._countMutants(filteredSets1)
        return mut_count_post/num_seqs_post, self.seqSets
    
    def runSingles(self):
        """
        Remove singletons only - ignore predefined cutoff
        :return mut_count_post/num_seqs_post: mutant sequences proportion
        :return filteredSets1: sequences set filtered based on mutant positions
            frequency cutoff
        """
        mut_count_pre, num_seqs_pre = self._countMutants(self.seqSets)
        print("           Before filtering: ")
        print("              Number of sequences: ",num_seqs_pre)
        print("              Number of mutant sequences: ",mut_count_pre)
        self.freqCutoff = 1
        filteredSets1 = self._filterSingletons()
        self.seqSets = filteredSets1
       
        mut_count_post, num_seqs_post = self._countMutants(filteredSets1)
        return mut_count_post/num_seqs_post, self.seqSets
    
    def runMut(self):
        """
        Run only calculation of overall mutant proportion
        :return trajectories: dictionary of mutant base trajectories
        :return allele_trajectories: dict of mutant sequences trajectories
        """
        trajectories = self._calculateMutantTrajectories(self.seqSets)
        allele_trajectories = self._calculateAlleleTrajectories(self.seqSets)
        return trajectories, allele_trajectories
    
    def runHapl(self):
        """
        Run haplotype clustering
        :return labels: list - number of cluster sequence belongs to
        """
        labels = self._haplotypeClustering()
        return labels
    
  