#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 18:50:49 2020

@author: mariatrofimova
"""

# Bootstrapping sequences in bins
import random

import numpy as np

import math

from modular_theta_from_dict import analyzeTrajectory

from parameter_est import parameterEstimation

class BootstrapFP:
    
    def __init__(self, seqSets, num_iter, delta_days, reference, mut_proportion):
        self.seqSets = seqSets
        self.num_iter = num_iter
        self.mut_proportion = mut_proportion
        self.delta_days = delta_days
        self.reference = reference

        
    def _callBootstrap(self):
        # Create a big bootstrap set to sample from
        # Resample individual bins in blocks
        bootstrapSeqSet = []
        sizes = []
        for t, seqSet in enumerate(self.seqSets):
            bootstrap = []
            #binsignatures = []
            sizes = []
            mutant_set = []
            if len(seqSet)!=0:
                for i,seq in enumerate(seqSet):
                    split = seq[1].split("-")
                    sizes.append(len(split))
                    for s in split:
                        if s!='':
                            mutant_set.append(s)
            mutants = list(set(mutant_set))
            #print(mutants)
            avrg_size = sum(sizes)/len(sizes)
            if len(seqSet)!=0:
                for i in range(len(seqSet)):
                    new_seq_n = math.floor(np.random.exponential(avrg_size))
                    new_set = random.choices(mutants,k=new_seq_n)
                    new_seq = '-'.join(np.unique(new_set))
                    bootstrap.append((seqSet[i][0],new_seq))
            bootstrapSeqSet.append(bootstrap)
                    
            #signatures.append(binsignatures)
#        # Average size
#        avrg_s = sum(sizes)/len(sizes)
#        # Make blocks - stationary bootstrap
#        block_signatures = []
#        for i in range(len(signatures)-1):
#            size = math.ceil(np.random.exponential(avrg_s))
#            if (i+size)<len(signatures):
#                block = signatures[i:(i+size)]
#                block_signatures.append(block)
#            else:
#                block = signatures[i:-1]
#                block_signatures.append(block)
#        # Iterate again to sample randomly from signatures
#        for t, seqSet in enumerate(self.seqSets):
#            size = len(seqSet)
#            bootstrap = []
#            if len(seqSet)!=0:
#                new_set = random.choices(block_signatures, k=1)[0]
#                for i, seq in enumerate(new_set):
#                    #print(bseq)
#                    bootstrap.append(seq)
#            bootstrapSeqSet.append(bootstrap)
        #print(bootstrapSeqSet)
        return bootstrapSeqSet
    
    def run(self):
        # Run simulations
        current_iter = 0
        estimates_theta = []
        estimates_mu = []
        while current_iter<self.num_iter:
            print("Resampling step: ",current_iter+1)
            bootstrapSeqSet = self._callBootstrap()
            #filt = parameterEstimation(bootstrapSeqSet,bootstrapSeqSet,self.reference)
            #mut_proportion, filtered_seqset = filt.run()
            #analyze = analyzeTrajectory(filtered_seqset, self.mut_proportion, self.reference)
            #thetas, variance, variance_size, num_seqs, num_mut, origins = analyze.analyzeBinsMLE()
            analyze1 = analyzeTrajectory(bootstrapSeqSet, self.mut_proportion, self.reference)
            thetas2, mus, num_seqs, variances = analyze1.analyzeBinsF()
            thetas, mu = analyze1.calculateParamsFromLS(thetas2,mus,self.delta_days)
            estimates_theta.append(thetas)
            estimates_mu.append(mu)
            current_iter += 1
        return estimates_theta, estimates_mu
            
            