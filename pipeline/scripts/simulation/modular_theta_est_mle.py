#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 11:56:45 2020

@author: Maria Trofimova
"""

"""
Metrics from bins
"""
import numpy as np

import random

import math

from scipy import optimize


class analyzeTrajectory:

    def __init__(self, trajectory, initSeq):
        self.trajectory = trajectory
        self.initSeq = initSeq

    """
    Make origins from bins
    """
    def _makeOriginsHIDE(self, mutSeq, mutantsCount):
        '''
        Find first occurence of mutant - origins
        '''
        origins = []

        firstOcc = np.zeros(len(mutSeq))

        for i,(key,val) in enumerate(mutSeq.items()):
            ind = np.where(val==1)[0][0]
            firstOcc[i] = ind

        for i in range(len(mutantsCount)):
            c = len(np.where(firstOcc==(i))[0])
            origins.append(c)

        return origins

    """
    Make origins from bins - no track of seen mutants
    """
    def _findPosition(self, index, occList):
        # List of indices of lists where the index was found (bleh)
        bin_indices = []
        for i, pos in enumerate(occList):
            val_list = np.where(pos==index)[0]
            if val_list.size!=0:
                bin_indices.append(i)
        return bin_indices

    def _makeOrigins(self, mutSeq, mutantsCount):
        '''
        Find ALL occurences of mutant - originsAll
        '''
        origins = []
        occurence = [[] for i in range(len(mutSeq))]
        for i,(key,val) in enumerate(mutSeq.items()):
            ind = np.where(val==1)[0]
            occurence[i] = ind

        for i in range(len(mutantsCount)):
            ind_list = self._findPosition(i,occurence)
            origins.append(len(ind_list))

        return origins

    """
    Function to optimize
    """
    def _fmle(self,x,nu,ns):
        '''
        Function of the estimate for root finding
        '''
        # print("Iteration: number of origins: ", nu)
        # print("           number of mutants: ", ns)
        # print("           theta:             ", x)

        a = float(x*np.log(1+ns/x))
        #print("a: ",a)
        b = np.math.factorial(nu)
        #print("b: ",b)
        c = np.exp(-(x*np.log(1+ns/x)))
        #print("c: ",c)
        d = nu*math.log(a)
        dd = math.log(b)
        e = x*math.log(1+ns/x)
        return -(d-dd-e)

    """
    Optimization routine to call in each sampling method
    """
    def _optimize(self, nu, ns):
        # Constraint
        con1 = {'type': 'ineq', 'fun': lambda x: x}
        # If no data available -- set estimate to zero
        if nu==0 or ns==0:
            return 0
        # Else optimize with trust-constr option to keep evaluation in feasible region
        else:
            x0 = 1
            sol = optimize.minimize(self._fmle, x0=x0, method='trust-constr', args=(nu,ns), constraints=con1)
            return sol.x[0]

    """
    Analyze - no subsampling
    """
    def analyzeBinsNS(self):
        # Number of (mutated) sequences
        num_mut = []
        num_seqs = []
        # Dictionary of mutants
        mutSeqDict = dict()

        #Length of evolution
        traj_length = len(self.trajectory)

        # Filling the mutants dict and presence array
        for t, seqSet in enumerate(self.trajectory):
            mut_count = 0
            num_seqs.append(len(seqSet))
            for i in range(len(seqSet)):
                if (seqSet[i]!=self.initSeq):
                    mut_count += 1
                    if not seqSet[i] in mutSeqDict:
                        times = np.zeros(traj_length)
                        times[t] = 1
                        mutSeqDict[seqSet[i]] = times
            num_mut.append(mut_count)


        origins = self._makeOrigins(mutSeqDict, num_mut)
        '''
        Estimate effective population size from data
        '''
        thetas = []
        # For the initial generation
        thetas.append(0)
        #Variance
        variance = []
        ## instead of zero, put 1/N
        #variance.append(1)

        for i in range(1,len(origins)):
            sol = self._optimize(origins[i], num_mut[i])
            thetas.append(sol)

        days_since = np.zeros(len(thetas))
        curr_d = 0
        for i in range(len(days_since)):
            variance.append(1/num_seqs[i])
            days_since[i] = curr_d
            curr_d = curr_d + 7

        # for i in range(1,len(thetas)):
        #     if thetas[i]!=0:
        #         var = thetas[i]*math.log(1+num_mut[i]/thetas[i])
        #         variance.append(var)
        #     else:
        #         variance.append(1)

        return days_since, thetas, variance, num_seqs, origins

    """
    Analyze - with subsampling
    """
    def analyzeBinsWS(self):
        # Number of (mutated) sequences
        num_mut = []
        num_seqs = []
        # Dictionary of mutants
        mutSeqDict = dict()

        #Length of evolution
        traj_length = len(self.trajectory)

        # Sizes of subsamples
        subsample_sizes = []

        # Variance
        variance = []

        # Filling the mutants dict and presence array
        for t, seqSet in enumerate(self.trajectory):
            #Subsample - at least 4 sequences of the whole sample
            #if set is smaller or equal to 4, then take the whole set
            subsample_size = len(seqSet)
            if len(seqSet) > 4:
                while (subsample_size<=4 or subsample_size==len(seqSet)):
                    subsample_size = random.randint(0,len(seqSet))
            else:
                subsample_size = len(seqSet)
            # Subsampled set
            print("Subsampled ",subsample_size," of ",len(seqSet)," sequences.")
            sseqSet = random.sample(seqSet,subsample_size)

            if not len(sseqSet) == 0:
                mut_count = 0
                num_seqs.append(len(seqSet))
                subsample_sizes.append(subsample_size)
                variance.append(1 / subsample_size)
                for i in range(len(sseqSet)):
                    if (sseqSet[i]!=self.initSeq):
                        mut_count += 1
                        if not sseqSet[i] in mutSeqDict:
                            times = np.zeros(traj_length)
                            times[t] = 1
                            mutSeqDict[sseqSet[i]] = times
                num_mut.append(mut_count)
            else:
                break


        origins = self._makeOrigins(mutSeqDict, num_mut)
        '''
        Estimate effective population size from data
        '''
        thetas = []
        # For the initial generation
        thetas.append(0)

        for i in range(1,len(origins)):
            sol = self._optimize(origins[i], num_mut[i])
            thetas.append(sol)

        days_since = np.zeros(len(thetas))
        curr_d = 0
        for i in range(len(days_since)):
            days_since[i] = curr_d
            curr_d = curr_d + 7

        # for i in range(1,len(thetas)):
        #     if thetas[i]!=0:
        #         var = thetas[i]*math.log(1+num_mut[i]/thetas[i])
        #         variance.append(var)
        #     else:
        #         variance.append(1)


        return days_since, thetas, variance, num_seqs, subsample_sizes, origins

    """
    Analyze - subsampling equal amount of sequences
    """
    def analyzeBinsES(self):
        # Number of (mutated) sequences
        num_mut = []
        num_seqs = []
        # Dictionary of mutants
        mutSeqDict = dict()

        #Length of evolution
        traj_length = len(self.trajectory)

        # Sizes of subsamples
        subsample_sizes = []

        # Variance
        variance = []

        # Filling the mutants dict and presence array
        for t, seqSet in enumerate(self.trajectory):
            #Subsample - at least 4 sequences of the whole sample
            #if set is smaller or equal to 4, then take the whole set
            subsample_size = 100

            if len(seqSet)<subsample_size:
                subsample_size = len(seqSet)-1
            else:
                subsample_size = 50
            # Subsampled set
            print("Subsampled ",subsample_size," of ",len(seqSet)," sequences.")
            sseqSet = random.sample(seqSet,subsample_size)

            if not len(sseqSet) == 0:
                mut_count = 0
                num_seqs.append(len(seqSet))
                subsample_sizes.append(subsample_size)
                variance.append(1 / subsample_size)
                for i in range(len(sseqSet)):
                    if (sseqSet[i]!=self.initSeq):
                        mut_count += 1
                        if not sseqSet[i] in mutSeqDict:
                            times = np.zeros(traj_length)
                            times[t] = 1
                            mutSeqDict[sseqSet[i]] = times
                num_mut.append(mut_count)
            else:
                break


        origins = self._makeOrigins(mutSeqDict, num_mut)
        '''
        Estimate effective population size from data
        '''
        thetas = []
        # For the initial generation
        thetas.append(0)

        for i in range(1,len(origins)):
            sol = self._optimize(origins[i], num_mut[i])
            thetas.append(sol)

        days_since = np.zeros(len(thetas))
        curr_d = 0
        for i in range(len(days_since)):
            days_since[i] = curr_d
            curr_d = curr_d + 7

        # for i in range(1,len(thetas)):
        #     if thetas[i]!=0:
        #         var = thetas[i]*math.log(1+num_mut[i]/thetas[i])
        #         variance.append(var)
        #     else:
        #         variance.append(1)


        return days_since, thetas, variance, num_seqs, subsample_sizes, origins

    """
    Analyze - merge bins
    """
    def analyzeBinsMB(self, N_para):
        # Number of (mutated) sequences
        num_mut = []
        num_seqs = []
        bin_sizes = []
        # Dictionary of mutants
        mutSeqDict = dict()

        N = 0

        # Infer the average date of the bin
        days_since = []

        #Merge trajectory by N bins
        m_trajectory = []
        for i, seqSet in enumerate(self.trajectory):
            if N == 0:
                m_trajectory.append([])
                days_since.append(i * len(seqSet))
                # If N is zero, take random bin sizes between 1 and 5
                if N_para == 0:
                    N = random.randint(1, 5)
                else:
                    N = N_para
                bin_sizes.append(N)
                for seq in seqSet:
                    m_trajectory[-1].append(seq)
            else:
                days_since[-1] += i * len(seqSet)
                for seq in seqSet:
                    m_trajectory[-1].append(seq)
            N=N-1
        # correct the last bin size, if N was bigger than the remaining number of time steps
        bin_sizes[-1] = bin_sizes[-1] - N

        #Length of evolution
        traj_length = len(m_trajectory)

        # Filling the mutants dict and presence array
        for t, seqSet in enumerate(m_trajectory):
            mut_count = 0
            num_seqs.append(len(seqSet))

            # normalize bin size by number of sequences in bin
            days_since[t] = round(days_since[t]/len(seqSet)*7)

            for i in range(len(seqSet)):
                if (seqSet[i]!=self.initSeq):
                    mut_count += 1
                    if not seqSet[i] in mutSeqDict:
                        times = np.zeros(traj_length)
                        times[t] = 1
                        mutSeqDict[seqSet[i]] = times

            num_mut.append(mut_count)


        origins = self._makeOrigins(mutSeqDict, num_mut)
        '''
        Estimate effective population size from data
        '''
        thetas = []
        # For the initial generation
        #Variance
        variance = []

        for i in range(len(origins)):
            sol = self._optimize(origins[i], num_mut[i])
            # normalise by bin size
            thetas.append(sol / bin_sizes[i])
            variance.append(1/num_seqs[i])
            # if thetas[i]!=0:
            #     var = thetas[i]*math.log(1+num_mut[i]/thetas[i])
            #     variance.append(var)
            # else:
            #     variance.append(1)


        return days_since, thetas, variance, num_seqs, bin_sizes, origins

    """
    Analyze - merge bins with subsampling
    """

    def analyzeBinsWSMB(self, N_para):
        # Number of (mutated) sequences
        num_mut = []
        num_seqs = []
        # sampled bin sizes
        bin_sizes = []
        # Sizes of subsamples
        subsample_sizes = []


        # Dictionary of mutants
        mutSeqDict = dict()

        N = 0


        # Infer the average date of the bin
        # For WSMB: collect all dates and retrieve tha subsampled sequence date in the next loop
        seq_dates = []

        # Merge trajectory by N bins
        m_trajectory = []
        for i, seqSet in enumerate(self.trajectory):
            if N == 0:
                m_trajectory.append([])
                seq_dates.append([])
                # If N is zero, take random bin sizes between 1 and 5
                if N_para == 0:
                    N = random.randint(1, 5)
                else:
                    N = N_para
                bin_sizes.append(N)
                for seq in seqSet:
                    m_trajectory[-1].append(seq)
                    seq_dates[-1].append(i)
            else:
                for seq in seqSet:
                    m_trajectory[-1].append(seq)
                    seq_dates[-1].append(i)
            N=N-1

        # correct the last bin size, if N was bigger than the remaining number of time steps
        bin_sizes[-1] = bin_sizes[-1]-N

        # Length of evolution
        traj_length = len(m_trajectory)

        # Infer the average date of the bin
        days_since = np.zeros(traj_length)

        #Variance
        variance = []

        # Filling the mutants dict and presence array
        for t, seqSet in enumerate(m_trajectory):
            #Subsample - at least 4 sequences of the whole sample
            #if set is smaller or equal to 4, then take the whole set
            subsample_size = len(seqSet)

            if len(seqSet) > 4:
                while (subsample_size<=4 or subsample_size==len(seqSet)):
                    subsample_size = random.randint(0,len(seqSet))
            else:
                subsample_size = len(seqSet)
            # Subsampled set
            print(" Subsampled ",subsample_size," of ",len(seqSet)," sequences.")
            index_value = random.sample(list(enumerate(seqSet)),subsample_size)
            sseqSet = []
            for idx, val in index_value:
                days_since[t] += seq_dates[t][idx]*7
                sseqSet.append(val)

            # compute average date of subsampled bin
            days_since[t] = round(days_since[t]/subsample_size)

            if not len(sseqSet) == 0:
                mut_count = 0
                num_seqs.append(len(seqSet))
                subsample_sizes.append(subsample_size)
                variance.append(1 / subsample_size)
                for i in range(len(sseqSet)):
                    if (sseqSet[i]!=self.initSeq):
                        mut_count += 1
                        if not sseqSet[i] in mutSeqDict:
                            times = np.zeros(traj_length)
                            times[t] = 1
                            mutSeqDict[sseqSet[i]] = times
                num_mut.append(mut_count)
            else:
                break


        origins = self._makeOrigins(mutSeqDict, num_mut)
        '''
        Estimate effective population size from data
        '''
        thetas = []

        for i in range(len(origins)):
            sol = self._optimize(origins[i], num_mut[i])
            #normalise by bin size
            thetas.append(sol/bin_sizes[i])

        # for i in range(len(thetas)):
        #     if thetas[i]!=0:
        #         var = thetas[i]*math.log(1+num_mut[i]/thetas[i])
        #         variance.append(var)
        #     else:
        #         variance.append(1)


        return days_since, thetas, variance, num_seqs, bin_sizes, subsample_sizes, origins
