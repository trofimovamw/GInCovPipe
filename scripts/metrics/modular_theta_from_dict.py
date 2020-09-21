#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 11:56:45 2020

@author: Maria Trofimova

Based on: Bhavin S Khatri, Austin Burt, Robust Estimation of Recent
Effective Population Size from Number of Independent Origins in Soft
Sweeps, Molecular Biology and Evolution,
Volume 36, Issue 9, September 2019,
Pages 2040–2052, https://doi.org/10.1093/molbev/msz081
"""

"""
Metrics from bins
"""
import numpy as np

import math

from scipy import optimize


class analyzeTrajectory:

    def __init__(self, dict_traj, initSeq):
        self.dict_traj = dict_traj
        self.initSeq = initSeq

    """
    Make origins from bins
    """
    def _makeOrigins(self, mutSeq, mutantsCount):
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



    def _makeOriginsHIDE(self, mutSeq, mutantsCount):
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
    def _f(self,x,nu,ns):
        '''
        Function of the estimate for root finding
        '''
        return (x*np.log(1+ns/x)-nu)**2+(1/13)*(x**2)
        # Squared function

    """
    Analyze bins - sequences list as tuples: (name, positions_set)
    """
    def analyzeBins(self):
        # Number of (mutated) sequences
        num_mut = []
        num_seqs = []
        # Dictionary of mutants
        mutSeqDict = dict()

        #Length of evolution
        traj_length = len(self.dict_traj)

        #Variance
        variance = []
        # Filling the mutants dict and presence array
        for t, seqSet in enumerate(self.dict_traj):
            mut_count = 0
            num_seqs.append(len(seqSet))
            variance.append(1/len(seqSet))
            for i in range(len(seqSet)):
                if (seqSet[i][1]!=[]):
                    mut_count += 1
                    if not seqSet[i][1] in mutSeqDict:
                        times = np.zeros(traj_length)
                        times[t] = 1
                        mutSeqDict[seqSet[i][1]] = times
                    else:
                        times = mutSeqDict[seqSet[i][1]]
                        if times[t]==1:
                            pass
                        elif times[t]==0:
                            times[t] = 1
                            mutSeqDict[seqSet[i][1]] = times
            num_mut.append(mut_count)

        origins = self._makeOrigins(mutSeqDict, num_mut)
        '''
        Estimate effective population size from data
        '''
        thetas = []

        for i in range(len(origins)):
            nu = origins[i]
            ns = num_mut[i]
            # Start at the left side of the interval
            x0 = 1
            # For the case nu==ns there is no solution -- outliers
            sol = optimize.fsolve(self._f, x0=x0, args=(nu,ns,))
            thetas.append(sol[0])


        return thetas, variance, num_seqs, num_mut, origins

    """
    Analyze bins - MLE origins
    """
    def _optimize(self, nu, ns):
        """
        Optimization function to call
        """
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

    def _fmle(self,x,nu,ns):
        '''
        Function of the estimate for root finding
        '''
#        print("Iteration: number of origins: ", nu)
#        print("           number of mutants: ", ns)
#        print("           theta:             ", x)

        a = float(x*np.log(1+ns/x))
        b = np.math.factorial(nu)
        c = np.exp(-(x*np.log(1+ns/x)))
        d = nu*math.log(a)
        dd = math.log(b)
        e = x*math.log(1+ns/x)
        return -(d-dd-e)

    def analyzeBinsMLE(self):
        # Number of (mutated) sequences
        num_mut = []
        num_seqs = []
        # Dictionary of mutants
        mutSeqDict = dict()

        # Length of evolution
        traj_length = len(self.dict_traj)

        # Filling the mutants dict and presence array
        variance_size = []
        for t, seqSet in enumerate(self.dict_traj):
            mut_count = 0
            if len(seqSet)!=0:
                print("Sample size: ",len(seqSet))
                num_seqs.append(len(seqSet))
                variance_size.append(1/len(seqSet))
                for i in range(len(seqSet)):
                    if (seqSet[i][1]!=[]):
                        mut_count += 1
                        if not seqSet[i][1] in mutSeqDict:
                            times = np.zeros(traj_length)
                            times[t] = 1
                            mutSeqDict[seqSet[i][1]] = times
                        else:
                            times = mutSeqDict[seqSet[i][1]]
                            if times[t]==1:
                                pass
                            elif times[t]==0:
                                times[t] = 1
                                mutSeqDict[seqSet[i][1]] = times
                num_mut.append(mut_count)
            else:
                num_seqs.append(0)
                variance_size.append(1)
                num_mut.append(mut_count)


        origins = self._makeOrigins(mutSeqDict, num_mut)


        '''
        Estimate effective population size from data
        '''
        thetas = []
        #Variance
        variance = []

        for i in range(len(origins)):
            sol = self._optimize(origins[i], num_mut[i])
            thetas.append(sol)

        for i in range(len(thetas)):
            if thetas[i]!=0:
                var = thetas[i]*math.log(1+num_mut[i]/thetas[i])
                variance.append(var)
            else:
                variance.append(0)

        return thetas, variance, variance_size, num_seqs, num_mut, origins


    """
    Analyze bins - MLE origins, not keeping track of seen origins
    """

    def analyzeBinsMLENT(self):
        # Number of (mutated) sequences
        num_mut = []
        num_seqs = []
        # Dictionary of mutants
        mutSeqDict = dict()

        # Length of evolution
        traj_length = len(self.dict_traj)

        # Filling the mutants dict and presence array
        variance_size = []
        for t, seqSet in enumerate(self.dict_traj):
            mut_count = 0
            if len(seqSet)!=0:
                num_seqs.append(len(seqSet))
                variance_size.append(1/len(seqSet))
                for i in range(len(seqSet)):
                    if (seqSet[i][1]!=[]):
                        mut_count += 1
                        if not seqSet[i][1] in mutSeqDict:
                            times = np.zeros(traj_length)
                            times[t] = 1
                            mutSeqDict[seqSet[i][1]] = times
                        else:
                            times = mutSeqDict[seqSet[i][1]]
                            if times[t]==1:
                                pass
                            elif times[t]==0:
                                times[t] = 1
                                mutSeqDict[seqSet[i][1]] = times
                num_mut.append(mut_count)
            else:
                num_seqs.append(0)
                variance_size.append(1)
                num_mut.append(mut_count)


        origins = self._makeOrigins(mutSeqDict, num_mut)


        '''
        Estimate effective population size from data
        '''
        thetas = []
        #Variance
        variance = []

        for i in range(len(origins)):
            sol = self._optimize(origins[i], num_mut[i])
            thetas.append(sol)

        for i in range(len(thetas)):
            if thetas[i]!=0:
                var = thetas[i]*math.log(1+num_mut[i]/thetas[i])
                variance.append(var)
            else:
                variance.append(0)

        return thetas, variance, variance_size, num_seqs, num_mut, origins



    """
    Analyze bins - from segregating sites
    """
    def _harm(self,n):
        an = 0
        for i in range(1,n-1):
            an += 1/i
        return an

    def analyzeBinsSegrS(self):
        thetas = []
        for t, seqSet in enumerate(self.dict_traj):
            sitesSet = []
            for seq, i in seqSet:
                for site in seq:
                    sitesSet.append(site)
            sitesSetDict = list(dict.fromkeys(sitesSet))
            K = len(sitesSetDict)
            an = self._harm(len(sitesSet))
            thetas.append(K/an)
        return thetas
