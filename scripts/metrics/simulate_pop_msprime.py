#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 10:30:57 2020

@author: mariatrofimova
"""

import msprime
import random
import math
import numpy as np
import copy
from scipy.stats import gamma
from scipy.stats import beta
from scipy.stats import uniform
from scipy.stats import truncnorm
from scipy.stats import norm
from Bio import SeqIO



# Simulate populations using msprime package 
# Msprime assumes diploid case so Ne for comparison must be Ne/2 and for simulation Ne*2



class ABCclass:

    def __init__(self,dict_traj,thetas,mu,num_seqs,days_between_bins,initSeq):
        self.dict_traj = dict_traj
        self.thetas = thetas
        self.mu = mu
        self.num_seqs = num_seqs
        self.days_between_bins = days_between_bins
        self.initSeq = initSeq
        
    
    def metric1(self,seqSet):
        seqTypes = set(seqSet)
        sum_ = []
        for i in range(len(seqTypes)):
            for j in range(i+1,len(seqTypes)):
                a = seqTypes[i]
                b = seqTypes[j]
                diff = sum(1 for la, lb in zip(a, b) if la != lb)
                sum_.append(diff)
        return sum_/len(sum_)
    
    #Eg segregating sites
    def metric2(self,seqSet):
        seqTypes = set(seqSet)
        sum_ = []
        for i in range(len(seqTypes)):
            for j in range(i+1,len(seqTypes)):
                a = seqTypes[i]
                b = seqTypes[j]
                diff = [i for i in range(len(a)) if a[i] != b[i]]
                for ii in diff:
                    sum_.append(ii)
        return len(set(sum_))
    
    def _harmonic(self,n,power):
        an = 0
        for i in range(1,n):
            an += 1/(i**power)
        return an
    
    """
    Watterson, and Tajimas D
    """
    def analyzeBinsSegrS(self):
        # Count segregating sites - Watterson statistic
        segr_sites = []
        thetas_wat = []
        variance_wat = []
        # Count pairwise differences - Pi/Tajimas statistic
        dij_bins = []
        thetas_pi = []
        variance_pi = []
        tajimas_d = []
        # Count genotypes
        genotypes = []
        # Count number of mutations
        num_mut_sites = []
        for t, seqSet in enumerate(self.dict_traj):
            num_mut_sites_bin = []
            segr_sites_set = []
            gen = []
            for seq in seqSet:
                if seq[1] != '':
                    if not seq[1] in gen:
                        gen.append(seq[1])
                    split_str = seq[1].split('-')
                    for pos in split_str:
                        posit = pos.split('>')[0]
                        segr_sites_set.append(int(posit))
                    num_mut_sites_bin.append(len(split_str))
            num_mut_sites.append(num_mut_sites_bin)
            genotypes.append(len(gen))
            if len(seqSet)>1:
                dij = []
                for i in range(len(seqSet)-1):
                    for j in range(i+1,len(seqSet)):
                        # Set 1
                        seq1 = seqSet[i][1]
                        seq1set = seq1.split("-")
                        # Set 2
                        seq2 = seqSet[j][1]
                        seq2set = seq2.split("-")
                        difference = max(len(set(seq1set).difference(set(seq2set))),len(set(seq2set).difference(set(seq1set))))
                        dij.append(difference)
                theta_t = sum(dij)/(len(seqSet)*(len(seqSet)-1)/2)
                thetas_pi.append(theta_t)
                dij_bins.append(dij)
                segr_sites.append(len(list(set(segr_sites_set))))
                # a1
                a = self._harmonic(len(seqSet),1)
                # a2
                a_ = self._harmonic(len(seqSet),2)
                # b1
                b = (len(seqSet)+1)/(3*(len(seqSet)-1))
                # b2
                b_ = (2*(len(seqSet)**2+len(seqSet)+3))/(9*len(seqSet)*(len(seqSet)-1))
                # c1 
                c = b - 1/a
                # c2
                c_ = b_ - (len(seqSet)+2)/(a*len(seqSet)) + a_/(a**2)
                # e1 
                e = c/a
                # e2
                e_ = c_/(a**2+a_)
                # Variance
                S = len(list(set(segr_sites_set)))
                var = np.sqrt(e*S + e_*S*(S-1))
                theta_w = S/a
                thetas_wat.append(theta_w)
                D = (theta_t - theta_w)/var
                tajimas_d.append(D)
                # Variance
                variance_wat.append((a_*theta_w**2)/a**2)
                variance_pi.append((b_*theta_t**2))
            else:
                
                thetas_wat.append(0)
                dij_bins.append(0)
                thetas_pi.append(0)
                tajimas_d.append(0)
                # Variance
                variance_wat.append(0)
                variance_pi.append(0)
        
        return thetas_pi, thetas_wat, tajimas_d, variance_wat, variance_pi, segr_sites, genotypes
    
    
    # Sum of summary statistrics - by type or by bin?
    def distance(self,stat_trajectory):
        return 0        
    
#    def weightsNe(self,values,old_weights,Ne,tau):
#        weights = []
#        bins = np.linspace(min(values), max(values),50)
#        digitized = np.digitize(values, bins)
#        probs = []
#        sum_old_weights = []
#        for i,value in enumerate(values):
#            prob = gamma.pdf(value,Ne)
#            probs.append(prob)
#            old_weight = old_weights[i]
#            K = 
#            sum_old_weights.append(old_weight*K)
#        
#        
#    
    def calcHeterozygocityS(self,genotypes):
        ugeno = dict()
        for g in genotypes:
            print(g)
            if not g[1] in ugeno.keys():
                ugeno[g[1]] = 1
            else:
                ugeno[g[1]] += 1
        sumfreq = []
        for key,value in ugeno.items():
            sumfreq.append((value/len(genotypes))**2)
        H = 1 - sum(sumfreq)
        return H
    
    def calcHeterozygocityM(self,genotypes):
        ugeno = dict()
        for g_ in genotypes:
            gs = [str(i) for i in g_] 
            g = "".join(gs)
            if not g in ugeno.keys():
                ugeno[g] = 1
            else:
                ugeno[g] += 1
        sumfreq = []
        for key,value in ugeno.items():
            sumfreq.append((value/len(genotypes))**2)
        H = 1 - sum(sumfreq)
        return H
        
    def priorPDFmu(self,mu):
        if 0<mu<=1000: #Ner*0.1<Ne<=(Ner*100) 
            return 1
        else:
            return 0
    
    def priorPDFgr(self,gr):
        if 0<gr<=1000: #Ner*0.1<Ne<=(Ner*100) 
            return 1
        else:
            return 0
        
    def prior(self,Ne,mu,Ner):
        if Ner*0.1<Ne<=(Ner*5) and 0<mu<=0.1:
            return 1
        else:
            return 0
    
    def runSMC(self):        
        ref_fasta = SeqIO.parse(open(self.initSeq),'fasta')
        ref_seq = ''
        for fasta in ref_fasta:
            ref_seq = str(fasta.seq)
        lref = len(ref_seq)
        thetas = self.thetas
        print(thetas)
        mu = self.mu
        print(mu)
        thetas_pi, thetas_wat, tajimas_d, variance_wat, variance_pi, segr_sites, genotypes = self.analyzeBinsSegrS()
        
        #thetas_pi, thetas_wat, tajimas_d, variance_wat, variance_pi, segr_sites, genotypes, num_mut_seq = self.analyzeBinsSegrS(initSeq,trajectory)
        acceptedMusg = []
        for i in range(len(self.dict_traj)):
            acceptedMus = []
            # Statistic to compare with
            tajimaR = tajimas_d[i]
            # Individual bin intial guess
            mu_ = thetas[i]*mu
            Ne = 1#thetas[i]
            #while mu_==0:
            #    mu_ = random.uniform(0,0.1)
            #print("Scaled theta: ",mu_)
            N = len(self.dict_traj[i])
            L = lref
            
            H = self.calcHeterozygocityS(self.dict_traj[i])
            print("Heterozygocity ",H)
            print("Sample size N: ",N)
            print("Mutation rate: ",mu)
            print("Tajimas D: ",tajimaR)
            print("thetas_pi ",thetas_pi[i])
            print("Genotypes: ",genotypes[i])
            #genotypes = len(np.unique(np.array(trajectory[i])))
            print("Segregating sites: ",segr_sites[i])
            print("Length: ",L)
            print(len(thetas))
            print(len(self.dict_traj))
            
            eps1 = [80,30,10,4]
            
            # Run first sampling
            initMu = []
            #initGr = []
            while len(initMu)!=50:
                mup = norm.rvs(thetas[i]*mu, scale=2)
                while mup<=0:
                    mup = norm.rvs(thetas[i]*mu, scale=2)
#                Nep = norm.rvs(thetas[i], scale=3)
#                while Nep<=0:
#                    Nep = norm.rvs(thetas[i], scale=3)
#                sigma = 0.01
#                a,b = ((0 - mu_) / sigma), ((0.1 - mu_) / sigma)
#                mup = truncnorm.rvs(a, b, loc=mu_, scale=sigma, size=1)[0] #beta.rvs(a,b,size=1)[0]
#                while mup<=0:
#                   mup = truncnorm.rvs(a,b,mu_, scale=sigma)
                print(mup)
                evol_run = msprime.simulate(sample_size=N, Ne=1/2, 
                                            length=L, recombination_rate=0, 
                                            mutation_rate=mup,num_replicates = 1)
                #tmrcas = np.zeros(10)
                #tajimas = np.zeros(10)
                pis = np.zeros(1)
                segrs = np.zeros(1)
                gs = np.zeros(1)
                hs = np.zeros(1)
                for i, tree in enumerate(evol_run):
                    #tmrca = ftree.time(ftree.root)
                    #tmrcas[i] = tmrca
                    
                    #tajima = tree.Tajimas_D()
                    #tajimas[i] = tajima
                    
                    pi = tree.pairwise_diversity()
                    pis[i] = pi
                    
                    s = tree.num_sites
                    segrs[i] = s
                    
                    m = tree.genotype_matrix()
                    gt = np.transpose(m)
                    g_ = []
                    if gt!=[]:
                        g_ = np.unique(gt,axis=0)
                    gs[i] = len(g_)
                    h_ = self.calcHeterozygocityM(gt)
                    hs[i] = h_
                # posterior based on summary statistics
                #tmrca = np.mean(tmrcas)
                #tajima = np.mean(tajimas)
                #pi = np.mean(pis)
                #print(pi)
                #print("Mean pi ",pi)
                #s = np.mean(segrs)
                #print(s)
                g = np.mean(gs)
                hz = np.mean(hs)
                print(hz)
                #print("Heteroz. ", hz)
                #thetas_pi[-1],len(segr_sites[-1]),
                #pi,s
                r = np.array([genotypes[i],H])
                p = np.array([g,hz])
                l2 = np.linalg.norm(r-p)
                print(l2)
                if l2 < eps1[0]:
                    initMu.append(mup)
                    #initGr.append(Nep)
                    print("Num. accepted ",len(initMu))
                    print(mup)
                    mu_ = mup
                    #Ne = Nep
            t = 1
            # Assign weights
            weightsMu = np.full(len(initMu),1/len(initMu))
            #weightsGr = np.full(len(initGr),1/len(initGr))
            
            t_final = len(eps1)-1
            for t in range(1,t_final):
                print("SMC step: ",t)
                priorMu = initMu
                #priorGr = initGr
                stepMu = []
                #stepGr = []
                nweightsMu = []
                #nweightsGr = []
                while len(stepMu)<50:
                    #Sample from previous sample
                    from numpy.random import choice
                    mus = choice(priorMu, 1,
                                  p=weightsMu)[0]
#                    grs = choice(priorGr, 1,
#                                  p=weightsGr)[0]
#                    # Perturbation
#                    a1,b1 = ((0 - mus) / 0.1), ((1 - mus) / 0.1)
                    musp = norm.rvs(mus, scale=5)
#                    a2,b2 = ((0 - grs) / 1), ((3 - grs) / 1)
#                    grsp = norm.rvs(grs, scale=3)#truncnorm.rvs(a2,b2,grs,scale=1)
                    while musp<=0:
                        musp = norm.rvs(mus, scale=5)
                    # Simulate
                    evol_run = msprime.simulate(sample_size=N, Ne=1/2, 
                                            length=L, recombination_rate=0, 
                                            mutation_rate=musp,num_replicates = 1)
                    #tmrcas = np.zeros(10)
                    #tajimas = np.zeros(10)
                    pis = np.zeros(1)
                    segrs = np.zeros(1)
                    gs = np.zeros(1)
                    hs = np.zeros(1)
                    for i, tree in enumerate(evol_run):
                        #ftree = tree.first()
                        
                        pi = tree.pairwise_diversity()
                        pis[i] = pi
                        
                        s = tree.num_sites
                        segrs[i] = s
                        
                        m = tree.genotype_matrix()
                        gt = np.transpose(m)
                        g_ = []
                        if gt!=[]:
                            g_ = np.unique(gt,axis=0)
                        gs[i] = len(g_)
                        h_ = self.calcHeterozygocityM(gt)
                        hs[i] = h_
                    # posterior based on summary statistics
                    #tmrca = np.mean(tmrcas)
                    #tajima = np.mean(tajimas)
                    pi = np.mean(pis)
                    s = np.mean(segrs)
                    hz = np.mean(hs)
                    r = np.array([thetas_pi[-1],segr_sites[i],H])
                    p = np.array([pi,s,hz])
                    l2 = np.linalg.norm(r-p)
                    if l2 < eps1[t]:
                        print("Accepted ",len(stepMu))
                        stepMu.append(musp)
                        print(musp)
                        #print(grsp)
                        #stepGr.append(grsp)
                        kweightmu = []
                        #kweightgr = []
                        for j in range(len(weightsMu)):
                            wjm = weightsMu[j]
                            #wjg = weightsGr[j]
#                            a1,b1 = ((0 - mus) / 0.1), ((1 - mus) / 0.1)
#                            musp_ = truncnorm.pdf(musp,a1,b1,mus,scale=0.1)
#                            a2,b2 = ((0 - grs) / 1), ((3 - grs) / 1)
                            musp_ = norm.pdf(musp,mus, scale=5)#truncnorm.pdf(grsp,a2,b2,grs,scale=1)
                            while musp_ <= 0:
                                musp_ = norm.pdf(musp,mus, scale=5)
                            Kjm = wjm*musp_
                            #Kjg = wjg*grsp_
                            kweightmu.append(Kjm)
                            #kweightgr.append(Kjg)
                        wmu = self.priorPDFmu(musp_)/sum(kweightmu)
                        #wgr = self.priorPDFgr(grsp)/sum(kweightgr)
                        nweightsMu.append(wmu)
                        #nweightsGr.append(wgr)
                #t = t+1
                initMu = stepMu
                #initGr = stepGr
                nweightsMunorm = []
                #nweightsGrnorm = []
                # Normalize weights
                for el in nweightsMu:
                    nweightsMunorm.append(el/sum(nweightsMu))
#                for el in nweightsGr:
#                    nweightsGrnorm.append(el/sum(nweightsGr))
                
                weightsMu = nweightsMunorm
                #weightsGr = nweightsGrnorm
                acceptedMus.append(copy.copy(initMu))
                #acceptedGrs.append(copy.copy(initGr))
            acceptedMusg.append(acceptedMus[-1])
            #acceptedGrsg.append(acceptedGrs[-1])
                
        return acceptedMusg
            
    
    
    def runSMCs(self):        
        ref_fasta = SeqIO.parse(open(self.initSeq),'fasta')
        ref_seq = ''
        for fasta in ref_fasta:
            ref_seq = str(fasta.seq)
        lref = len(ref_seq)
        thetas = self.thetas
        print(thetas)
        mu = self.mu
        print(mu)
        thetas_pi, thetas_wat, tajimas_d, variance_wat, variance_pi, segr_sites, genotypes = self.analyzeBinsSegrS()
        # Fit parameters
        #thetas, mu = analyze1.calculateParamsFromLS(thetas2,mus,np.arange(1,len(thetas2)))
        # Last sample
        thetas_o = 1
        acceptedMus = []
        acceptedGrs = []
        for i in range(0,len(self.dict_traj)):
            # Statistic to compare with
            tajimaR = tajimas_d[-1]
            # Individual bin intial guess
            mu_ = random.uniform(0,0.1) #thetas[i]*mu
            gr_ = random.uniform(0,5)
            while mu_==0:
                mu_ = random.uniform(0,0.1)
            #print("Scaled theta: ",mu_)
            N = len(self.dict_traj[-1])
            L = lref
            
            #acceptedNe = []
            #acceptedMu = []
            H = self.calcHeterozygocityS(self.dict_traj[-1])
            print("Heterozygocity ",H)
            print("Sample size N: ",N)
            print("Mutation rate: ",mu)
            print("Tajimas D: ",tajimaR)
            print("thetas_pi ",thetas_pi[-1])
            print("Genotypes: ",genotypes[-1])
            #genotypes = len(np.unique(np.array(trajectory[i])))
            print("Segregating sites: ",segr_sites[-1])
            print("Length: ",L)
            
            eps1 = [20,17,14,4,0.1]
            
            # Run first sampling
            initMu = []
            initGr = []
            while len(initMu)!=50:
                mup = norm.rvs(thetas[i]*mu, scale=0.2)
                while mup<=0:
                    mup = norm.rvs(thetas[i]*mu, scale=0.2)
                # Simulate with replicates
                evol_run = msprime.simulate(
                        population_configurations=[
                                msprime.PopulationConfiguration(sample_size=N,
                                                                initial_size=1/2,
                                                                growth_rate=gr)],
                                                                length=L,
                                                                mutation_rate=mup,
                                                                num_replicates = 1)
                #tmrcas = np.zeros(1)
                #tajimas = np.zeros(1)
                pis = np.zeros(1)
                segrs = np.zeros(1)
                gs = np.zeros(1)
                hs = np.zeros(1)
                for i, tree in enumerate(evol_run):
                    ftree = tree.first()
                    tmrca = ftree.time(ftree.root)
                    tmrcas[i] = tmrca
                    
                    #tajima = tree.Tajimas_D()
                    #tajimas[i] = tajima
                    
                    pi = tree.pairwise_diversity()
                    pis[i] = pi
                    
                    s = tree.num_sites
                    segrs[i] = s
                    
                    m = tree.genotype_matrix()
                    gt = np.transpose(m)
                    g_ = []
                    if gt!=[]:
                        g_ = np.unique(gt,axis=0)
                    gs[i] = len(g_)
                    h_ = self.calcHeterozygocityM(gt)
                    hs[i] = h_
                # posterior based on summary statistics
                tmrca = np.mean(tmrcas)
                tajima = np.mean(tajimas)
                pi = np.mean(pis)
                #print("Mean pi ",pi)
                s = np.mean(segrs)
                #g = np.mean(gs)
                hz = np.mean(hs)
                #print("Heteroz. ", hz)
                #thetas_pi[-1],len(segr_sites[-1]),
                #pi,s
                r = np.array([thetas_pi[-1],len(segr_sites[i]),H])
                p = np.array([pi,s,hz])
                l2 = np.linalg.norm(r-p)
                if l2 < eps1[0]:
                    initMu.append(mup)
                    initGr.append(gr)
                    print("Num. accepted ",len(initMu))
                    mu_ = mup
                    gr_ = gr
            t = 1
            # Assign weights
            weightsMu = np.full(len(initMu),1/len(initMu))
            weightsGr = np.full(len(initGr),1/len(initGr))
            
            t_final = len(eps1)-1
            for t in range(1,t_final):
                print("SMC step: ",t)
                priorMu = initMu
                priorGr = initGr
                stepMu = []
                stepGr = []
                nweightsMu = []
                nweightsGr = []
                while len(stepMu)<50:
                    #Sample from previous sample
                    from numpy.random import choice
                    mus = choice(priorMu, 1,
                                  p=weightsMu)[0]
                    grs = choice(priorGr, 1,
                                  p=weightsGr)[0]
                    # Perturbation
                    a1,b1 = ((0 - mus) / 0.1), ((1 - mus) / 0.1)
                    musp = truncnorm.rvs(a1,b1,mus,scale=0.1)
                    a2,b2 = ((0 - grs) / 1), ((10 - grs) / 1)
                    grsp = truncnorm.rvs(a2,b2,grs,scale=1)
                    # Simulate
                    evol_run = msprime.simulate(
                        population_configurations=[msprime.PopulationConfiguration(sample_size=N,
                                                                initial_size=1/2,
                                                                growth_rate=grsp)],
                                                                length=L,
                                                                mutation_rate=musp,
                                                                num_replicates = 10)
                    #tmrcas = np.zeros(10)
                    #tajimas = np.zeros(10)
                    pis = np.zeros(10)
                    segrs = np.zeros(10)
                    gs = np.zeros(10)
                    hs = np.zeros(10)
                    for i, tree in enumerate(evol_run):
                        ftree = tree.first()
                        #tmrca = ftree.time(ftree.root)
                        #tmrcas[i] = tmrca
                        
                        #tajima = tree.Tajimas_D()
                        #tajimas[i] = tajima
                        
                        pi = tree.pairwise_diversity()
                        pis[i] = pi
                        
                        s = tree.num_sites
                        segrs[i] = s
                        
                        m = tree.genotype_matrix()
                        gt = np.transpose(m)
                        g_ = []
                        if gt!=[]:
                            g_ = np.unique(gt,axis=0)
                        gs[i] = len(g_)
                        h_ = self.calcHeterozygocityM(gt)
                        hs[i] = h_
                    # posterior based on summary statistics
                    #tmrca = np.mean(tmrcas)
                    #tajima = np.mean(tajimas)
                    pi = np.mean(pis)
                    s = np.mean(segrs)
                    hz = np.mean(hs)
                    r = np.array([thetas_pi[-1],len(segr_sites[i]),H])
                    p = np.array([pi,s,hz])
                    l2 = np.linalg.norm(r-p)
                    if l2 < eps1[t]:
                        print("Accepted ",len(stepMu))
                        stepMu.append(musp)
                        print(musp)
                        print(grsp)
                        stepGr.append(grsp)
                        kweightmu = []
                        kweightgr = []
                        for j in range(len(weightsMu)):
                            wjm = weightsMu[j]
                            wjg = weightsGr[j]
                            a1,b1 = ((0 - mus) / 0.1), ((1 - mus) / 0.1)
                            musp_ = truncnorm.pdf(musp,a1,b1,mus,scale=0.1)
                            a2,b2 = ((0 - grs) / 1), ((10 - grs) / 1)
                            grsp_ = truncnorm.pdf(grsp,a2,b2,grs,scale=1)
                            Kjm = wjm*musp_
                            Kjg = wjg*grsp_
                            kweightmu.append(Kjm)
                            kweightgr.append(Kjg)
                        wmu = self.priorPDFmu(musp)/sum(kweightmu)
                        wgr = self.priorPDFgr(grsp)/sum(kweightgr)
                        nweightsMu.append(wmu)
                        nweightsGr.append(wgr)
                #t = t+1
                initMu = stepMu
                initGr = stepGr
                nweightsMunorm = []
                nweightsGrnorm = []
                # Normalize weights
                for el in nweightsMu:
                    nweightsMunorm.append(el/sum(nweightsMu))
                for el in nweightsGr:
                    nweightsGrnorm.append(el/sum(nweightsGr))
                
                weightsMu = nweightsMunorm
                weightsGr = nweightsGrnorm
                acceptedMus.append(copy.copy(initMu))
                acceptedGrs.append(copy.copy(initGr))
        return acceptedMus[-1], acceptedGrs[-1]
            
                

        
    def runMCMCavrg(self):
        # Initial guess - LS
        thetas = self.thetas[0]
        print(thetas)
        mu = self.mu
        num_seqs = self.num_seqs
        print(mu)
        thetas_pi, thetas_wat, tajimas_d, variance_wat, variance_pi, segr_sites, genotypes = self.analyzeBinsSegrS()
        thetas_pim = sum(thetas_pi)/len(thetas_pi)
        #thetas_watm = sum(thetas_wat)/len(thetas_wat)
        segr_sitesm = math.floor(sum(segr_sites)/len(segr_sites))
        genotypesm = math.floor(sum(genotypes)/len(genotypes))
        acceptedNes = []
        acceptedMus = []
        print(num_seqs)
        # Statistic to compare with
        #tajimaR = tajimas_d[i]
        # Individual bin intial guess
        print(mu)
        print(thetas)
        mu_ = mu#*thetas
        N = math.floor(sum(num_seqs)/len(num_seqs))
        Ne = thetas
        L = 29903
        eps1 = 20
        acceptedNe = []
        acceptedMu = []
        print("Sample size N: ",N)
        print("Ne: ",Ne)
        print("Mutation rate: ",mu)
        print("Tajimas pi: ",thetas_pim)
        print("Genotypes ",genotypesm)
        print("Segr. sites: ",segr_sitesm)
        t = 0
        # Initial sampling from prior
        for i in range(len(thetas_o)):
            # Statistic to compare with
            tajimaR = tajimas_d[-1]
            # Individual bin intial guess
            mu_ = random.uniform(0,0.1) #thetas[i]*mu
            gr_ = 1
            while mu_==0:
                mu_ = random.uniform(0,0.1)
            #print("Scaled theta: ",mu_)
            N = len(self.dict_traj[-1])#num_seqs[i]
            Ne = thetas_o #thetas_o[i]/mu
            L = len(initSeq)
            eps1 = 30
            acceptedNe = []
            acceptedMu = []
            print("Sample size N: ",N)
            print("Ne: ",Ne)
            print("Mutation rate: ",mu)
            print("Tajimas D: ",tajimaR)
            print("thetas_pi ",thetas_pi[-1])
            print("Genotypes: ",genotypes[-1])
            #genotypes = len(np.unique(np.array(trajectory[i])))
            print("Segregating sites: ",segr_sites[-1])
            print("Length: ",L)
            print("Esti time to root: ",10)
            tmrcar = len(days_between_bins)*1.6
            
            t = 0
            # Initial sampling from prior
            while len(acceptedMu)<500:
                #Nep = gamma.rvs(Ne, scale=5)
                #while Nep<=0:
                #    Nep = gamma.rvs(Ne, scale=5)
                sigma = 0.1
                a,b = ((0 - mu_) / sigma), ((1 - mu_) / sigma)
                mup = truncnorm.rvs(a, b, loc=mu_, scale=sigma, size=1)[0] #beta.rvs(a,b,size=1)[0]
                while mup<=0:
                   mup = truncnorm.rvs(a,b,mu_, scale=sigma)
                sigma = 2
                a,b = ((0 - gr_) / sigma), ((1 - gr_) / sigma)
                gr = truncnorm.rvs(a, b, loc=gr_, scale=sigma, size=1)[0]
                while gr==0:
                    gr = truncnorm.rvs(a, b, loc=gr_, scale=sigma, size=1)[0]
                print(mup)
                #print(Nep)
#                sigma = 0.1
#                a,b = ((0 - mu_) / sigma), ((1 - mu_) / sigma)
#                mup = truncnorm.rvs(a, b, loc=mu_, scale=sigma, size=1)[0] #beta.rvs(a,b,size=1)[0]
                evol_run = msprime.simulate(
                        population_configurations=[
                                msprime.PopulationConfiguration(sample_size=N,
                                                                initial_size=1/2,
                                                                growth_rate=gr)],
                                                                length=L,
                                                                mutation_rate=mup,
                                                                num_replicates = 1)
                tmrcas = np.zeros(1)
                tajimas = np.zeros(1)
                pis = np.zeros(1)
                segrs = np.zeros(1)
                gs = np.zeros(1)
                for i, tree in enumerate(evol_run):
                    ftree = tree.first()
                    tmrca = ftree.time(ftree.root)
                    tmrcas[i] = tmrca
                    
                    tajima = tree.Tajimas_D()
                    tajimas[i] =tajima
                    
                    pi = tree.pairwise_diversity()
                    pis[i] = pi
                    
                    s = tree.num_sites
                    segrs[i] = s
                    
                    m = tree.genotype_matrix()
                    gt = np.transpose(m)
                    g_ = []
                    if gt!=[]:
                        g_ = np.unique(gt,axis=0)
                    gs[i] = len(g_)
                # posterior based on summary statistics
                tmrca = np.mean(tmrcas)
                tajima = np.mean(tajimas)
                pi = np.mean(pis)
                s = np.mean(segrs)
                g = np.mean(gs)
 
                r = np.array([thetas_pi[-1],len(segr_sites[-1]),genotypes[-1],tajimaR])
                p = np.array([pi,s,g,tajima])
                l2 = np.linalg.norm(r-p)
                
                print("Distances: ",l2)
                if l2 < eps1:
                    print("Distances: ",l2)
                    #upp1 = gamma.pdf(Nep,Ne, scale=5)
                    a1,b1 = ((0 - mu_) / sigma), ((1 - mu_) / sigma)
                    upp2 = truncnorm.pdf(mup,a1,b1,loc=mu_, scale=sigma) #beta.pdf(mup,a,b)
                    #down1 = gamma.pdf(Ne,Nep, scale=5)
                    a2,b2 = ((0 - mup) / sigma), ((1 - mup) / sigma)
                    down2 = truncnorm.pdf(mu_,a2,b2,loc=mup, scale=sigma) #beta.pdf(mu_,a,b)
                    priorup = self.prior(mu)
                    print(priorup)
                    acc_ratio = min(1,(upp2*priorup)/(down2))
                    u = random.uniform(0,1)
                    if u < acc_ratio:
                        print("    Distance")
                        print(l2)
                        print("Segr. ",s)
                        print("Geno. ",g)
                        print("    Accepted")
                        print(mup)
                        acceptedMu.append(mup)
                        acceptedNe.append(gr)
                        print("Num. accepted ",len(acceptedMu))
                        mu_ = mup
                        gr_ = gr
                        #Ne = Nep
                t = t+1
                if t >= 100000:
                    print("Could not converge in 100000 steps...")
                    break
#            tau = 
#            t = 1
#            acceptedNei = []
#            acceptedMui = []
#            eps = [9,8,7,6,5,4,3,2]
#            for e in eps:
#                pickNe = random.choices(acceptedNe, weights=np.full(len(init), 1/len(init)),k=3)
            print("Done in ",t," steps")
            #acceptedNes.append(acceptedNe)
            acceptedMus.append(acceptedMu)
        return acceptedMus#, thetas2, mu, trajectory, tajimas_d, thetas_pi, segr_sites, num_seqs, num_mut_seq
                
        
    def runMCMC(self):
        # Initial guess - LS
        thetas = self.thetas
        print(thetas)
        mu = self.mu
        num_seqs = self.num_seqs
        print(mu)
        thetas_pi, thetas_wat, tajimas_d, variance_wat, variance_pi, segr_sites, genotypes = self.analyzeBinsSegrS()
        acceptedNes = []
        acceptedMus = []
        print(num_seqs)
        for i in range(len(tajimas_d)):
            if num_seqs[i]>5:
                # Statistic to compare with
                #tajimaR = tajimas_d[i]
                # Individual bin intial guess
                mu_ = mu*thetas[i]
                N = num_seqs[i]
                Ne = thetas[i]
                L = 29992
                eps1 = 20
                acceptedNe = []
                acceptedMu = []
                print("Sample size N: ",N)
                print("Ne: ",Ne)
                print("Mutation rate: ",mu)
                #print("Tajimas D: ",tajimaR)
                print("Genotypes ",len(genotypes[i]))
                t = 0
                # Initial sampling from prior
                while len(acceptedNe)<500:
                    Nep = norm.rvs(Ne, scale=3)
                    while Nep<=0:
                        Nep = norm.rvs(Ne, scale=3)
                    sigma = 0.01
                    a,b = ((0 - mu_) / sigma), ((0.1 - mu_) / sigma)
                    mup = truncnorm.rvs(a, b, loc=mu_, scale=sigma, size=1)[0] #beta.rvs(a,b,size=1)[0]
                    while mup<=0:
                       mup = truncnorm.rvs(a,b,mu_, scale=sigma)
                    evol_run = msprime.simulate(sample_size=N, Ne=Nep, 
                                                length=L, recombination_rate=0, 
                                                mutation_rate=mup)
                    # posterior based on summary statistics
                    #tajima = evol_run.Tajimas_D()
                    pi = evol_run.pairwise_diversity()
                    s = evol_run.num_sites
                    m = evol_run.genotype_matrix()
                    gt = np.transpose(m)
                    g = []
                    if gt!=[]:
                        g = np.unique(gt,axis=0)
                    #print(pi)
                    #print(thetas_pi[i])
    #                d1 = 0
    #                if tajima>=0 and tajimaR>=0:
    #                    d1 = abs(tajima-tajimaR)
    #                elif tajima>=0 and tajimaR<0:
    #                    d1 = tajima + abs(tajimaR)
    #                else:
    #                    d1 = abs(abs(tajima)-abs(tajimaR))
#                    d2 = abs(abs(pi)-abs(thetas_pi[i]))
#                    d3 = abs(len(segr_sites[i])-s)
#                    d4 = abs(len(genotypes[i])-len(g))
    #                print("Segregating sites theor.: ",s)
    #                print("Segregating sites real.: ",len(segr_sites[i]))
    #                d4 = abs(genotypes-len(g))
#                    d = np.array([d2,d3,d4])
#                    dt = np.transpose(d)
#                    am = np.array([[0.4,0,0],[0,0.3,0],[0,0,0.3]])
#                    l1 = np.dot(d,am)
#                    l2 = np.dot(l1,dt)
                    r = np.array([thetas_pi[i],len(segr_sites[i]),len(genotypes[i])])
                    p = np.array([pi,s,len(g)])
                    l2 = np.linalg.norm(r-p)
                    if l2 < eps1:
                        upp1 = norm.pdf(Nep,Ne, scale=3)
                        a1,b1 = ((0 - mu_) / sigma), ((0.1 - mu_) / sigma)
                        upp2 = truncnorm.pdf(mup,a1,b1,loc=mu_, scale=sigma) #beta.pdf(mup,a,b)
                        down1 = norm.pdf(Ne,Nep, scale=3)
                        a2,b2 = ((0 - mup) / sigma), ((0.1 - mup) / sigma)
                        down2 = truncnorm.pdf(mu_,a2,b2,loc=mup, scale=sigma) #beta.pdf(mu_,a,b)
                        priorup = self.prior(Nep,mu,thetas[i])
                        
                        acc_ratio = min(1,(upp1*upp2*priorup)/(down1*down2*1))
                        u = random.uniform(0,1)
                        if u < acc_ratio:
                            print("    Distance")
                            print(l2)
                            print("Segr. ",s)
                            print("Geno. ",len(g))
                            print("    Accepted")
                            print(Nep)
                            print(mup)
                            acceptedNe.append(Nep)
                            acceptedMu.append(mup)
                            print("Num. accepted ",len(acceptedNe))
                            mu_ = mup
                            Ne = Nep
                    t = t+1
                    if t >= 100000:
                        print("Could not converge in 100000 steps...")
                        break

                print("Done in ",t," steps")
                acceptedNes.append(acceptedNe)
                acceptedMus.append(acceptedMu)
            else:
                acceptedNes.append([])
                acceptedMus.append([])
        return acceptedNes, acceptedMus
  