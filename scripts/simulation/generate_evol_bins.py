#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 11:58:05 2020

@authors: Maria Trofimova, Maureen Smith
"""


from random import choice, choices
import numpy as np
import csv
import os
import datetime
import random
import math

from pathlib import Path

strptime = datetime.datetime.strptime


class generateEvolPoiss:

    def __init__(self, length, p_repl, p_repl2, p_mut, N, t_start, t_final, t_switch, out, sim, total_sim, time_delta, init_seq):
        self.length = length
        self.p_repl = p_repl
        self.p_repl2 = p_repl2
        self.p_mut = p_mut
        self.N = N
        self.init = self._initSet(init_seq)
        self.t_start = t_start
        self.t_final = t_final
        self.t_switch = t_switch
        self.outdir = out
        self.sim = sim
        self.total_sim = total_sim
        self.time_delta = time_delta

    def _initSet(self,init_seq):
        '''
        Create initial sequence, if parameter is None, or use predefined initial sequence
        '''
        if init_seq=='none':
            alphabet = 'ACGT'
            sequences = []
            seq = ''.join(choice(alphabet) for i in range(self.length))
            sequences = [seq] * self.N
            return sequences
        else:
            init_seq_set = [init_seq] * self.N
            return init_seq_set

    def _mutateBase(self,base,initSeqBase):
        '''
        Mutate individual sequences according
        to Kimura model, fixed transition probability
        '''
        # Mutate if the base didnt mutate yet
        if base==initSeqBase:
            stay = 0
            alpha = 0.2
            beta = 0.2
            gamma = 0.6

            transitions = {'A':'T','T':'A','C':'G','G':'C'}
            transversions_alpha = {'A':'C','C':'A','T':'G','G':'T'}
            transversions_beta = {'A':'G','G':'A','T':'C','C':'T'}

            r = np.random.uniform(0,1)

            if stay <= r < (stay+gamma):
                return transitions[base]
            elif (stay+gamma) <= r < (stay+gamma+alpha):
                return transversions_alpha[base]
            elif (stay+gamma+alpha) <= r < (stay+gamma+alpha+beta):
                return transversions_beta[base]
            return base
        else:
            return base

    def _mutateBasenorm(self,base):
        '''
        Mutate individual sequences according
        to Kimura model, fixed transition probability
        '''
        stay = 0
        alpha = 0.2
        beta = 0.2
        gamma = 0.6

        transitions = {'A':'T','T':'A','C':'G','G':'C'}
        transversions_alpha = {'A':'C','C':'A','T':'G','G':'T'}
        transversions_beta = {'A':'G','G':'A','T':'C','C':'T'}

        r = np.random.uniform(0,1)

        if stay <= r < (stay+gamma):
            return transitions[base]
        elif (stay+gamma) <= r < (stay+gamma+alpha):
            return transversions_alpha[base]
        elif (stay+gamma+alpha) <= r < (stay+gamma+alpha+beta):
            return transversions_beta[base]
        return base

    def _writeFile(self, t, foldername, init, curr_time):
        '''
        Write the alignment at time t to file
        '''
        FILEPATH = Path(self.outdir)
        tname = str(t)
        tot_t_order = len(str(self.t_final))
        diff_t_tfinal = len(str(t))
        n_zero = tot_t_order-diff_t_tfinal
        zeros = "0" * n_zero
        tname = zeros+str(t)
        name = 't'+tname+'.fasta'

        file1 = FILEPATH / foldername
        file2 = file1 / 'bins'
        file = file2 / 'per_gen'
        file_name = file / name
        outfile = str(file_name)
        os.makedirs(file, exist_ok=True)

        with open(outfile, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            for i in range(len(init)):
                writer.writerow([">alignment_"+str(i)+'_'+str(curr_time.strftime("%Y-%m-%d"))])
                writer.writerow([init[i]])


    def mutate(self):
        '''
        Mutate initial sequence set over the course of t generations
        '''
        init = self.init
        initSeq = init[0]
        t = self.t_start
        time_steps = []
        time_steps.append(init)
        # Fake date to start with
        dt = "20/03/2015"
        curr_time = strptime(dt,"%d/%m/%Y").date()

        # Create folder for the simulation
        sim_name = self.sim
        tot_sim_order = len(str(self.total_sim))
        diff_sim_atot = len(str(self.sim))
        n_zero = tot_sim_order-diff_sim_atot
        zeros = "0" * n_zero
        sim_name = zeros+str(self.sim)
        foldername = 'timesteps_p_repl_' + str(self.p_repl) + '_p_mut_' + str(self.p_mut) + '_L_' + str(L) + '_N_' + str(self.N) + '_sim_' + sim_name

        trajectory = []

        curr_p_repl = self.p_repl

        while t<self.t_final:
            # Switch the replication rate somewhere in the middle if that is the correct mode (p_repl2!=0)
            if (self.p_repl2!=0) and (t>=self.t_switch):
                curr_p_repl = self.p_repl2

            new_set = []

            #num_repl = len(init)+1
            #while num_repl >= len(init):
            #    num_repl = np.random.poisson(curr_p_repl*len(init))

            # draw number of sequences of the next generation
            n_new = np.random.poisson(curr_p_repl*len(init))

            print("Sampled number of sequences in next generation: ",n_new)

            if n_new > 0:

                # draw sequences randomly with replacement
                new_set = random.choices(init, k=n_new)

                #for seq in repl_seqs:
                    #seq1 = seq
                    #seq2 = seq
                    #new_set.append(seq1)
                    #new_set.append(seq2)

                #print(len(new_set))

                num_mut_sites = len(new_set)*L + 1
                while num_mut_sites >= len(new_set)*L:
                    num_mut_sites = np.random.poisson(self.p_mut*len(new_set)*L)

                #print("Number of mutating sites: ",num_mut_sites)

                ind_mut_sites = np.arange(0,len(new_set)*L)

                #print(ind_mut_sites)

                # draw mutation sites without replacement
                mut_sites = random.sample(ind_mut_sites.tolist(),num_mut_sites)

                if num_mut_sites > 0:
                    for i in mut_sites:
                        seq_ind = math.floor(i/L)
                        seq_pos = i % L
                        seq = list(new_set[seq_ind])
                        #seq[seq_pos] = self._mutateBase(new_set[seq_ind][seq_pos],initSeq[seq_pos])
                        seq[seq_pos] = self._mutateBasenorm(new_set[seq_ind][seq_pos])
                        s = "".join(seq)
                        new_set[seq_ind] = s

                print("Number of sequences: ",len(new_set)," sequences")
                print(" ")
                print(" ")

            else:
                new_set = []

            '''
            Write file
            '''
            #self._writeFile(t, foldername, init, curr_time)
            trajectory.append(init)

            init = new_set
            t += 1
            #curr_time = curr_time+datetime.timedelta(days=self.time_delta)

            if new_set == []:
                break

        return trajectory, initSeq
