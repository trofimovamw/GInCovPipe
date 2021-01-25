#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 18:426:09 2021

@authors: Maureen Smith, Maria Trofimova
"""

import sys
import argparse
import os
import random
import math
import copy
import numpy as np

from sequence_evolution import sequenceEvol
from output_writer import writer



# initialize argument parser
parser = argparse.ArgumentParser(description='Simulate sequence evolution.')

# Required arguments
parser.add_argument('-file_prefix', required=True,
                    help='file_prefix for fasta file to write the simulated sequences to')

parser.add_argument('-output', required=True,
                    help='output path were the simulated fasta and table files are saved')

parser.add_argument('-L', type=int, required=True,
                    help='sequence length')

parser.add_argument('-N_init', type=int, required=True,
                    help='number of initial sequences')

parser.add_argument('-p_mut', type=float, required=True,
                    help='mutation rate')

parser.add_argument('-t_final', type=int, required=True,
                    help='number of time steps')

parser.add_argument('-p_death', type=float, nargs=1, required=False,
                    help='death rate')


# Optional arguments
parser.add_argument('-p_rep', type=float,
                    help='(optional) replication rate')

parser.add_argument('-p_rep2', type=float,
                    help='(optional) second replication rate')

# Switch
parser.add_argument('--switch_orig', action='store_true', required=False,
					help='(optional) In the case of introductions and two replication rates: '
                         'If true: switch to secon replication rate at the same time point as the original outbreak. '
                         'If false: switch is at the half of the respective time frame.')

parser.add_argument('-init_seq', nargs='+',
                    help='(optional) Initial sequence. If not given, a random sequence is created.')

parser.add_argument('-sub_rel', type=float, nargs='+', required=False,
                    help='(optional) list of relative subsampling of the complete sequence set')

parser.add_argument('-sub_abs', type=int, nargs='+', required=False,
                    help='(optional) list of absolute amount of subsampled sequences per day')

parser.add_argument('-intros', type=int, nargs='+', required=False,
                    help='(optional) list of number of introductions')


args = parser.parse_args()

print("*"*100)
print("Running simulation of sequence evolution with argruments\n")
for a in vars(args):
    print(a, getattr(args, a))
    #print(' {} {}'.format(a, getattr(args, a) or ''))
print("*"*100)


evol_sim = sequenceEvol(length=args.L,
                            p_repl=args.p_rep,
                            p_repl2=args.p_rep2,
                            p_death=args.p_death,
                            p_mut=args.p_mut,
                            N_init=args.N_init,
                            t_start=0,
                            t_final=args.t_final,
                            t_switch=math.floor(args.t_final / 2)
                        )

# run simulation
time_trajectory = evol_sim.evolve_poi()

wr = writer(outputpath=args.output, file_prefix=args.file_prefix)

# write initial sequence as reference
wr.write_reference(evol_sim.init_seq)

# always run with 0 introductions
num_intros = [0]
if args.intros is not None:
    num_intros = num_intros + args.intros

ts = range(args.t_final)

for ni in num_intros:
    # number of introduction per generation
    num_intros = np.zeros(args.t_final)
    num_intro_seq = np.zeros(args.t_final)
    file_suffix_intro = ""

    # needs to be copied, otherwise the original one is overwritten
    trajectory_withIntroduction = copy.deepcopy(time_trajectory)
    #print("Numer of introductions added: " + str(ni))
    if ni != 0:
        file_suffix_intro="_intros_"+str(ni)
        for t_start in random.choices(range(1, args.t_final - 1), k=ni):
            # random number of copies
            n_intro = random.randrange(11)
            evol_sim_intro = sequenceEvol(length=args.L,
                                    p_repl=args.p_rep,
                                    p_repl2=args.p_rep2,
                                    p_mut=args.p_mut,
                                    N_init=n_intro,
                                    t_start=t_start,
                                    t_final=args.t_final)

            trajectory_intro = evol_sim_intro.evolve_poi()
            for t in range(t_start, args.t_final):
                # if not empty
                if trajectory_intro[t-t_start]:
                    # add evolving introduction species to initial outbreak
                    for s, ns in trajectory_intro[t-t_start].items():
                        trajectory_withIntroduction[t][s] = trajectory_withIntroduction[t].get(s, 0) + ns
                        # number of new introduced sequences at time t
                        num_intro_seq[t] += ns

            # number of new introduced sequence type at time t
            num_intros[t_start] += 1

    # header=hCoV-19/Italy/LAZ-INMI1-isl/2020|EPI_ISL_410545|2020-01-29
    header_prefix = ">NS|"+file_suffix_intro
    file_suffix = "NS"+file_suffix_intro
    # write fasta with all sequences
    df_NS = wr.write_fasta(file_suffix=file_suffix,
                           header_prefix=header_prefix,
                           species_dict=trajectory_withIntroduction)

    df_NS["true_N"] = df_NS["sampled_N"]
    # add intro column
    df_NS["numIntros"] = num_intros
    df_NS["numIntrosSeqs"] = num_intro_seq

    wr.write_table(table=df_NS, file_suffix=file_suffix)
    wr.write_config_yaml(file_suffix=file_suffix)

    # write fasta with all absolute subsample
    if args.sub_abs is not None:
        for s_abs in args.sub_abs:
            print("---  Subsample sequence set taking " + str(s_abs) + " ---")
            header_prefix = header_prefix=">WS|" + str(s_abs) + "|" + file_suffix_intro
            file_suffix = "WSABS_"+ str(s_abs) + file_suffix_intro

            subsampled_time_trajectory = []
            for t in ts:
                num_seq = sum(trajectory_withIntroduction[t].values())
                time_trajectory_sub = {}
                # take abolsute subsample or N(t) if less
                if(num_seq > s_abs):
                    # sampling the particular sequences for each time point without replacemenet
                    seq_subset = random.sample(list(trajectory_withIntroduction[t].keys()), counts=trajectory_withIntroduction[t].values(), k=s_abs)
                    # counts sampled sequences
                    for s in seq_subset:
                        time_trajectory_sub[s] = time_trajectory_sub.get(s, 0) + 1
                else:
                    time_trajectory_sub = trajectory_withIntroduction[t]
                subsampled_time_trajectory.append(time_trajectory_sub)

            df_WS_abs = wr.write_fasta(file_suffix=file_suffix, header_prefix=header_prefix, species_dict=subsampled_time_trajectory)
            df_WS_abs["true_N"] = df_NS["true_N"]
            df_WS_abs["numIntros"] = num_intros
            df_WS_abs["numIntrosSeqs"] = num_intro_seq

            wr.write_table(table=df_WS_abs, file_suffix=file_suffix)
            wr.write_config_yaml(file_suffix=file_suffix)


    # write fasta with all relative subsample
    if args.sub_rel is not None:
        for s_rel in args.sub_rel:
            print("--- Subsample sequence set taking " + str(s_rel) + " ---")
            header_prefix = header_prefix=">WS|" + str(s_rel) + "|" + file_suffix_intro
            file_suffix = "WSREL_" + str(s_rel) + file_suffix_intro
            # total sample set size
            subsample_size = round(sum(df_NS["true_N"]) * s_rel)

            # sampling without replacement from which time point the sequences are coming (weighted by the number seqs)
            time_subset = random.sample(ts, k=subsample_size, counts=df_NS["true_N"])
            subsampled_time_trajectory = []
            for t in ts:
                #sampling the particular sequences for each time point without replacemenet
                seq_subset = random.sample(list(trajectory_withIntroduction[t].keys()), counts=trajectory_withIntroduction[t].values(), k=time_subset.count(t))

                # counts sampled sequences
                time_trajectory_sub = {}
                for s in seq_subset:
                    time_trajectory_sub[s] = time_trajectory_sub.get(s, 0) + 1
                subsampled_time_trajectory.append(time_trajectory_sub)

            df_WS_rel = wr.write_fasta(file_suffix=file_suffix, header_prefix=header_prefix, species_dict=subsampled_time_trajectory)
            df_WS_rel["true_N"] = df_NS["true_N"]
            df_WS_rel["numIntros"] = num_intros
            df_WS_rel["numIntrosSeqs"] = num_intro_seq

            wr.write_table(table=df_WS_rel, file_suffix=file_suffix)
            wr.write_config_yaml(file_suffix=file_suffix)
