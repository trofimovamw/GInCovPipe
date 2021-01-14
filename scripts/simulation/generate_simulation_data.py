#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 18:426:09 2021

@authors: Maureen Smith, Maria Trofimova
"""

import sys
import argparse

import numpy as np

import os
import csv
import copy

import random
import string
import math

from pathlib import Path

from generate_evol_bins_gillespie import generateEvolGill
from simulate_sequence_evo import sequenceEvol

from modular_theta_est_mle import analyzeTrajectory

from datetime import datetime, timedelta, date



import matplotlib.pyplot as plt



# initialize argument parser
parser = argparse.ArgumentParser(description='Simulate sequence evolution.')

# Required arguments
parser.add_argument('-fasta', required=True,
                    help='fasta file to write the simulated sequences to')

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

parser.add_argument('-subsampling', type=float,
                    help='(optional) subsampling of all samples. if zero, take random subsample')

# Switch
parser.add_argument('--switch_orig', action='store_true', required=False,
					help='(optional) In the case of introductions and two replication rates: '
                         'If true: switch to secon replication rate at the same time point as the original outbreak. '
                         'If false: switch is at the half of the respective time frame.')

parser.add_argument('-init_seq', nargs='?',
                    help='(optional) Initial sequence. If not given, a random sequence is created.')


DAY0 = date.fromisoformat("2020-01-01")

args = parser.parse_args()

print("*"*100)
print("Running simulation of sequence evolution with argruments\n")
for a in vars(args):
    print(a, getattr(args, a))
    #print(' {} {}'.format(a, getattr(args, a) or ''))
print("*"*100)

evol_run = sequenceEvol(length=args.L,
                            p_repl=args.p_rep,
                            p_repl2=args.p_rep2,
                            p_death=args.p_death,
                            p_mut=args.p_mut,
                            N_init=args.N_init,
                            t_start=0,
                            t_final=args.t_final,
                            t_switch=math.floor(args.t_final / 2),
                            outFile=args.fasta)
# Sampling on tips only
#trajectory_indiv, time_trajectory, initSeq, tipsList = evol_run.mutate()
# Bin sampled tips
# Sizes of each sampled tip set
#sampling_sizes = [len(s) for s in tipsList]
