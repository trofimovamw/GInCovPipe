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
from sequence_evolution import sequenceEvol

from modular_theta_est_mle import analyzeTrajectory

from datetime import datetime, timedelta, date



import matplotlib.pyplot as plt



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

parser.add_argument('-subsampling', type=float,
                    help='(optional) subsampling of all samples. if zero, take random subsample')

# Switch
parser.add_argument('--switch_orig', action='store_true', required=False,
					help='(optional) In the case of introductions and two replication rates: '
                         'If true: switch to secon replication rate at the same time point as the original outbreak. '
                         'If false: switch is at the half of the respective time frame.')

parser.add_argument('-init_seq', nargs='?',
                    help='(optional) Initial sequence. If not given, a random sequence is created.')


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
                            file_prefix=args.file_prefix)

# run simulation
time_trajectory = evol_run._evolve_poi()

# header=hCoV-19/Italy/LAZ-INMI1-isl/2020|EPI_ISL_410545|2020-01-29
header_prefix=">NS|"
file_suffix = "_NS"

file_str = args.output + "/" + args.file_prefix + file_suffix +".fasta"

# write fasta with all sequences
evol_run._write_fasta(outputfile=file_str, header_prefix=header_prefix, species_dict=time_trajectory)

# create files for subsampling
#for s_abs in


# Bin sampled tips
# Sizes of each sampled tip set
#sampling_sizes = [len(s) for s in tipsList]
