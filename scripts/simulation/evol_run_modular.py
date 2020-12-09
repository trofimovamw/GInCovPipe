#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 12:28:07 2020

@authors: Maureen Smith, Maria Trofimova
"""

import sys

import numpy as np

import os
import csv
import copy

import random
import string
import math

from pathlib import Path

from generate_evol_bins_gillespie import generateEvolGill

from modular_theta_est_mle import analyzeTrajectory

from datetime import datetime, timedelta, date

import matplotlib.pyplot as plt


FILEPATH = Path(__file__).parent

LEN_SEQ = 300 
NUM_SEQ_INIT = 101
P_MUT = 0.0001
P_REPL=1.15
T_FINAL = 3
SIM_ORDER = 20
TOTAL_SIM = 1
SAMPLING_RATE = 0.2
# all bin sizes 2 to N_BIN (by 2) are evaluated
N_BIN = 2
# zero, if there should be no change of growth rate, otherwise something like 0.49
P_REPL2 = 0.9
#number of introduction of random sequence (same start replication rate as original run)
NUM_INTRO = round(1/T_FINAL)
# replication rate switch for the introduced sequence
INTRO_REPL2 = 0.9
# if true: switch to repl2 at the same time point as the original outbreak, false: switch is at the half  the respective time frame
T_SWITCH_ORIG = True
# subsampling of all samples. if zero, take random subsample
SUBSAMPLE = 0.0

INIT_SEQ = "ACGACACGACACCAATTTAGATTTAGACCATTAGCACAAGACCGTAAAACGTCGCT"

DAY0 = date.fromisoformat("2020-01-01")

# casting a string to boolean
def str_to_bool(s):
    if s.lower() == 'true':
         return True
    else :
         return False

print("Number of arguments: " + str(len(sys.argv)))

if(len(sys.argv) > 1):
    FILEPATH = Path(sys.argv[1])
    print("OUTDIR "+ sys.argv[1])
if(len(sys.argv) > 2):
    P_REPL = float(sys.argv[2])
    print("P_REPL " + sys.argv[2])
if(len(sys.argv) > 3):
    P_MUT = float(sys.argv[3])
    print("P_MUT " + sys.argv[3])
if(len(sys.argv) > 4):
    LEN_SEQ = int(sys.argv[4])
    print("LEN_SEQ " + sys.argv[4])
if (len(sys.argv) > 5):
    NUM_SEQ_INIT = int(sys.argv[5])
    print("NUM_SEQ_INIT " + sys.argv[5])
if (len(sys.argv) > 6):
    T_FINAL = int(sys.argv[6])
    print("T_FINAL " + sys.argv[6])
if (len(sys.argv) > 7):
    SIM_ORDER = int(sys.argv[7])
    print("SIM_ORDER " + sys.argv[7])
if (len(sys.argv) > 8):
    TOTAL_SIM = int(sys.argv[8])
    print("TOTAL_SIM " + sys.argv[8])
if (len(sys.argv) > 9):
    N_BIN = int(sys.argv[9])
    print("N_BIN " + sys.argv[9])
if (len(sys.argv) > 10):
    P_REPL2 = float(sys.argv[10])
    print("P_REPL2 " + sys.argv[10])
if (len(sys.argv) > 11):
    NUM_INTRO = int(sys.argv[11])
    print("NUM_INTRO " + sys.argv[11])
if (len(sys.argv) > 12):
    INTRO_REPL2 = float(sys.argv[12])
    print("INTRO_REPL2 " + sys.argv[12])
if (len(sys.argv) > 13):
        T_SWITCH_ORIG = str_to_bool(sys.argv[13])
        print("T_SWITCH_ORIG " + sys.argv[13])
if (len(sys.argv) > 14):
        SUBSAMPLE = float(sys.argv[14])
        print("SUBSAMPLE " + sys.argv[14])

if (len(sys.argv) > 15):
    INIT_SEQ = str(sys.argv[15])
    print("INIT_SEQ " + sys.argv[15])


if (len(sys.argv) <= 16):
    
    sim_folder = "simulation_p_repl_"+str(P_REPL)+"_p_mut_"+str(P_MUT)

    OUTDIR = FILEPATH / sim_folder
    
    time_delta = 7
    p_death = P_REPL-0.2
    # If no change of replication rate is wanted, set p_repl2 to 0
    # If no particular initial sequence is wanted, set init_seq to 'None'
    evol_run = generateEvolGill(length=LEN_SEQ,
                                 p_repl=P_REPL,
                                 p_repl2=P_REPL2,
                                 p_mut=P_MUT,
                                 p_death = p_death,
                                 N=NUM_SEQ_INIT,
                                 t_start=0,
                                 t_final=T_FINAL,
                                 t_switch=math.floor(T_FINAL/2),
                                 out=OUTDIR,
                                 sim=SIM_ORDER,
                                 total_sim=TOTAL_SIM,
                                 #time_delta=time_delta,
                                 init_seq='none')
    # Sampling on tips only
    trajectory_indiv, time_trajectory, initSeq, tipsList = evol_run.mutate()
    # Bin sampled tips
    # Sizes of each sampled tip set
    sampling_sizes = [len(s) for s in tipsList]
    # Prefered binning
    NUM_BINS = 20
    binsize = math.ceil(sum(sampling_sizes)/NUM_BINS)
    # Make bins of equal size
    bins = [[] for i in range(NUM_BINS)]
    sampling_merged = [j for i in tipsList for j in i]
    for i in range(len(sampling_merged)):
        bins[math.floor(i / binsize)].append(sampling_merged[i])
    # Binned trajectory
    trajectory = bins
    timestepd = OUTDIR.parent

    p_out = "python_out_p_repl_" + str(P_REPL) + "_repl2_" + str(P_REPL2) + "_p_mut_" + str(
        P_MUT) + "_switchInit_" + str(T_SWITCH_ORIG)
    print("Writing results to "+ p_out)

    t = timestepd / p_out

    os.makedirs(str(t), exist_ok=True)


    # if introductions are added only compute NS and WSMB with intros, otherwise compute
    if NUM_INTRO > 0 :

        '''
        Add introduced sequences
        '''
        # number of introduction per generation
        numIntros = np.zeros(len(trajectory))

        trajectory_withIntroduction = copy.deepcopy(trajectory)
        print("TOTAL INTROS ", NUM_INTRO)
        for t_start in  np.random.randint(1, T_FINAL-1, size=NUM_INTRO):
            #TODO: for now random, but maybe decide on similarity to existing sequences
            intro_sequence = ''.join(random.choices("ACGT", k=LEN_SEQ))
            #t_start = round((T_FINAL*ni) /(NUM_INTRO+1))
            t_switch = math.floor((T_FINAL - t_start)/2)
            #random number of copies
            n_intro = np.random.randint(1,11)
            #n_intro=50
            if T_SWITCH_ORIG:
                t_switch = evol_run.t_switch
            
            evol_run = generateEvolGill(length=LEN_SEQ,
                                 p_repl=P_REPL,
                                 p_repl2=P_REPL2,
                                 p_mut=P_MUT,
                                 p_death = p_death,
                                 N=NUM_SEQ_INIT,
                                 t_start=0,
                                 t_final=T_FINAL,
                                 t_switch=math.floor(T_FINAL/2),
                                 out=OUTDIR,
                                 sim=SIM_ORDER,
                                 total_sim=TOTAL_SIM,
                                 #time_delta=time_delta,
                                 init_seq=intro_sequence)
            trajectory_intro, time_trajectory_intro, initSeq_intro, tipsList_intro  = evol_intro.mutate()
            # Bin sampled tips
            # Sizes of each sampled tip set
            sampling_sizes = [len(s) for s in tipsList_intro]
            # Prefered binning
            NUM_BINS = 20
            binsize = math.ceil(sum(sampling_sizes)/NUM_BINS)
            # Make bins of equal size
            bins = [[] for i in range(NUM_BINS)]
            sampling_merged = [j for i in tipsList for j in i]
            for i in range(len(sampling_merged)):
                bins[math.floor(i / binsize)].append(sampling_merged[i])
            # Binned trajectory
            numIntros[t_start] += 1
            trajectory_intro = bins
            print(" t ", t_start)
            for i in range(len(trajectory_intro)) :
                trajectory_withIntroduction[t_start+i].extend(trajectory_intro[i])
                print("time ", t_start+i, " seq ", len(trajectory_intro[i]))

        analyze_run_intro = analyzeTrajectory(trajectory_withIntroduction, initSeq)

        ### Output with introduction

        """
        Analyze without subsampling
        """
        days_since1_intro, thetas1_intro, variance1_intro, num_seqs1_intro, origins1_intro, wattersonEst, avrgSegrSites = analyze_run_intro.analyzeBinsNS()
        n1 = "NS_theta_origins_sim_" + str(SIM_ORDER) + "_intro_" + str(NUM_INTRO) + ".tsv"
        name = t / n1
        with open(str(name), 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(
                ["t", "value", "variance", "origins", "trueN", "sampled_N", "bin_size", "meanBinDate", "numIntros"])
            for i, d in enumerate(days_since1_intro):
                date = (DAY0 + timedelta(days=d)).strftime("%Y-%m-%d")
                writer.writerow(
                    [days_since1_intro[i], thetas1_intro[i], variance1_intro[i], origins1_intro[i], num_seqs1_intro[i],
                     num_seqs1_intro[i], 1, date, numIntros[i]])

        """
        Analyze with bin merging and subsampling
        """

        for nb in range(2, N_BIN + 1, 2):

            days_since4_intro, thetas4_intro, variance4_intro, num_seqs4_intro, bin_sizes4_intro, subsample_sizes4_intro, origins4_intro \
                = analyze_run_intro.analyzeBinsWSMB(nb)

            # get number of introductions per bin
            numIntros_bin = np.zeros(len(days_since4_intro))
            idx = 0
            nb_flag = bin_sizes4_intro[idx]
            for numIn in numIntros:
                if nb_flag == 0:
                    idx += 1
                    nb_flag = bin_sizes4_intro[idx]
                numIntros_bin[idx] += numIn
                nb_flag -= 1

            n4_intro = "WSMB_theta_origins_sim_" + str(SIM_ORDER) + "_numBin_" + str(nb) + "_intro_" + str(
                NUM_INTRO) + ".tsv"
            name = t / n4_intro
            with open(str(name), 'w+', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(
                    ["t", "value", "variance", "origins", "trueN", "sampled_N", "bin_size", "meanBinDate", "numIntros"])
                for i, d in enumerate(days_since4_intro):
                    date = (DAY0 + timedelta(days=d)).strftime("%Y-%m-%d")
                    writer.writerow(
                        [days_since4_intro[i], thetas4_intro[i], variance4_intro[i], origins4_intro[i],
                         num_seqs4_intro[i],
                         subsample_sizes4_intro[i], bin_sizes4_intro[i], date, numIntros_bin[i]])

    else :

        analyze_run = analyzeTrajectory(trajectory, initSeq)

        """
        Analyze without subsampling
        """
        print("Analyze without subsampling")
        days_since1, thetas1, variance1, num_seqs1, origins1, wattersonEst1, avrgSegrSites1 = analyze_run.analyzeBinsNS()

        n1 = "NS_theta_origins_sim_"+str(SIM_ORDER)+".tsv"
        name = t / n1
        with open(str(name), 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["t","value","variance","origins","trueN","sampled_N","bin_size",  "meanBinDate"])
            for i, d in enumerate(days_since1):
                date=(DAY0 + timedelta(days=d)).strftime("%Y-%m-%d")
                writer.writerow([days_since1[i],thetas1[i],variance1[i],origins1[i],num_seqs1[i],num_seqs1[i],1,date])


        """
        Analyze with subsampling
        """
        print("Analyze with subsampling")
        days_since2, thetas2, variance2, num_seqs2, subsample_sizes2, origins2 = analyze_run.analyzeBinsWS()

        n2 = "WS_theta_origins_sim_"+str(SIM_ORDER)+".tsv"
        name = t / n2
        with open(str(name), 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["t","value","variance","origins","trueN","sampled_N", "bin_size", "meanBinDate"])
            for i,d in enumerate(days_since2):
                date = (DAY0 + timedelta(days=d)).strftime("%Y-%m-%d")
                writer.writerow([days_since2[i],thetas2[i],variance2[i],origins2[i],num_seqs2[i],subsample_sizes2[i],1,date])



        for nb in range(2, N_BIN+1, 2):

            """
            Analyze with bin merging
            """
            print("Analyze with bin merging")
            days_since3, thetas3, variance3, num_seqs3, bin_sizes3, origins3 = analyze_run.analyzeBinsMB(nb)

            n3 = "MB_theta_origins_sim_" + str(SIM_ORDER) + "_numBin_" + str(nb) + ".tsv"
            name = t / n3
            with open(str(name), 'w+', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(["t", "value", "variance", "origins", "trueN", "sampled_N", "bin_size", "meanBinDate"])
                for i, d in enumerate(days_since3):
                    date = (DAY0 + timedelta(days=d)).strftime("%Y-%m-%d")
                    writer.writerow(
                        [days_since3[i], thetas3[i], variance3[i], origins3[i], num_seqs3[i], num_seqs3[i], bin_sizes3[i],
                         date])


            """
            Analyze with bin merging and subsampling
            """
            print("Analyze with subsampling and bin merging")
            days_since4, thetas4, variance4, num_seqs4, bin_sizes4, subsample_sizes4, origins4 = analyze_run.analyzeBinsWSMB(nb)

            n4 = "WSMB_theta_origins_sim_"+str(SIM_ORDER)+"_numBin_" + str(nb) + ".tsv"
            name = t / n4
            with open(str(name), 'w+', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t',
                                        quotechar='|', quoting=csv.QUOTE_MINIMAL)
                writer.writerow(["t","value","variance","origins","trueN","sampled_N", "bin_size", "meanBinDate"])
                for i,d in enumerate(days_since4):
                    date = (DAY0 + timedelta(days=d)).strftime("%Y-%m-%d")
                    writer.writerow([days_since4[i],thetas4[i],variance4[i],origins4[i],num_seqs4[i],subsample_sizes4[i],bin_sizes4[i],date])


        """
        Analyze with equal sample sizes (100, if sample size is >100, else sample size)
        """
        print("Analyze with subsampling")
        days_since5, thetas5, variance5, num_seqs5, subsample_sizes5, origins5 = analyze_run.analyzeBinsES()

        n2 = "ES_theta_origins_sim_"+str(SIM_ORDER)+".tsv"
        name = t / n2
        with open(str(name), 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["t","value","variance","origins","trueN","sampled_N", "bin_size", "meanBinDate"])
            for i,d in enumerate(days_since2):
                date = (DAY0 + timedelta(days=d)).strftime("%Y-%m-%d")
                writer.writerow([days_since2[i],thetas2[i],variance2[i],origins2[i],num_seqs2[i],subsample_sizes2[i],1,date])


        """
        Plot thetas together to check trajectory
        """
        plt.plot(days_since1, thetas1, 'o', color='green', label="OriginsNoSubsampling")
        plt.plot(days_since2, thetas2, 'o', color='blue', label="OriginsSubsampling")
        plt.plot(days_since3, thetas3, 'o', color='red', label="OriginsMergingBins")
        plt.plot(days_since4, thetas4, 'o', color='yellow', label="OriginsSubsamplingMergingBins")
        plt.plot(days_since5, thetas5, 'o', color='orange', label="OriginsSubsamplingEqualSize")

        plt.xlabel('Time')
        plt.ylabel('Theta')
        plt.legend()
        n = "python_plot_p_repl_" + str(P_REPL) + "_p_mut_" + str(P_MUT) + "_sim_" + str(SIM_ORDER) + "_thetas.png"
        name = t / n
        plt.savefig(str(name))
    
