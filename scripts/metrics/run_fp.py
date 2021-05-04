#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 16:44:06 2020

@author: Maria Trofimova
"""

import os
from pathlib import Path
import numpy as np
import csv
import matplotlib.pyplot as plt
import pandas as pd
import datetime

from bam_to_fingerprints import SAMtoFP
from modular_theta_from_dict import analyzeTrajectory
from parameter_est import parameterEstimation


# Output path
out_dir = Path(snakemake.output[0]).parent
RESULT_PATH = Path(os.getcwd()) / "results"
bins_dir = RESULT_PATH / "bins"
head_dir = bins_dir

# Reference location
reference = str(snakemake.params.ref)
## Name of reference sequence

with open(str(reference), "r") as file:
    header = file.readline()
refname = header.strip(">")
refname = refname.strip("\n")

file_suffix = snakemake.params.group

# Base frequency cutoff
freqCutoff = snakemake.params.cutoff

# Filtering and transformation parameters
min_bin_size = snakemake.params.min_bin_size
min_days_span = snakemake.params.min_days_span
max_days_span = snakemake.params.max_days_span


list_binnings = os.listdir(str(bins_dir))
binnings = []
for file in list_binnings:
    if not file.startswith(".") and not file.endswith(".tsv"):
        binnings.append(file)
binnings.sort()

# Collect all trajectories to merge into one data set
bin_merging_data = []

#List individual binning directories
for folder in binnings:
    binnings_dir = str(bins_dir) + '/' + folder
    headers_dir = str(head_dir) + '/' + folder

    list_files = os.listdir(binnings_dir)
    bam_files = []

    directory = os.listdir(binnings_dir)
    files = []
    for file in directory:
        if file.endswith(".bam"):
            files.append(file)
    files.sort()

    #Initialize results arrays
    seq_list_base_complete = []
    seq_list_pos_complete = []
    seq_list_pairs_complete = []

    seq_dict_int = []
    print('-' * 80)
    print("Current binning folder: ",folder)
    print("      Starting to record mutant positions from CIGAR strings...")
    for filename in files:
        path = binnings_dir+'/'+filename
        sam_to_fp = SAMtoFP(path,reference,refname)
        seq_lbase, seq_posbase, lref, mutants_pairs_list = sam_to_fp.writeFP()
        seq_list_base_complete.append(seq_lbase)
        seq_list_pos_complete.append(seq_posbase)
        seq_list_pairs_complete.append(mutants_pairs_list)

    print("      Done.\n")
    print("      Filtering variants; cutoff of <= %d sequences" % freqCutoff)
    filt = parameterEstimation(seq_list_base_complete,seq_list_pairs_complete,reference,freqCutoff)
    mut_proportion, filtered_seqset = filt.run()
    print("      Done.\n")

    # Time span calculation
    # 1. Get the names of reads from headers and dates
    list_headers = os.listdir(headers_dir)
    headers = []
    for file in list_headers:
        if file.startswith("header"):
            headers.append(file)
    headers.sort()

    # Get the mean date of bin from the sequences
    mean_header_bin = []
    # List of dates that have sequences in the bin
    lists_of_dates = []
    # Indices for bin numbering for internal plotting
    times = []
    # Set of dates that have sequences in the bin
    dates_per_bin = []

    median_dates = []
    max_dates = []
    min_dates = []
    i = 0
    # How many days do sequences from one bin span? - from parameter in config
    # later check if the bin is "long enough"
    num_days_per_bin = []
    for header_file in headers:
        path = headers_dir + '/' + header_file
        table = pd.read_table(path, header=0)
        if not table.empty:
            dates = table['date'].tolist()
            dates_dt = np.array(dates, dtype='datetime64[s]')
            delta_days = (max(dates_dt) - min(dates_dt)) / np.timedelta64(1, 'D') + 1
            num_days_per_bin.append(delta_days)
            mean = (np.array(dates, dtype='datetime64[s]').view('i8').mean().astype('datetime64[s]'))
            # collect also median and max date of bin
            median = np.median(np.array(dates, dtype='datetime64[s]').view('i8')).astype('datetime64[s]')
            minDate = (np.array(dates, dtype='datetime64[s]').view('i8').min().astype('datetime64[s]'))
            maxDate = (np.array(dates, dtype='datetime64[s]').view('i8').max().astype('datetime64[s]'))
            median_dates.append(str(median)[:10])
            max_dates.append(str(maxDate)[:10])
            min_dates.append(str(minDate)[:10])
            # print(mean)
            lists_of_dates.append(dates)
            dates_per_bin.append(set(dates))
            mean_header_bin.append(str(mean)[:10])
            times.append(i)
            i += 1

    print("      Starting to compute optimal metric parameters...")
    #Theta from origins - MLE
    analyze = analyzeTrajectory(filtered_seqset, mut_proportion, num_days_per_bin, '')
    thetas, variance, variance_size, num_seqs, num_mut, origins = analyze.analyzeBinsMLE()
    weeks = np.arange(0,len(thetas))
    print("      Done.\n")
    print("      Starting to write tables and making complementary plots...")

    # Plots for different variables and outcomes of individual binnings
    # a) Bin sizes bar plot
    plot_title = "Bin sizes for %s" % folder
    fig, ax1 = plt.subplots()
    ax1.set_title("%s" % plot_title)
    ax1.set_xlabel('Mean bin date')
    ax1.set_ylabel("Bin size", color='crimson')
    ax1.bar(mean_header_bin, num_seqs, color='crimson')
    ax1.tick_params(axis='y', labelcolor='crimson')
    ax1.set_xticks(mean_header_bin[::10])
    ax1.set_xticklabels(mean_header_bin[::10])
    fig.tight_layout()
    name = str(out_dir)+"/plot_"+folder+"_bin_sizes.png"
    fig.savefig(str(name),dpi=300)
    plt.clf()

    # b) Number of origins vs number of mutants
    plot_title = "Origins and mutants for %s" % folder
    fig, ax1 = plt.subplots()
    ax1.set_title("%s" % plot_title)
    ax1.plot(mean_header_bin, origins, '-', color='green', label="Origins")
    ax1.plot(mean_header_bin, num_mut, 'o', color='blue', label="Number of mutants")
    ax1.set_xlabel('Mean bin date')
    ax1.set_ylabel('Count')
    ax1.set_xticks(mean_header_bin[::10])
    ax1.set_xticklabels(mean_header_bin[::10])
    ax1.legend()
    name = str(out_dir)+"/plot_"+folder+"_originsvsmut.png"
    fig.savefig(str(name),dpi=300)
    plt.clf()

    # c) Theta estimates normalized by time delta
    plot_title = "Population size estimate %s" % folder
    fig, ax1 = plt.subplots()
    ax1.set_title("%s" % plot_title)
    ax1.plot(mean_header_bin, np.array(thetas), '-', color='royalblue')
    ax1.set_xlabel('Mean bin date')
    ax1.set_ylabel(r'$\phi_{est}$')
    ax1.set_xticks(mean_header_bin[::10])
    ax1.set_xticklabels(mean_header_bin[::10])
    name = str(out_dir)+"/plot_"+folder+"_thetas.png"
    fig.savefig(str(name),dpi=300)
    fig.clf()

    # Write bin file
    name_table = "table_"+file_suffix+"_"+folder+"_phi_estimate_var_from_size.tsv"
    table_path = str(out_dir) + '/' + name_table
    with open(table_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["date","value","variance_size","variance_mle","num_seqs","times",
                         "num_days_per_bin", "median_date", "min_date", "max_date"])
        for i in range(len(weeks)):
            writer.writerow([mean_header_bin[i], thetas[i], variance_size[i], variance[i], num_seqs[i], times[i],
                             num_days_per_bin[i], median_dates[i], min_dates[i], max_dates[i]])
    print("      Done.\n")
    # Write to merged bins dataset
    for i, date in enumerate(mean_header_bin):
        if not (folder.startswith("fuzzy")):
            if not i+1 == len(mean_header_bin):
                bin_merging_data.append((mean_header_bin[i], thetas[i], variance_size[i], variance[i], num_seqs[i],
                                         times[i], num_days_per_bin[i], median_dates[i], max_dates[i], folder))


print("\n      Making the final results table...")
# Make merged dataset sorted by date
bin_merging_data_ = sorted(bin_merging_data)
times = np.arange(0, len(bin_merging_data_))
name_table = "table_merged_phi_estimates_var_from_size.tsv"

table_path = str(out_dir) + '/' + name_table

with open(table_path, 'w+', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["t","value","variance","meanBinDate","sampleSize",
                     "daysPerBin", "medianBinDate", "maxBinDate", "binning"])
    for i in range(len(times)):
        # If bin size==1/variance is bigger or equal to min_bin_size
        # if maxsize, minsize satisfied
        if bin_merging_data_[i][2]<=1/int(min_bin_size) and bin_merging_data_[i][6]>=min_days_span and bin_merging_data_[i][6]<=max_days_span:
            writer.writerow([times[i],bin_merging_data_[i][1],bin_merging_data_[i][2], bin_merging_data_[i][0],
                            bin_merging_data_[i][4],bin_merging_data_[i][6],bin_merging_data_[i][7],
                            bin_merging_data_[i][8],bin_merging_data_[i][9]])

print("Done.\n")
