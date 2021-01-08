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

from bam_to_fingerprints import SAMtoFP

from modular_theta_from_dict import analyzeTrajectory

from parameter_est import parameterEstimation

from msn import MSN

from lookup_tree import treeLookup

import matplotlib.pyplot as plt

import pandas as pd

import datetime

from Bio import SeqIO

import copy



#plt.style.use('ggplot')
#font = {'family' : 'sans-serif',
#        'weight' : 'normal',
#        'size'   : 10}
#
#matplotlib.rc('font', **font)


FILEPATH = Path(__file__).parent
FILEPATH_interm = FILEPATH.parent
bins_dir = FILEPATH_interm.parent / "results" / "bins"
head_dir = FILEPATH_interm.parent / "results" / "bins"
# Output path
od = snakemake.output[0]
od_split = od.split("/table")
out_dir = FILEPATH_interm.parent / od_split[0]

# Reference location
ref = snakemake.params.ref + '.fasta'
reference = FILEPATH_interm.parent / 'consensus' / ref

# Length of reference
ref_fasta = SeqIO.parse(open(reference),'fasta')
ref_seq = ''
for fasta in ref_fasta:
    ref_seq = str(fasta.seq)
lref = len(ref_seq)

# Reported cases data
table_name = snakemake.params.rep_cases[0]
table_delim = snakemake.params.rep_cases[1]
table_date_col = snakemake.params.rep_cases[2]
table_active_col = snakemake.params.rep_cases[3]
table_date_format = snakemake.params.rep_cases[4]
table_path = FILEPATH_interm.parent / "reported_cases" / str(table_name)

# Filtering and transformation parameters
min_bin_size = snakemake.params.min_bin_size
min_days_span = snakemake.params.min_days_span
max_days_span = snakemake.params.max_days_span
smoothing = snakemake.params.smoothing


WORKING_PATH = FILEPATH_interm.parent

list_binnings = os.listdir(str(bins_dir))
binnings = []
for file in list_binnings:
    if not file.startswith(".") and not file.endswith(".tsv"):
        binnings.append(file)
binnings.sort()

#Load real data table - from config
table = pd.read_table(str(table_path),delimiter=str(table_delim),header=0)
dates = table[table_date_col].tolist()
new_positives = table[table_active_col].tolist()

# Reformat dates
dates_nt = []
for i, date in enumerate(dates):
    date_r = datetime.datetime.strptime(date, table_date_format).strftime('%Y-%m-%d')
    dates_nt.append(date_r)

bin_merging_data = []

# Save clusters from first run
iteration = 0
components_names_global = []

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
    print("Doing folder ",folder)
    print("      Starting to record mutant positions from CIGAR strings...")
    for filename in files:
        path = binnings_dir+'/'+filename
        sam_to_fp = SAMtoFP(path,reference)
        seq_lbase, seq_posbase, lref, mutants_pairs_list = sam_to_fp.writeFP()
        seq_list_base_complete.append(seq_lbase)
        seq_list_pos_complete.append(seq_posbase)
        seq_list_pairs_complete.append(mutants_pairs_list)
        
    print("      Starting to filter out low frequency variant positions and inferring total mutant population... ")
    
    filt = parameterEstimation(seq_list_base_complete,seq_list_pairs_complete,reference)
    mut_proportion, filtered_seqset = filt.run()
    print("           Calculated overall mutant proportion: ",mut_proportion)
    
    print("      Writing allele frequencies")
    basesDict = dict()
    for i, seqSet in enumerate(seq_list_pairs_complete):
        for j,seq in enumerate(seqSet):
            for bases in seq:
                if bases[0] in basesDict.keys():
                    basesDict[bases[0]][i] += 1
                else:
                    traj = np.zeros(len(seq_list_pairs_complete))
                    traj[i] += 1
                    basesDict[bases[0]] = traj
                    
    n1 = "allele_frequencies_"+folder+".tsv"
    table_path = str(out_dir) + '/' + n1
    header_row = list(np.arange(len(basesDict)))
    with open(table_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(list(basesDict.keys()))
        for i, d in enumerate(seq_list_pairs_complete):
            row = []
            for (key,value) in basesDict.items():
                row.append(value[i])
            writer.writerow(row)

    print("      Done.\n")
    
    print("      Clustering and MSN")
    components_names = []
    if iteration==0:
    
        msn = MSN(filtered_seqset, ref_seq)
        msn.buildGraph() 
        result, ids_dict = msn.KruskalMST()
        # Connected components with graph IDs and dictionary with corresponding sequence IDs
        # and clustered network edges 
        components, components_dict, clusters_edges = msn.KruskalClustering(result)
    
        inv_ids_dict = {v: k for k, v in ids_dict.items()}
        #Make dict from graph ID to FP
        fp_seqname_dict = dict()
        for ii in range(len(seq_list_base_complete)):
            for jj in range(len(seq_list_base_complete[ii])):
                fp_seqname_dict[seq_list_base_complete[ii][jj][0]] = seq_list_base_complete[ii][jj][1]
        print("inv_seqnames_dict")
        # Clusters of sequence identifiers
        for i in range(len(components)):
            c = []
            for j in range(len(components[i])):
                c.append(inv_ids_dict[components[i][j]])
            components_names.append(c)
            
        components_names_global = components_names
        # get sequence names by ID from components_dict
        components_seqnames = []
        for i in range(len(components)):
            cc = []
            for j in range(len(components[i])):
                name = components_dict[components[i][j]]
                cc.append(name)
            components_seqnames.append(cc)
        components_fps = []
        for i in range(len(components_seqnames)):
            cc = []
            for j in range(len(components_seqnames[i])):
                name = inv_ids_dict[components[i][j]]
                for ii in range(len(seq_list_base_complete)):
                    for jj in range(len(seq_list_base_complete[ii])):
                        if seq_list_base_complete[ii][jj][0]==name:
                            cc.append(seq_list_base_complete[ii][jj][1])
            components_fps.append(cc)
            
        name_table = "table_msn_"+folder+"_edges_thetas_var_from_size.tsv"
        table_path = str(out_dir) + '/' + name_table
        nodes_names = dict()
        with open(table_path, 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(["from","to","weight"])
            for i in range(len(result)):
                fr = inv_ids_dict[result[i][0]]
                to = inv_ids_dict[result[i][1]]
                nodes_names[fr] = components_dict[result[i][0]]
                nodes_names[to] = components_dict[result[i][1]]
                w = result[i][2]
                writer.writerow([fr,to,w])
                
        name_table = "table_msn_"+folder+"_nodes_thetas_var_from_size.tsv"
        table_path = str(out_dir) + '/' + name_table
        with open(table_path, 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(["node","label"])
            for (key,value) in nodes_names.items():
                writer.writerow([key,value])
                
        
        cl_edges = []
        name_table = "table_msn_"+folder+"_clusters_edges_thetas_var_from_size.tsv"
        table_path = str(out_dir) + '/' + name_table
        nodes_cl = dict()
        with open(table_path, 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(["from","to","weight","l1","l2"])
            for i in range(len(result)):
                fr = inv_ids_dict[result[i][0]]
                to = inv_ids_dict[result[i][1]]
                nodes_cl[fr] = components_dict[result[i][0]]
                nodes_cl[to] = components_dict[result[i][1]]
                w = result[i][2]
                if not (nodes_cl[fr]==nodes_cl[to]):
                    writer.writerow([nodes_cl[fr],nodes_cl[to],w])
                    cl_edges.append([nodes_cl[fr],nodes_cl[to],w])
                 
        cl_nodes = []
        name_table = "table_msn_"+folder+"_clusters_nodes_thetas_var_from_size.tsv"
        table_path = str(out_dir) + '/' + name_table
        with open(table_path, 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(["node","size"])
            for i in range(len(components)):
                writer.writerow([i,len(components[i])])
                cl_nodes.append([i,len(components[i])])
                
        # Get all FPs in cluster
        cl_fps = []
        name_table = "table_msn_"+folder+"_clusters_sequences_thetas_var_from_size.tsv"
        table_path = str(out_dir) + '/' + name_table
        with open(table_path, 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(["node","sequence","locations"])
            for i in range(len(components_fps)):
                mutations = []
                locs = []
                for j in range(len(components_fps[i])):
                    split = components_fps[i][j].split("-")
                    name = inv_ids_dict[components[i][j]]
                    locs.append(name[16:18])
                    for s in split:
                        mutations.append(s)
                umutations = np.unique(np.array(sorted(mutations)))
                mutantsstrbase = "-".join(umutations)
                locations = "-".join(locs)
                cl_fps.append(mutantsstrbase)
                writer.writerow([i,mutantsstrbase,locations])

        # Find diconnected variants
        #fps_unique, fps_components_clusters = msn.getConnectedMutClusters(cl_edges,cl_nodes,cl_fps)
        print("Found connected mutant clusters:")
    else:
        components_names = components_names_global
        
    #print("           Density based clustering")
    #clustering = parameterEstimation(seq_list_base_complete,seq_list_pairs_complete,reference)
    #components_names = clustering.runHapl()
            
    print("           Calculated overall mutant proportion: ",mut_proportion)
    
    print("      Done.\n")
    print("      Starting to compute optimal metric parameters...")
    #Theta from origins - MLE
    analyze = analyzeTrajectory(filtered_seqset, mut_proportion, '')
    #thetas, variance, variance_size, num_seqs, num_mut, origins = analyze.analyzeBinsMLE()
    
    thetas, variance, variance_size, num_seqs, num_mut, origins, thetas_e = analyze.analyzeBinsMLEcl(components_names)
    # Theta from phylogeny - homoplasic sites
    #analyze = analyzeTrajectory(filtered_seqset, mut_proportion, '', tree)
    #thetas, variance, variance_size, num_seqs, num_mut, origins = analyze.analyzeBinsMLE()
    
    weeks = np.arange(0,len(thetas))
    
    print("      Done.\n")
    
    print("      Calculating mutant positions trajectories... ")
    filt = parameterEstimation(seq_list_base_complete,seq_list_pairs_complete,reference)
    trajectories, allele_trajectories = filt.runMut()
    print("      Done.\n")  
    
            
    
    print("      Starting to write tables and making complementary plots...")
    #1. Get the names of reads from headers and dates
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
    i = 0
    # How many days do sequences from one bin span? - from parameter in config 
    # later check if the bin is "long enough"
    num_days_per_bin = []
    for header_file in headers:
        path =  headers_dir+ '/'+header_file
        table = pd.read_table(path, header=0)
        if not table.empty:
            dates = table['date'].tolist()
            dates_dt = np.array(dates, dtype='datetime64[s]')
            delta_days = (max(dates_dt)-min(dates_dt))/np.timedelta64(1,'D')
            num_days_per_bin.append(delta_days)
            mean = (np.array(dates, dtype='datetime64[s]').view('i8').mean().astype('datetime64[s]'))
            #print(mean)
            lists_of_dates.append(dates)
            dates_per_bin.append(set(dates))
            mean_header_bin.append(str(mean)[:10])
            times.append(i)
            i+=1

    cases_on_mean_date = []
    for i, date in enumerate(mean_header_bin):
        if date in dates_nt:
            index = dates_nt.index(date)
            cases_on_mean_date.append(new_positives[index])
        else:
            cases_on_mean_date.append(0)
    
    
    # 2) Get the total number of cases on dates that have sequences
    cases_on_all_dates_in_seqbin = []
    for i, dates in enumerate(dates_per_bin):
        cases_on_all_dates_in_seqbin_tosum = []
        for j, date in enumerate(dates):
            if date in dates_nt:
                index = dates_nt.index(date)
                cases_on_all_dates_in_seqbin_tosum.append(new_positives[index])
            else:
                cases_on_all_dates_in_seqbin_tosum.append(0)
        cases_on_all_dates_in_seqbin.append(sum(cases_on_all_dates_in_seqbin_tosum))
        
    # Write table for SNP trajectories
    traj_snps = []
    for (key,value) in trajectories.items():
        traj = []
        traj.append(key)
        for i in range(len(value)):
            traj.append(value[i])
        traj_snps.append(traj)
    #print(traj_snps)
    header_alleles = copy.deepcopy(mean_header_bin)
    header_alleles.insert(0, "snp_name")
    name_table = "table_snps_"+folder+"_thetas_var_from_size.tsv"
    table_path = str(out_dir) + '/' + name_table
    with open(table_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header_alleles)
        for i in range(len(traj_snps)):
            wr = traj_snps[i]
            #wr.insert(0, alleles_fps[i])
            writer.writerow(wr)
        

    # PLots for different variables and outcomes of individual binnings
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
    plot_title = "Mutants for %s" % folder
    fig, ax1 = plt.subplots()
    ax1.set_title("%s" % plot_title)
    #ax1.plot(mean_header_bin, origins, '-', color='green', label="Origins")
    ax1.plot(mean_header_bin, num_mut, 'o', color='blue', label="Number of mutants")
    ax1.set_xlabel('Mean bin date')
    ax1.set_ylabel('Count')
    ax1.set_xticks(mean_header_bin[::10])
    ax1.set_xticklabels(mean_header_bin[::10])
    ax1.legend()
    name = str(out_dir)+"/plot_"+folder+"_mut.png"
    fig.savefig(str(name),dpi=300)
    plt.clf()
    
    plot_title = "Origins for %s" % folder
    fig, ax1 = plt.subplots()
    ax1.set_title("%s" % plot_title)
    ax1.plot(mean_header_bin, origins, '-', color='green', label="Origins")
    #ax1.plot(mean_header_bin, num_mut, 'o', color='blue', label="Number of mutants")
    ax1.set_xlabel('Mean bin date')
    ax1.set_ylabel('Count')
    ax1.set_xticks(mean_header_bin[::10])
    ax1.set_xticklabels(mean_header_bin[::10])
    ax1.legend()
    name = str(out_dir)+"/plot_"+folder+"_origins.png"
    fig.savefig(str(name),dpi=300)
    plt.clf()
    
    # c) Theta estimates normalized by time delta
    plot_title = "Population size estimate %s" % folder
    fig, ax1 = plt.subplots()
    ax1.set_title("%s" % plot_title)
    ax1.plot(mean_header_bin, np.array(thetas), '-', color='royalblue')
    ax1.set_xlabel('Mean bin date')
    ax1.set_ylabel(r'$\theta_{est}$')
    ax1.set_xticks(mean_header_bin[::10])
    ax1.set_xticklabels(mean_header_bin[::10])
    name = str(out_dir)+"/plot_"+folder+"_thetas.png"
    fig.savefig(str(name),dpi=300)
    fig.clf()  
    
    # Write bin file
    name_table = "table_"+folder+"_thetas_var_from_size.tsv"
    table_path = str(out_dir) + '/' + name_table
    with open(table_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["date","value","variance_size","variance_mle","origins_from_clusters","num_seqs","times","cases_on_mean_date","cases_on_all_dates_in_seqbin","num_days_per_bin","thetaEwens"])
        for i in range(len(weeks)):
            writer.writerow([mean_header_bin[i], thetas[i], variance_size[i], variance[i], origins[i], num_seqs[i], times[i], cases_on_mean_date[i], cases_on_all_dates_in_seqbin[i], num_days_per_bin[i],thetas_e[i]])
    print("      Done.\n")
    # Write to merged bins dataset
    for i, date in enumerate(mean_header_bin):
        if not (folder.startswith("fuzzy")):
            if not i+1 == len(mean_header_bin):
                bin_merging_data.append((mean_header_bin[i], thetas[i], variance_size[i], variance[i], num_seqs[i], times[i], cases_on_mean_date[i], cases_on_all_dates_in_seqbin[i], num_days_per_bin[i]))
    
    iteration += 1

# Make merged dataset sorted by date
bin_merging_data_ = sorted(bin_merging_data)
times = np.arange(0, len(bin_merging_data_))
# Write arrays just in case for individual plots and files
#TODO: write file directly from tuple
datesm = []
num_days_per_binm = []
thetasm = []
variance_sizem = []
variancem = []
num_seqsm = []
timesm = []
cases_on_mean_datem = []
cases_on_all_dates_in_seqbinm = []
num_days_per_binm = []


for i, mbin in enumerate(bin_merging_data_):
    datesm.append(mbin[0])
    thetasm.append(mbin[1])
    variance_sizem.append(mbin[2])
    variancem.append(mbin[3])
    num_seqsm.append(mbin[4])
    timesm.append(mbin[5])
    cases_on_mean_datem.append(mbin[6])
    cases_on_all_dates_in_seqbinm.append(mbin[7])
    num_days_per_binm.append(mbin[8])


print("\nMaking the final results table...")
name_table = "table_merged_thetas_var_from_size.tsv"

table_path = str(out_dir) + '/' + name_table

with open(table_path, 'w+', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["t","value","value_nonorm","variance","trueN","meanBinDate","sampleSize"])
    for i in range(len(times)):
        # If bin size==1/variance is bigger or equal to min_bin_size
        #if maxsize, minsize satisfied 
        if variance_sizem[i]<=1/int(min_bin_size) and num_days_per_binm[i]>=min_days_span and num_days_per_binm[i]<=max_days_span:
            writer.writerow([times[i],thetasm[i],thetasm[i],variance_sizem[i],cases_on_mean_datem[i],datesm[i],num_seqsm[i]])
        
print("Done.\n")