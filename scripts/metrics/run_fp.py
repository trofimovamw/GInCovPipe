#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 16:44:06 2020

@author: Maria Trofimova
"""

import os

import re

from pathlib import Path

import numpy as np

import csv

from bam_to_fingerprints import SAMtoFP

from modular_theta_from_dict import analyzeTrajectory

import matplotlib.pyplot as plt

import pandas as pd

from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

import datetime


plt.style.use('ggplot')
#font = {'family' : 'sans-serif',
#        'weight' : 'normal',
#        'size'   : 10}
#
#matplotlib.rc('font', **font)


FILEPATH = Path(__file__).parent
FILEPATH_interm = FILEPATH.parent
bins_dir = FILEPATH_interm.parent / "results" / "fixed_cigars_bins"
head_dir = FILEPATH_interm.parent / "results" / "bins"
# Output path
od = snakemake.output[0]
od_split = od.split("/table")
out_dir = FILEPATH_interm.parent / od_split[0]
# Reference location
ref = snakemake.params.ref[0]
reference = FILEPATH_interm.parent / ref
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
smoothing = snakemake.params.smoothing
# Metadata
#metadata = snakemake.params.meta
#meta_path = FILEPATH_interm.parent

WORKING_PATH = FILEPATH_interm.parent


list_binnings = os.listdir(str(bins_dir))
binnings = []
for file in list_binnings:
    if not file.startswith(".") and not file.endswith(".tsv"):
        binnings.append(file)
binnings.sort()


'''
Load real data table - from config
'''
table = pd.read_table(str(table_path),delimiter=str(table_delim),header=0)
dates = table[table_date_col].tolist()
new_positives = table[table_active_col].tolist()
total_positives = table[table_active_col].tolist()
# Reformat dates
dates_nt = []
for i, date in enumerate(dates):
    date_r = datetime.datetime.strptime(date, table_date_format).strftime('%Y-%m-%d')
    dates_nt.append(date_r)

bin_merging_data = []

# Parameters for smoothing

window = 15
poly = 2

"""
List individual binning directories
"""
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

    """
    Initialize results arrays
    """
    
    seq_list_base_complete = []
    

    seq_dict_int = []
    print("Doing folder ",folder)
    for filename in files:
        path = binnings_dir+'/'+filename
        sam_to_fp = SAMtoFP(path)
        seq_lbase = sam_to_fp.writeFP()
        seq_list_base_complete.append(seq_lbase)
        

    '''
    Theta from origins - MLE
    '''

    analyze = analyzeTrajectory(seq_list_base_complete, '')
    thetas, variance, variance_size, num_seqs, num_mut, origins = analyze.analyzeBinsMLE()
    weeks = np.arange(0,len(thetas))

    '''
    Theta from segregating sites
    '''
    #analyze2 = analyzeTrajectory(seq_dict_int, '')
    #thetas = analyze2.analyzeBinsSegrS()

    """
    Plot thetas together to check trajectory

    1. Get the names of reads from headers and dates
    """
    list_headers = os.listdir(headers_dir)
    headers = []
    for file in list_headers:
        if file.startswith("header"):
            headers.append(file)
    headers.sort()

    mean_header_bin = []
    lists_of_dates = []
    times = []
    dates_per_bin = []
    i = 0
    #How many days do sequences from one bin span? - for parameter in config
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
    

    """
    2. Get three counts from real data:
        a) total_positive on mean date
    """

    rep_cases_a = []
    control_dates_a = []
    for i, date in enumerate(mean_header_bin):
        if date in dates_nt:
            index = dates_nt.index(date)
            rep_cases_a.append(total_positives[index])
            control_dates_a.append(dates_nt[index])
        else:
            rep_cases_a.append(0)

    """
        b) new_positive on mean date
    """
    rep_cases_b = []
    control_dates_b = []
    for i, date in enumerate(mean_header_bin):
        if date in dates_nt:
            index = dates_nt.index(date)
            rep_cases_b.append(new_positives[index])
            control_dates_b.append(dates_nt[index])
        else:
            rep_cases_b.append(0)

    """
        c) new_positive sum on all dates on bin
    """
    rep_cases_c = []
    control_dates_c = []
    for i, ldates in enumerate(lists_of_dates):
        cases_list = []
        cdates = []
        for j,date in enumerate(ldates):
            if date in dates_nt:
                index = dates_nt.index(date)
                cases_list.append(new_positives[index])
                cdates.append(dates_nt[index])
            else:
                cases_list.append(0)
        sum_cases = np.sum(np.array(cases_list))
        rep_cases_c.append(sum_cases)
    """
        d) total_positive per bin from change in total_positive number
    """
    rep_cases_d = []
    delta_days = []
    for i, ldates in enumerate(lists_of_dates):
        cases_list = []
        cdates = []
        sorted_dates = ldates.sort()
        # Get range of dates for bin
        start_date = datetime.datetime.strptime(ldates[0], '%Y-%m-%d')
        end_date = datetime.datetime.strptime(ldates[-1], '%Y-%m-%d')
        delta = end_date-start_date
        delta_days.append(int(delta.days)+1)
        step = datetime.timedelta(days=1)
        dates_range = []
        dates_range.append(start_date.strftime('%Y-%m-%d'))
        while start_date < end_date:
            start_date += step
            dates_range.append(start_date.strftime('%Y-%m-%d'))

        # Total positive per bin == mean of total_positives on dates in range
        total_pos = []
        for m,date in enumerate(dates_range):
            if date in dates_nt:
                index = dates_nt.index(date)
                new_tot_pos = total_positives[index]
                total_pos.append(new_tot_pos)
        rep_cases_d.append(np.sum(np.array(total_pos))/len(total_pos))

    """
        e) total_positive per bin from smoothed interval of real cases on
           dates between first and last bin date
    """
    rep_cases_smooth_range = []
    # Smooth through the dates that are in the list and take value at mean
    for i, ldates in enumerate(lists_of_dates):
        cases_list = []
        cdates = []
        sorted_dates = ldates.sort()
        # Get range of dates for bin
        start_date = datetime.datetime.strptime(ldates[0], '%Y-%m-%d')
        end_date = datetime.datetime.strptime(ldates[-1], '%Y-%m-%d')
        delta = end_date-start_date
        delta_days.append(int(delta.days)+1)
        step = datetime.timedelta(days=1)
        dates_range = []
        dates_range.append(start_date.strftime('%Y-%m-%d'))
        while start_date < end_date:
            start_date += step
            dates_range.append(start_date.strftime('%Y-%m-%d'))

        # Total positive per bin == positive on date in range
        total_pos = []
        ex_dates = []
        for m,date in enumerate(dates_range):
            if date in dates_nt:
                index = dates_nt.index(date)
                new_tot_pos = total_positives[index]
                total_pos.append(new_tot_pos)
                ex_dates.append(date)
        if len(ex_dates) > 1:
            x = np.arange(0,len(ex_dates))
            xx = np.linspace(x.min(),x.max(), 1000)
            itp = interp1d(x,total_pos, kind='linear')
            window_size, poly_order = 5, 3
            smooth_total_pos = savgol_filter(itp(xx),window_size, poly_order)

            smooth_mean_tot_pos = 0
            if len(ex_dates)>1:
                for j in range(1,len(ex_dates)):
                    if mean_header_bin[i] >= ex_dates[j-1] and mean_header_bin[i] < ex_dates[j]:
                        smooth_mean_tot_pos = (smooth_total_pos[j-1]+smooth_total_pos[j])/2
            else:
                smooth_mean_tot_pos = smooth_total_pos[0]
            rep_cases_smooth_range.append(smooth_mean_tot_pos)
        elif len(ex_dates) == 1:
            rep_cases_smooth_range.append(total_pos[0])
        else:
            rep_cases_smooth_range.append(0)

    """
        f) total_positive per bin from smoothed interval the complete bin,
            as defined by binning routines - from files with "range" prefix
    """
    # TODO:
    # check-check-check
    rep_cases_smooth_wbin = []

    # Smooth the complete data - 7 smoothing window
    smoothed_total_positives = []
    if smoothing==True:
        smoothed_total_positives = savgol_filter(total_positives, window, poly)
    else:
        smoothed_total_positives = total_positives


    if folder.startswith("eq_size"):
        rep_cases_smooth_wbin = rep_cases_smooth_range

    else:
        list_ranges = os.listdir(headers_dir)
        ranges = []
        for file in list_ranges:
            if file.startswith("range"):
                ranges.append(file)
        ranges.sort()

        for d,range_file in enumerate(ranges):
            path =  headers_dir+'/'+range_file
            table = pd.read_table(path, header=0)
            if not table.empty:
                days_range = []

                # List all dates in range
                start_date = table['start_day'].tolist()[0]
                end_date = table['end_day'].tolist()[0]
                delta = datetime.datetime.strptime(end_date, '%Y-%m-%d')-datetime.datetime.strptime(start_date, '%Y-%m-%d')

                for i in range(delta.days + 1):
                    days_range.append((datetime.datetime.strptime(start_date, '%Y-%m-%d') + datetime.timedelta(days=i)).strftime('%Y-%m-%d'))

                mean = str(np.array(days_range, dtype='datetime64[s]').view('i8').mean().astype('datetime64[s]'))[:10]

                # Get the available data points on qeach date in the bin range
                smoothing_data = []
                ex_dates = []

                for m, date in enumerate(days_range):
                    if date in dates_nt:
                        index = dates_nt.index(date)
                        new_tot_pos = total_positives[index]
                        smoothing_data.append(new_tot_pos)
                        ex_dates.append(date)

                index = np.where(np.array(dates_nt)==mean)[0]
                if len(index)>0:
                    mean_rep_cases = smoothed_total_positives[index[0]]
                    rep_cases_smooth_wbin.append(mean_rep_cases)
                else:
                    rep_cases_smooth_wbin.append(0)

        print(rep_cases_smooth_wbin)


    """
    3. Plot two axes on the same graphic
    """
    print("Plotting")
    if not os.path.exists(str(out_dir)):
        os.makedirs(str(out_dir))
    
    

    plot_title = folder
    if folder.startswith('cal_week'):
        plot_title = "Binning by calendar week"
    if folder.startswith('eq_days'):
        m = re.search(r'\d+', folder)
        plot_title = "Binning by "+m[0]+" days"
    if folder.startswith('eq_size'):
        m = re.search(r'\d+', folder)
        plot_title = "Binning by "+m[0]+" sequences"
    if folder.startswith('fuzzy_days'):
        m = re.search(r'\d+', folder)
        plot_title = "Fuzzy binning by "+m[0]+" days"
        
    # Date formatting for plots
    dates = np.array([datetime.datetime.strptime(d,"%Y-%m-%d").date() for d in mean_header_bin])

    #locator = mdates.DayLocator()
    
    fig, ax1 = plt.subplots()
    ax1.set_title("%s" % plot_title)
    ax1.set_xlabel('Mean bin date')
    ax1.set_ylabel(r'$\theta_{est}$', color='crimson')
    ax1.errorbar(mean_header_bin, thetas, yerr=np.array(variance_size)*10, fmt='.', elinewidth=0.7, color='crimson')
    ax1.tick_params(axis='y', labelcolor='crimson')
    ax1.set_xticks(mean_header_bin[::3])
    ax1.set_xticklabels(mean_header_bin[::3],rotation=45)
    #ax1.xaxis.set_major_locator(locator)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Sample size', color='royalblue')
    ax2.plot(weeks, num_seqs, color='royalblue')
    ax2.tick_params(axis='y', labelcolor='royalblue')

    fig.tight_layout()
    #plt.show()

    name = str(out_dir)+"/plot_"+folder+"_thetas.png"
    fig.savefig(str(name),dpi=300)
    plt.clf()

    fig, ax1 = plt.subplots()
    ax1.set_title("%s" % plot_title)
    ax1.set_xlabel('Mean bin date')
    ax1.set_ylabel(r'$\theta_{est}$', color='crimson')
    ax1.errorbar(mean_header_bin, thetas, yerr=np.array(variance_size)*10, fmt='.', elinewidth=0.7, color='crimson')
    ax1.tick_params(axis='y', labelcolor='crimson')
    ax1.set_xticks(mean_header_bin[::3])
    ax1.set_xticklabels(mean_header_bin[::3],rotation=45)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Active cases (smoothed bin)', color='forestgreen')
    ax2.plot(weeks, rep_cases_smooth_range, color='forestgreen')
    ax2.tick_params(axis='y', labelcolor='forestgreen')

    fig.tight_layout()
    #plt.show()

    name = str(out_dir)+"/plot_"+folder+"_totpos_cases_smooth_bin_range.png"
    fig.savefig(str(name),dpi=300)
    plt.clf()

    fig, ax1 = plt.subplots()
    ax1.set_title("%s" % plot_title)
    ax1.set_xlabel('Mean bin date')
    ax1.set_ylabel(r'$\theta_{est}$', color='crimson')
    ax1.errorbar(mean_header_bin, thetas, yerr=np.array(variance_size)*10, fmt='.', elinewidth=0.7, color='crimson')
    ax1.tick_params(axis='y', labelcolor='crimson')
    ax1.set_xticks(mean_header_bin[::3])
    ax1.set_xticklabels(mean_header_bin[::3],rotation=45)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Active cases (%s-day window, poly-order %s)' % (window,poly), color='forestgreen')
    ax2.plot(weeks, rep_cases_smooth_wbin, color='forestgreen')
    ax2.tick_params(axis='y', labelcolor='forestgreen')

    fig.tight_layout()
    #plt.show()

    file_insert = "_totpos_cases_on_smooth_poly_window_%s_poly_%s_.png" % (window,poly)
    name = str(out_dir)+"/plot_"+folder+file_insert
    fig.savefig(str(name),dpi=300)
    plt.clf()

    plt.title("%s" % plot_title)
    plt.plot(weeks, origins, '-', color='green', label="Origins")
    plt.plot(weeks, num_mut, 'o', color='royalblue', label="Number of mutants")
    plt.xlabel('Mean bin date')
    plt.ylabel('Count')
    plt.xticks(range(0,len(mean_header_bin)), mean_header_bin)
    plt.xticks(rotation=45)
    plt.legend()

    name = str(out_dir)+"/plot_"+folder+"_originsvsmut.png"
    plt.savefig(str(name),dpi=300)
    plt.cla()




    """
    4. Write table - with mean dates from headers; to compare with reported cases
    """

    name_table = "table_"+folder+"_thetas_var_from_size.tsv"
    table_path = str(out_dir) + '/' + name_table
    with open(table_path, 'w+', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["t","value","variance","binSize","daysPerBin","numberOfOrigins","meanBinDate","newCases_total_positive_on_mean","newCases_total_positive_mean_on_all","newCases_new_positive_on_mean","newCases_new_positive_on_all","newCases_smooth_bin","listOfDates"])
        for i in range(len(weeks)):
            writer.writerow([weeks[i],thetas[i],variance_size[i],num_seqs[i],delta_days[i],origins[i],mean_header_bin[i],rep_cases_a[i],rep_cases_d[i],rep_cases_b[i],rep_cases_c[i],rep_cases_smooth_range[i],lists_of_dates[i]])

    """
    Bin merging
    """
    for i, date in enumerate(mean_header_bin):
        #if not folder.startswith("eq_size"):
        bin_merging_data.append((date, thetas[i], variance_size[i], variance[i], num_seqs[i], times[i], rep_cases_smooth_range[i], rep_cases_a[i], rep_cases_smooth_wbin[i], num_days_per_bin[i], folder))

"""
Sort and write the merged bins array
"""
bin_merging_data_ = sorted(bin_merging_data)
weeks = np.arange(0, len(bin_merging_data_))
variance_m = []
variance_sm = []
num_seqs_m = []
thetas_m = []
dates_m = []
times = []
rep_smooth = []
totpos_onmean = []
rep_cases_wbin = []
num_days_wbin = []
folders = []

for i, mbin in enumerate(bin_merging_data_):
    variance_m.append(mbin[3])
    variance_sm.append(mbin[2])
    num_seqs_m.append(mbin[4])
    thetas_m.append(mbin[1])
    dates_m.append(mbin[0])
    times.append(mbin[-6])
    rep_smooth.append(mbin[-5])
    totpos_onmean.append(mbin[-4])
    rep_cases_wbin.append(mbin[-3])
    num_days_wbin.append(mbin[-2])
    folders.append(mbin[-1])

name_table = "table_merged_thetas_var_from_size.tsv"

table_path = str(out_dir) + '/' + name_table


with open(table_path, 'w+', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["t","value","value_nonorm","variance","trueN","meanBinDate","sampleSize","binningMode","numDaysSeqSpan"])
    for i in range(len(weeks)):
        # If bin size==1/variance is bigger or equal to min_bin_size
        if variance_sm[i]<=1/int(min_bin_size):# and num_days_wbin[i]>=min_days_span:
            writer.writerow([weeks[i],thetas_m[i]/(num_days_wbin[i]+1),thetas_m[i],variance_sm[i],rep_cases_wbin[i],dates_m[i],num_seqs_m[i],folders[i][:-3],num_days_wbin[i]])
        


"""
Plot the merged estimates
"""

plt.clf()
fig, ax1 = plt.subplots(figsize=(15,7))
ax1.set_title("Estimates from all binning modes")
ax1.set_xlabel('Mean bin date')
ax1.set_ylabel(r'$\theta_{est}$', color='crimson')
ax1.errorbar(dates_m, thetas_m, yerr=np.array(variance_sm)*10, fmt='.', elinewidth=0.7, color='crimson')
ax1.tick_params(axis='y', labelcolor='crimson')
ax1.set_xticks(dates_m[::3])
ax1.set_xticklabels(dates_m[::3],rotation=45)
ax2 = ax1.twinx()
ax2.set_ylabel('Active cases (%s-day window, poly-order %s)' % (window,poly), color='royalblue')
ax2.plot(dates_m, rep_cases_wbin, 'o', color='royalblue')
ax2.tick_params(axis='y', labelcolor='royalblue')

fig.tight_layout()

name = str(out_dir)+"/python_merged_thetas_smooth.png"
fig.savefig(str(name),dpi=300)

plt.clf()
fig, ax1 = plt.subplots(figsize=(15,7))
ax1.set_title("Estimates from all binning modes")
ax1.set_xlabel('Mean bin date')
ax1.set_ylabel(r'$\theta_{est}$', color='crimson')
ax1.errorbar(dates_m, thetas_m, yerr=np.array(variance_sm)*10, fmt='.', elinewidth=0.7, color='crimson')
ax1.tick_params(axis='y', labelcolor='crimson')
ax1.set_xticks(dates_m[::3])
ax1.set_xticklabels(dates_m[::3],rotation=45)
ax2 = ax1.twinx()
ax2.set_ylabel('Active cases on mean bin date', color='forestgreen')
ax2.plot(dates_m, totpos_onmean, 'o', color='forestgreen')
ax2.tick_params(axis='y', labelcolor='forestgreen')

fig.tight_layout()

name = str(out_dir)+"/python_merged_thetas_on_mean.png"
fig.savefig(str(name),dpi=300)
