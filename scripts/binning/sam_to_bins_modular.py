#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 17:19:40 2020

@authors: Yannick Duport, Maria Trofimova

Getting base counts with pysam and binning routine example
"""
import csv
import datetime
from datetime import date
import math
import numpy as np
import pandas as pd
import os
import pysam
import re
from warnings import warn
strptime = datetime.datetime.strptime


class SAM:
    def __init__(self, samfile_path, bins_dir, num_per_bin, days_per_bin, seq_name):
        self.samfile_path = samfile_path
        self.seq_name = seq_name
        self.samfile = pysam.AlignmentFile(self.samfile_path, "rb")
        self.header = self.samfile.header
        self.sam_dict = self._index_dates()
        self.eq_num_param = num_per_bin
        self.eq_days_param = days_per_bin
        self.bins_dir = bins_dir

    def _to_dict(self, l):
        sam_dict = {}
        for i in range(len(l)):
            idx, date, name = l[i]
            sam_dict[i] = {
                'index': idx,
                'date': date,
                'header': name
            }
        return sam_dict
    
    def _index_dates(self):
        """Extracts headers, dates and index from a samfile
        Extracts headers, dates, index from samfile, sorts them by date and adds them to a dictionary
        Dictionary keys contain the index from the sorted dates, values are header, date, index
        :return: dict[sorted_index]={header:string, date:string, index:int}
        """
        idx_dates = []
        index = 0
        for read in self.samfile.fetch(self.seq_name):
            name = read.query_name
            # date_format = "[|]20[\d]{2}\-[\d]{2}\-[\d]{2}"
            date_format = r"[|]20[\d]{2}\-*\d*-*\d*"
            dt = re.search(date_format, name)
            if dt:
                dt = dt[0][1:]
                if len(dt) == 7:    # dates without days (format YYYY-mm)
                    # FIXME: include date in case for certain binnings (e.g. by month)
                    # dt = dt+"-15"
                    continue
                elif len(dt) == 4:  # dates without months and years (format YYYY)
                    # dt = f"{dt}-02-15"
                    continue
                try:
                    Date = strptime(dt, "%Y-%m-%d").date()
                except ValueError:
                    Date = strptime(dt, "%Y-%d-%m").date()
                    dt = f"{dt[0:5]}{dt[8:]}{dt[4:7]}"
                    warn(f"Date-format different from '%Y-%m-%d' in: {name}\n"
                         f"'%Y-%d-%m' used instead")
                Today = date.today()
                if Date > Today:
                    dt = f"{dt[0:5]}{dt[8:]}{dt[4:7]}"
                idx_dates.append((index, dt, name))
                index += 1
            else:
                warn(f"No date found in the following query name: {name}")
        idx_dates.sort(key=lambda x: x[1])
        sam_dict = self._to_dict(idx_dates)
        return sam_dict

    def _create_filenames(self, folder, n_bins, infix=None):
        """Creates folder, files and returns filenames

        :param method: string
        :param n_bins: int
        :return:
        """
        filepath = self.bins_dir / folder
        basename = os.path.basename(self.samfile_path)
        basename = os.path.splitext(basename)[0]

        os.makedirs(filepath, exist_ok=True)
        filenames_bin = sorted([filepath / f"bin_{basename}_{idx:04d}.bam" for idx in range(n_bins)])
        filenames_header = sorted([filepath / f"header_{basename}_{idx:04d}.tsv" for idx in range(n_bins)])
        filenames_range = sorted([filepath / f"range_{basename}_{idx:04d}.tsv" for idx in range(n_bins)])
        
        return filenames_bin, filenames_header, filenames_range

    def _write_bins(self, filenames, indices):
        bin_idx = 0
        for filename in filenames:
            print(f"Writing Bin {bin_idx}")
            names = [self.sam_dict[i]['header'] for i in indices[bin_idx]]
            with pysam.AlignmentFile(filename, "wb", header=self.header) as outfile:
                for read in self.samfile.fetch(self.seq_name):
                    # https://github.com/pysam-developers/pysam/issues/509!!!
                    read.tid = 0
                    if read.query_name in names:
                        outfile.write(read)
            bin_idx = bin_idx + 1

    def _write_header(self, filenames, indices):
        # uniquify header names
        for i in range(len(indices)):
            headers = []
            idx_list = []
            for idx in indices[i]:
                if self.sam_dict[idx]['header'] not in headers:
                   headers.append(self.sam_dict[idx]['header'])
                   idx_list.append(idx)
            indices[i] = idx_list
        fieldnames = ['header', 'date']
        for i in range(len(filenames)):
            with open(filenames[i], 'w+') as csvfile:
                writer = csv.DictWriter(csvfile,
                                        delimiter='\t',
                                        fieldnames=fieldnames,
                                        extrasaction='ignore')
                writer.writeheader()
                for j in indices[i]:
                    writer.writerow(self.sam_dict[j])
                    
    def _write_days_ranges(self, filenames, days_ranges):
        for i, filename in enumerate(filenames):
            with open(filename, 'w+') as csvfile:
                writer = csv.writer(csvfile,
                                    delimiter='\t')
                writer.writerow(['start_day','end_day'])
                writer.writerow([days_ranges[i][0],days_ranges[i][1]])
    

    def bin_eq_size_names(self):
        """Creates bins of equal size (n=10), but files with same name are treated as one

        :return:
        """
        print(f"Reads per bin: {self.eq_num_param}")

        binsize = self.eq_num_param
        n_reads = len(self.sam_dict)
        names = [self.sam_dict[i]['header'] for i in range(n_reads)]
        names_unique = pd.unique(names)
        n_bins = math.ceil(len(names_unique)/binsize)

        filenames_bin, filenames_header, filenames_range  = self._create_filenames(f"eq_size_names_{binsize}", n_bins)
        for filename in filenames_bin:
            open(filename, 'w+')
        
        bins_n = [[] for _ in range(n_bins)]
        indices = [[] for _ in range(n_bins)]
        for i in range(len(names_unique)):
            bins_n[math.floor(i / binsize)].append(names_unique[i])
        
        for i in range(len(bins_n)):
            for ii in range(len(bins_n[i])):
                for j in range(n_reads):
                    if self.sam_dict[j]['header'] == bins_n[i][ii]:
                        indices[i].append(j)

        self._write_bins(filenames_bin, indices)
        self._write_header(filenames_header, indices)

    def bin_eq_days(self, fuzzy=False):
        """Create bins containing equal number of days

        :param fuzzy:
        :return:
        """
        if fuzzy:
            folder = f"fuzzy_days"
        else:
            folder = f"eq_days"

        n_reads = len(self.sam_dict)
        date_format = "%Y-%m-%d"
        start_day = strptime(self.sam_dict[0]['date'], date_format).date()
        end_day = strptime(self.sam_dict[n_reads - 1]['date'], date_format).date()
        n_days = (end_day - start_day).days + 1
        days_per_bin = self.eq_days_param
        n_bins = math.ceil(n_days / days_per_bin)
        
        days_ranges = []
        
        for i in range(n_bins):
            days_ranges.append((start_day+datetime.timedelta(days=(days_per_bin*i)),start_day+datetime.timedelta(days=(days_per_bin*(i+1)-1))))
            

        print(f"total number of days: {n_days}")
        print(f"days per bin: {days_per_bin}")

        # Create empty files
        filenames_bin, filenames_header, filenames_range = self._create_filenames(f"{folder}_{days_per_bin}", n_bins)
        for filename in filenames_bin:
            open(filename, 'w+')

        # fill list of indices and bins
        indices = [[] for _ in range(n_bins)]
        curr_bin = 0
        for i in range(n_reads):
            curr_day = strptime(self.sam_dict[i]['date'], date_format).date()
            curr_n_days = (curr_day - start_day).days + 1
            if curr_n_days > days_per_bin:
                # r1 = np.random.uniform(low=0.0, high=1.0, size=1)
                # r2 = np.random.uniform(low=0.0, high=1.0, size=1)
                # if fuzzy and (r1 > r2):
                if fuzzy and np.random.randint(2):
                    pass
                else:
                    curr_bin += 1
                    start_day = strptime(self.sam_dict[i]['date'], date_format).date()
            indices[curr_bin].append(i)

        # write bam-files (1 file per bin)
        self._write_bins(filenames_bin, indices)
        self._write_header(filenames_header, indices)
        self._write_days_ranges(filenames_range, days_ranges)
        
 
    def bin_cal_week(self):
        """Binning by calendar week

        :return:
        """
        n_reads = len(self.sam_dict)
        date_format = "%Y-%m-%d"
        start_day = strptime(self.sam_dict[0]['date'], date_format).date()
        end_day = strptime(self.sam_dict[n_reads - 1]['date'], date_format).date()
        start_day_offset = start_day.isocalendar()[2]
        n_days = (end_day - start_day).days
        n_bins = math.ceil((n_days + start_day_offset) / 7)
        
        days_ranges = []
        
        for i in range(n_bins):
            days_ranges.append((start_day+datetime.timedelta(days=(7*i)),start_day+datetime.timedelta(days=(7*(i+1)-1))))
            

        filenames_bin, filenames_header, filenames_range = self._create_filenames('cal_week', n_bins)
        for filename in filenames_bin:
            open(filename, 'w+')

        indices = [[] for _ in range(n_bins)]
        counter = 0
        curr_year, curr_week = strptime(self.sam_dict[0]['date'], date_format).isocalendar()[:2]
        for i in range(n_reads):
            next_year, next_week = strptime(self.sam_dict[i]['date'], date_format).isocalendar()[:2]
            if next_week > curr_week or next_year > curr_year:
                counter += 1
                curr_year, curr_week = next_year, next_week
            indices[counter].append(i)

        # write bam-files (1 file per bin)
        self._write_bins(filenames_bin, indices)
        self._write_header(filenames_header, indices)
        self._write_days_ranges(filenames_range, days_ranges)
        

"""Deprecated

    def bin_eq_size(self, binsize=10):
        '''Creates bins of equal size (n=10)
        :return:
        '''
        print(f"Reads per bin: {binsize}")
        n_reads = len(self.sam_dict)
        n_bins = math.ceil(n_reads/binsize)
        filenames_bin, filenames_header = self._create_filenames(f"eq_size_{binsize}", n_bins)
        for i in range(n_bins):
            open(filenames_bin[i], 'w+')

        indices = [[] for _ in range(n_bins)]
        for i, value in self.sam_dict.items():
            indices[math.floor(i / binsize)].append(i)
    
        self._write_bins(filenames_bin, indices)
        self._write_header(filenames_header, indices)
"""