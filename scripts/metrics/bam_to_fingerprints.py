#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:05:41 2020

@author: Maria Trofimova
"""
import pysam


class SAMtoFP:

    def __init__(self, filename):
        self.filename = filename
    
    def _cigarToFP(self, seq, cigar, start, name):
        # Counter reference
        counter = start
        # Counter query
        counter_q = 0
        # Mutants list
        mutants_base = []

        for i in range(len(cigar)):
            if cigar[i][0]==0:
                counter += cigar[i][1]
                counter_q += cigar[i][1]
            elif cigar[i][0]==1:
                counter_q += cigar[i][1]
            elif cigar[i][0]==2:
                counter += cigar[i][1]
            elif cigar[i][0]==3:
                counter += cigar[i][1]
            elif cigar[i][0]==4:
                counter_q += cigar[i][1]
            elif  cigar[i][0]==7:
                counter += cigar[i][1]
                counter_q += cigar[i][1]
            # Record base that is a mismatch (X) on both query and reference counters
            elif cigar[i][0]==8:
                # Also write the mutation event in the string/fingerprint +ref_seq[counter]+'>'
                for j in range(cigar[i][1]):
                    alt_rec = str(counter)
                    alt_base = seq[counter_q]
                    if str(alt_base)!='N':
                        #mutants.append(alt_rec)
                        mutants_base.append(alt_rec+'>'+str(alt_base))
                    counter += 1
                    counter_q += 1
        
        mutantsstrbase = ''
        if mutants_base!=[]:
            mutantsstrbase = "-".join(mutants_base)
        
        return mutantsstrbase


    def writeFP(self):
        # Write pairs
        sequences_list_base = []
       
        file = pysam.AlignmentFile(self.filename)

        for read in file.fetch("EMBOSS_001"):
            name = read.query_name
            seq = str(read.query_sequence)
            start_ = read.get_reference_positions()
            start = start_[0]
            # Trim with cigar string
            cigar = read.cigartuples
            
            mutants_string = self._cigarToFP(seq, cigar, start, name)

            sequences_list_base.append((name,mutants_string))

        return sequences_list_base
