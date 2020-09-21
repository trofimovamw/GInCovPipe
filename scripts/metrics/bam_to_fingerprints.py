#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:05:41 2020

@author: Maria Trofimova
"""
import pysam

#from Bio import SeqIO


class SAMtoFP:

    def __init__(self, filename):
        self.filename = filename


    def writeFP(self):
        # Load reference
#        ref_fasta = SeqIO.parse(open(self.reference),'fasta')
#        ref_seq = ''
#        for fasta in ref_fasta:
#            ref_seq = str(fasta.seq)

        # Write dictionary
        sequences_dict_orf = []
        sequences_dict_ = dict()
        cigar_dict = dict()

        # Write pairs
        sequences_list = []
        sequences_list_base = []
        sequences_list_int = []
        sequences_list_ = []
        cigar_list = []
        starts = []

        file = pysam.AlignmentFile(self.filename)

        for read in file.fetch("EMBOSS_001"):
            name = read.query_name
            seq = str(read.query_sequence)
            start_ = read.get_reference_positions()
            start = start_[0]
            starts.append(start_[0])
            # Trim with cigar string
            cigar = read.cigartuples
            cigarstr = read.cigarstring
            #print(cigar)
            #print(cigarstr)
            cigar_dict[name] = cigarstr
            cigar_list.append((name,start,cigarstr))
            mutants = []
            mutants_base = []
            mutantsint = []
            mutantsint_orf1ab = []
            # Counter relative to reference
            counter = start
            # Counter relative to query
            counter_q = 0

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
                        mutants.append(alt_rec)
                        mutants_base.append(alt_rec+'>'+str(alt_base))
                        # Only look at ORF1lab as an option
                        # 255..21544 - from annotation table
                        if counter>253.0 and counter<21545.0:
                            mutantsint_orf1ab.append(alt_rec)
                        mutantsint.append(counter)
                        counter += 1
                        counter_q += 1
                #else:
                #    counter += cigar[i][1]
            mutantsstr = ''
            mutantsstrbase = ''
            if mutants!=[]:
                mutantsstr = "-".join(mutants)
            if mutants_base!=[]:
                mutantsstrbase = "-".join(mutants_base)

            mutantsstr_orf1ab = ''
            if mutants!=[]:
                mutantsstr_orf1ab = "-".join(mutantsint_orf1ab)

            sequences_dict_orf.append((name,mutantsstr_orf1ab))
            sequences_dict_[name] = mutants
            sequences_list.append((name,mutantsstr))
            sequences_list_int.append((name,mutantsint))
            sequences_list_base.append((name,mutantsstrbase))
            mismatches = []

            mismatchesstr = ''
            if mismatches!=[]:
                mismatchesstr = "-".join(str(mismatches))
            sequences_list_.append((name,mismatchesstr))

        return sequences_dict_orf, sequences_dict_, cigar_dict, sequences_list, sequences_list_int, sequences_list_, cigar_list, starts, sequences_list_base
