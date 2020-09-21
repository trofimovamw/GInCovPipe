#!/usr/bin/env python
'''
Binning factory
Creates bins of equal size, (equal number of days and equal number of days (fuzzy))
'''
from pathlib import Path
import sys

from sam_to_bins_modular import SAM
from generate_matrices import binnedBSAM
from sam_to_bins_modular_test import testDIR


FILEPATH = Path(__file__).parent
s = str(sys.argv[1])
SAM_PATH = s # Main big SAM or BAM file
print(SAM_PATH)
WORKING_PATH = FILEPATH.parent

sam = SAM(SAM_PATH)
#print('\n')
#print('-'*80)
#print("Bins of equal Size")
#print('-'*80)
#sam.bin_eq_size_names()
print('\n')
print('-'*80)
print("Bins with equal number of days")
print('-'*80)
sam.bin_eq_days()
print('\n')
print('-'*80)
print("Bins with equal number of days (fuzzy)")
print('-'*80)
sam.bin_eq_days(fuzzy=True)
print('\n')
print('-'*80)
print("Bins by calendar weeks")
print('-'*80)
sam.bin_cal_week()

#SAM_BINDIR1 = FILEPATH.parent / "eq_size"
#sam_dir = binnedBSAM(SAM_BINDIR1)
#sam_dir.count_matrices()

ANNOT_DIR = FILEPATH.parent / "reference_prep/annotation.txt"

SAM_BINDIR2 = FILEPATH.parent / "eq_days"
sam_dir = binnedBSAM(SAM_BINDIR2,ANNOT_DIR)
sam_dir.count_matrices()
sam_dir.truncate_matrix()

SAM_BINDIR3 = FILEPATH.parent / "fuzzy_days"
sam_dir = binnedBSAM(SAM_BINDIR3,ANNOT_DIR)
sam_dir.count_matrices()
sam_dir.truncate_matrix()

SAM_BINDIR4 = FILEPATH.parent / "cal_week"
sam_dir = binnedBSAM(SAM_BINDIR4,ANNOT_DIR)
sam_dir.count_matrices()
sam_dir.truncate_matrix()

#test = testDIR(SAM_BINDIR2)
#test.test_binning()

test = testDIR(SAM_BINDIR2)
print('\n')
print('-'*80)
print("Equal days test")
print('-'*80)
test.test_binning()

test = testDIR(SAM_BINDIR3)
print('\n')
print('-'*80)
print("Fuzzy days test")
print('-'*80)
test.test_binning()

test = testDIR(SAM_BINDIR4)
print('\n')
print('-'*80)
print("Calendar week test")
print('-'*80)
test.test_binning()



