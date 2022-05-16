#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Data Mining Methods -- Project
Jakub Ciemięga 2022
Warsaw University of Technology
"""

from data_manager import read_spmf_txt, load_test_data
from gsp_algorithm import gsp_algorithm
from gsp_test import TestGSP

data_sources = ['bike', 'retail', 'leviathan']
data_ind = 1

if __name__ == "__main__":

    # seqs, seqs_times = load_test_data()
    seqs, seqs_times = read_spmf_txt(data_sources[data_ind])
    # print(seqs[0:10])
    # print(seqs_times[0:10])
    gsp_algorithm(seqs, seqs_times)
    # TestGSP().test_supermarket()
