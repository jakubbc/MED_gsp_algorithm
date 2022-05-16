#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Data Mining Methods -- Project
Jakub Ciemięga 2022
Warsaw University of Technology
"""

from data_manager import read_spmf_txt

data_sources = ['bike', 'retail', 'leviathan']
data_ind = 0

if __name__ == "__main__":

    seqs, seqs_times = read_spmf_txt(data_sources[data_ind])
    print(seqs[0])
    print(seqs_times[0])
