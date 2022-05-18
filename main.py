#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Data Mining Methods -- Project
Jakub Ciemięga 2022
Warsaw University of Technology
"""
import time
from data_manager import read_spmf_txt, load_test_data
from gsp_algorithm import gsp_algorithm
import config
from unit_tests import test_generate_candidates, test_prune_candidates, test_seq_c_sup


def experiments():
    """Check algorithms results and performance for various parameter values and datasets"""
    config.data_ind = 0
    config.dname = config.data_sources[config.data_ind]

    config.minsup = 0.05
    config.wSize = 10
    config.minGap = 5
    config.maxGap = 20

    for d_ind in [0]:
        config.data_ind = d_ind
        config.dname = config.data_sources[config.data_ind]
        print(config.dname)
        for v in [10, 15]:
            config.maxGap = v
            print(v)
            seqs, seqs_times = read_spmf_txt(config.dname)
            start = time.time()
            freq_seqs = gsp_algorithm(seqs, seqs_times, config.dname, config.minsup)
            print(time.time() - start)
            print(f'\nNumber of sequences: {sum(len(v) for v in freq_seqs.values())}')
            # print(f'\nSequences found with supports:\n{freq_seqs}')


if __name__ == "__main__":

    if config.conduct_unit_tests:
        test_generate_candidates()
        test_prune_candidates()
        test_seq_c_sup()

    # seqs, seqs_times = load_test_data()
    seqs, seqs_times = read_spmf_txt(config.dname)
    freq_seqs = gsp_algorithm(seqs, seqs_times, config.dname, config.minsup)
    print(f'\nSequences found with supports:\n{freq_seqs}')
    print(f'\nNumber of sequences: {sum(len(v) for v in freq_seqs.values())}')
    
    # experiments()

