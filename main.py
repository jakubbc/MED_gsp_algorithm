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
    """Check algorithms results and time for various parameter values and datasets"""
    config.data_ind = 0
    config.dname = config.data_sources[config.data_ind]

    config.minsup = 0.12
    config.wSize = 10
    config.minGap = 5
    config.maxGap = 20

    for d_ind in [0, 1, 2]:
        config.data_ind = d_ind
        config.dname = config.data_sources[config.data_ind]
        print(config.dname)
        for v in [0, 10, 15]:
            # change desired parameter
            config.maxGap = v
            print(f'parameter value: {v}')
            test_seqs, test_seqs_times = read_spmf_txt(config.dname)
            start = time.time()
            test_freq_seqs = gsp_algorithm(test_seqs, test_seqs_times, config.dname, config.minsup)
            print(f'Time: {round(time.time() - start, 1)}s')
            print(f'Number of sequences: {sum(len(v) for v in test_freq_seqs.values())}\n')
            # print(f'Sequences found with supports:\n{test_freq_seqs}')


# def performance_experiment():
#     """Check algorithms results and time for number of Leviathan data"""
#
#     seqs, seqs_times = read_spmf_txt(config.dname)
#     for
#     freq_seqs = gsp_algorithm(seqs, seqs_times, config.dname, config.minsup)
#     print(f'\nSequences found with supports:\n{freq_seqs}')
#     print(f'\nNumber of sequences: {sum(len(v) for v in freq_seqs.values())}')


if __name__ == "__main__":

    if config.conduct_unit_tests:
        test_generate_candidates()
        test_prune_candidates()
        test_seq_c_sup()

    if config.conduct_experiments:
        experiments()

    # if config.conduct_performance_experiment:
    #     performance_experiment()

    # seqs, seqs_times = load_test_data()
    seqs, seqs_times = read_spmf_txt(config.dname)
    start = time.time()
    freq_seqs = gsp_algorithm(seqs, seqs_times, config.dname, config.minsup)
    print(f'Time: {round(time.time() - start, 1)}s')
    # print(f'\nSequences found with supports:\n{freq_seqs}')
    print(f'\nNumber of sequences: {sum(len(v) for v in freq_seqs.values())}')
