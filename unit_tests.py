#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" unit_tests
Author: Jakub Ciemięga
"""
import config
from gsp_algorithm import generate_candidates, prune_candidates, seq_c_sup


def test_generate_candidates():
    """Test function based on small datasets"""
    print(f'\nTest generate_candidates()')

    # test 1
    prev_seqs = [((3005, 1, 2, 3), 1000, 2000), ((1, 2), 3), ((1000, 1), (3030, 2),), (1, (3, 4)), (3042,),
                 (3069, (1, 2, 3),), (3,), (2, 3, 5)]
    good_candidates = [((1, 2), (3, 4)), ((1, 2), 3, 5), (3042, 3042), (3042, 3), ((3042, 3),), (3, 3042), (3, 3)]
    test_candidates = generate_candidates(prev_seqs)
    # print(test_candidates)
    # print(good_candidates)
    print(test_candidates == good_candidates)

    # test 2
    prev_seqs = [((1, 2), 3), ((1, 2), 4), (1, (3, 4)), ((1, 3), 5), (2, (3, 4)), (2, 3, 5)]
    good_candidates = [((1, 2), (3, 4)), ((1, 2), 3, 5)]
    test_candidates = generate_candidates(prev_seqs)
    print(test_candidates == good_candidates)

    # test 2
    prev_seqs = [(1,), (2,), (3,)]
    good_candidates = [(1, 1), (1, 2), ((1, 2),), (1, 3), ((1, 3),), (2, 1), (2, 2), (2, 3), ((2, 3),), (3, 1), (3, 2),
                       (3, 3)]
    test_candidates = generate_candidates(prev_seqs)
    print(test_candidates == good_candidates)


def test_prune_candidates():
    """Test function based on small datasets"""
    print(f'\nTest prune_candidates()')

    # test 1
    input_candidates = [((1, 2), (3, 4)), ((1, 2), 3, 5)]
    prev_seqs = [((1, 2), 3), ((1, 2), 4), (1, (3, 4)), ((1, 3), 5), (2, (3, 4)), (2, 3, 5)]
    good_candidates = [((1, 2), (3, 4))]
    test_candidates = prune_candidates(input_candidates, prev_seqs)
    # print(test_candidates)
    # print(good_candidates)
    print(test_candidates == good_candidates)

    # test 2
    input_candidates = [((1, 2), (3, 5)), ((1, 2, 3), 5)]
    prev_seqs = [((1, 2), 3), ((1, 2), 4), (1, (3, 4)), ((1, 3), 5), (2, (3, 4)), (2, 3, 5)]
    good_candidates = []
    test_candidates = prune_candidates(input_candidates, prev_seqs)
    print(test_candidates == good_candidates)

    # test 3
    input_candidates = [(1, 1), (1, 2), ((1, 2),), (1, 3), ((1, 3),), (2, 1), (2, 2), (2, 3), ((2, 3),), (3, 1), (3, 2),
                        (3, 3)]
    prev_seqs = [(1,), (2,), (3,)]
    good_candidates = [(1, 1), (1, 2), ((1, 2),), (1, 3), ((1, 3),), (2, 1), (2, 2), (2, 3), ((2, 3),), (3, 1), (3, 2),
                       (3, 3)]
    test_candidates = prune_candidates(input_candidates, prev_seqs)
    print(test_candidates == good_candidates)


def test_seq_c_sup():
    """Test function based on small datasets"""
    print(f'\nTest seq_c_sup()')
    # print('Warning: remember to set in config: wSize = 10, minGap = 5, maxGap = 20')
    old_params = config.wSize, config.minGap, config.maxGap
    config.wSize = 10
    config.minGap = 5
    config.maxGap = 20

    # test 1
    c = ((1, 2, 3), 1, 2, 4)
    seq_tids = {1: [10, 20, 21], 2: [15, 25, 30], 3: [5, 10, 15], 4: [45]}
    print(seq_c_sup(seq_tids, c) is True)

    # test 2
    c = ((1, 2, 3), 1, 2, 4)
    seq_tids = {1: [10, 20], 2: [15, 25, 30], 3: [5, 10, 15], 4: [45]}
    print(seq_c_sup(seq_tids, c) is False)

    # test 3
    c = ((1, 2, 3), 1, 2, 4)
    seq_tids = {1: [10, 20, 21], 2: [15, 25, 30], 3: [5, 10, 15], 4: [51]}
    print(seq_c_sup(seq_tids, c) is False)

    config.wSize, config.minGap, config.maxGap = old_params


