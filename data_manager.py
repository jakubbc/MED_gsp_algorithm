#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" data_manager
Author: Jakub Ciemięga
"""


def read_spmf_txt(fname: str, load_times: bool = False) -> list:
    """ Read txt file in an SPMF format

        :param fname: name of file to read from
        :type fname: str

        :param load_times: decides, if additional file with transaction times is provided by user. If not, generates
                            transaction times by transaction order
        :type load_times: bool

        :return (seqs, seq_times):
            seqs: list of sequences from file: seqs=[seq1, seq2, ...], seq1=[tran1, tran2, ...], tran1=[el1, el2, ...]
            seq_times: list of transaction times corresponding to seqs: seq_times=[seq1, seq2, ...], seq1=[tran_time1,
                        tran_time2, ...]
        :rtype: (list, list)
        """

    seqs = []
    seqs_times = []

    with open(f'data/{fname}.txt') as f:
        line = f.readline()
        while line:
            seq = line.split(' -1 ')
            seq.pop()
            seq = [tran.split(' ') for tran in seq]
            seq = [[int(el) for el in tran] for tran in seq]
            seqs.append(seq)

            if not load_times:
                seqs_times.append(list(range(1, len(seq) + 1)))

            line = f.readline()

    if load_times:
        with open(f'data/{fname}_times.txt') as f:
            line = f.readline()
            while line:
                seq_times = line.split(' -1 ')
                seq_times.pop()
                seq_times = [int(time) for time in seq_times]
                seqs_times.append(seq_times)
                line = f.readline()

    return seqs, seqs_times


def load_test_data() -> list:
    """ Generate simple data for easier testing
    :return (seqs, seq_times):
            seqs: list of sequences
            seq_times: list of transaction times corresponding to seqs
        :rtype: (list, list)
    """

    seqs = [
            [[1], [2]],
            [[1], [3], [4], [5], [1]],
            [[2], [3], [4], [6]],
            [[1], [2], [3], [4]],
            [[1], [2], [3], [6]]
        ]
    seqs_times = [
            [1, 2],
            [1, 2, 3, 4, 5],
            [1, 2, 3, 4],
            [1, 2, 3, 4],
            [1, 2, 3, 4]
        ]

    return seqs, seqs_times
