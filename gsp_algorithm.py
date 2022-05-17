#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" gsp_algorithm
Author: Jakub Ciemięga
"""

import time
import numpy as np


def common_element(list1: list, list2: list) -> bool:
    """ Return True if two lists have at least one common element and False otherwise
    """
    for x in list1:
        for y in list2:
            if x == y:
                return True
    return False


def generate_candidates(prev_seqs: dict) -> list:
    """ generate new candidates of length k based on frequent sequences k-1 long

        :param prev_seqs: frequent sequences of k-1 length
        :type prev_seqs: dict

        :return candidates: frequent sequences of k length
        :rtype: list
        """
    # prev_seqs = [((3005, 1, 2, 3), 1000, 2000), ((1, 2), 3), ((1000, 1), (3030, 2),), (1, (3, 4)), (3042,),
    #              (3069, (1, 2, 3),), (3,), (2, 3, 5)]
    print(prev_seqs)
    candidates = []
    seqs_no_first = {}
    seqs_no_last = {}
    for seq_n, seq in enumerate(prev_seqs):
        seq = list(seq)
        seqs_no_first[seq_n] = []
        seqs_no_last[seq_n] = []
        # remove first
        if type(seq[0]) is not tuple:
            seqs_no_first[seq_n].append(seq[1:])
        else:
            if len(seq[0]) == 2:
                for i in range(len(seq[0])):
                    seqs_no_first[seq_n].append([(seq[0][:i] + seq[0][i + 1:])[0]] + seq[1:])
            else:
                for i in range(len(seq[0])):
                    seqs_no_first[seq_n].append([seq[0][:i] + seq[0][i + 1:]] + seq[1:])
        # remove last
        if type(seq[-1]) is not tuple:
            seqs_no_last[seq_n].append(seq[:-1])
        else:
            if len(seq[-1]) == 2:
                for i in range(len(seq[-1])):
                    seqs_no_last[seq_n].append(seq[:-1] + [(seq[-1][:i] + seq[-1][i + 1:])[0]])
            else:
                for i in range(len(seq[-1])):
                    seqs_no_last[seq_n].append(seq[:-1] + [seq[-1][:i] + seq[-1][i + 1:]])

    print(f'{seqs_no_first=}')
    print(f'{seqs_no_last=}')

    # join phase -- join all seqs where no_fist==no_last by adding last symbol of s2 to s1
    for i in range(len(prev_seqs)):
        for j in range(len(prev_seqs)):
            if common_element(seqs_no_first[i], seqs_no_last[j]):
                # if sequences len = 1, add both as separate item and as itemset
                if len(prev_seqs[i]) == 1:
                    candidates.append((prev_seqs[i][0], prev_seqs[j][0]))
                    if i < j:
                        candidates.append(((prev_seqs[i][0], prev_seqs[j][0]),))
                # if last of seq2 is a single element, add it to seq1
                elif type(prev_seqs[j][-1]) is not tuple:
                    candidates.append(tuple(list(prev_seqs[i])+[prev_seqs[j][-1]]))
                else:
                    candidates.append(prev_seqs[i][:-1]+(prev_seqs[j][-1],))

    print(candidates)

    return []


def gsp_algorithm(seqs: list, seqs_times: list, minsup: float = 0.2):
    """ Run the gsp algorithm

        :param seqs: sequences to apply the algorithm to
        :type seqs: list

        :param seqs_times: sequences' transactions times
        :type seqs_times: list

        :param minsup: minimum sequence support defined as a fraction of sequence number
        :type minsup: float

        :return freq_seqs: frequent sequences found in seqs
        :rtype: # TODO add type
        """

    minsup *= len(seqs)
    # print(minsup)
    c_len = 1
    freq_seqs = {c_len: {}}

    candidates = list({el for seq in seqs for tran in seq for el in tran})
    candidates.sort()
    # print(candidates)

    # count support
    # start = time.time()
    tid_list = {}
    for sid, seq in enumerate(seqs):
        for tid, tran in enumerate(seq):
            for c in candidates:
                if c in tran:
                    if c not in tid_list:
                        tid_list[c] = {sid: [tid]}
                    elif sid not in tid_list[c]:
                        tid_list[c][sid] = [tid]
                    else:
                        tid_list[c][sid].append(tid)
    for c in candidates:
        sup = len(tid_list[c])
        if sup > minsup:
            freq_seqs[c_len][tuple([c])] = sup

    # old version, multiple passes of data set
    # for c in candidates:
    #     count = 0
    #     for seq in seqs:
    #         for tran in seq:
    #             if c in tran:
    #                 count += 1
    #                 break
    # count += tran.count(c)
    # if count > minsup:
    #     freq_seqs[tuple([c])] = count
    # # print(count)

    # print(freq_seqs)
    # print(tid_list)
    # print(time.time() - start)

    # repeat as long as frequent sequences generated
    while len(freq_seqs[c_len]) > 0:
        c_len += 1
        candidates = generate_candidates(list(freq_seqs[c_len - 1].keys()))
        freq_seqs[c_len] = candidates

# # the set of frequent 1-sequence: all singleton sequences
# # (k-itemsets/k-sequence = 1) - Initially, every item in DB is a
# # candidate
# candidates = self.unique_candidates
#
# # scan transactions to collect support count for each candidate
# # sequence & filter
# self.freq_patterns.append(self._support(candidates, minsup))
#
# # (k-itemsets/k-sequence = 1)
# k_items = 1
#
# self._print_status(k_items, candidates)
#
# # repeat until no frequent sequence or no candidate can be found
# while len(self.freq_patterns[k_items - 1]) and (k_items + 1 <= self.max_size):
#     k_items += 1
#
#     # Generate candidate sets Ck (set of candidate k-sequences) -
#     # generate new candidates from the last "best" candidates filtered
#     # by minimum support
#     items = np.unique(
#         list(set(self.freq_patterns[k_items - 2].keys())))
#
#     candidates = list(product(items, repeat=k_items))
#
#     # candidate pruning - eliminates candidates who are not potentially
#     # frequent (using support as threshold)
#     self.freq_patterns.append(self._support(candidates, minsup))
#
#     self._print_status(k_items, candidates)
# return self.freq_patterns[:-1]
