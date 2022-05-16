#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" gsp_algorithm
Author: Jakub Ciemięga
"""

import time


def gsp_algorithm(seqs: list, seqs_times: list, minsup: float = 0.2):
    """ Run

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
    print(minsup)
    freq_seqs = {}

    candidates = list({el for seq in seqs for tran in seq for el in tran})
    candidates.sort()
    # print(candidates)

    # count support
    start = time.time()
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
            freq_seqs[tuple([c])] = sup

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

    print(freq_seqs)
    # print(tid_list)
    print(time.time() - start)

# while

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
