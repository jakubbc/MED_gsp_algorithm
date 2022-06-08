#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" gsp_algorithm
Author: Jakub Ciemięga
"""

import pickle
import itertools
from copy import copy
from operator import itemgetter

import numpy as np

import config


def common_element(list1: list, list2: list) -> bool:
    """ Return True if two lists have at least one common element and False otherwise
    """
    for x in list1:
        for y in list2:
            if x == y:
                return True
    return False


def generate_candidates(prev_seqs: list) -> list:
    """ generate new candidates of length k based on frequent sequences k-1 long

        :param prev_seqs: frequent sequences of k-1 length
        :type prev_seqs: list

        :return candidates: candidates of k length for frequent sequences
        :rtype: list
    """

    candidates = []
    seqs_no_first = {}
    seqs_no_last = {}
    for seq_n, seq in enumerate(prev_seqs):
        seq = list(seq)
        seqs_no_first[seq_n] = []
        seqs_no_last[seq_n] = []
        # remove first
        # for single item elements
        if type(seq[0]) is not tuple:
            seqs_no_first[seq_n].append(seq[1:])
        # for multi item elements
        else:
            # different action for len=2 bacause it has to be a tuple
            if len(seq[0]) == 2:
                for i in range(len(seq[0])):
                    seqs_no_first[seq_n].append([(seq[0][:i] + seq[0][i + 1:])[0]] + seq[1:])
            else:
                for i in range(len(seq[0])):
                    seqs_no_first[seq_n].append([seq[0][:i] + seq[0][i + 1:]] + seq[1:])
        # remove last
        # for single item elements
        if type(seq[-1]) is not tuple:
            seqs_no_last[seq_n].append(seq[:-1])
        # for multi item elements
        else:
            if len(seq[-1]) == 2:
                for i in range(len(seq[-1])):
                    seqs_no_last[seq_n].append(seq[:-1] + [(seq[-1][:i] + seq[-1][i + 1:])[0]])
            else:
                for i in range(len(seq[-1])):
                    seqs_no_last[seq_n].append(seq[:-1] + [seq[-1][:i] + seq[-1][i + 1:]])

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
                    candidates.append(tuple(list(prev_seqs[i]) + [prev_seqs[j][-1]]))
                # if last of seq2 is a tuple, add seq1 with last element taken from seq2
                else:
                    candidates.append(prev_seqs[i][:-1] + (prev_seqs[j][-1],))

    # get unique values
    seen = set()
    candidates = [x for x in candidates if not (x in seen or seen.add(x))]

    return candidates


def prune_candidates(candidates: list, prev_seqs: list) -> list:
    """ prune candidates by removing those who have contiguous subsequences that are not in freq_seqs of k-1 length

        :param candidates: candidates to prune
        :type candidates: list

        :param prev_seqs: frequent sequences of k-1 length
        :type prev_seqs: list

        :return candidates: candidates of k length for frequent sequences
        :rtype: list
    """
    candidates_pruned = []

    # for tests in development
    # candidates = [((1, 2), (3, 4)), ((1, 2), 3, 5)]
    # prev_seqs = [((1, 2), 3), ((1, 2), 4), (1, (3, 4)), ((1, 3), 5), (2, (3, 4)), (2, 3, 5)]

    # for each candidate if any contiguous subsequences not in prev_seqs, don't append candidate
    for c in candidates:
        continue_flag = False
        # if first element not tuple, remove it
        if type(c[0]) is not tuple:
            con_subseq = c[1:]
            if con_subseq not in prev_seqs:
                continue
        # if last element not tuple, remove it
        if type(c[-1]) is not tuple:
            con_subseq = c[:-1]
            if con_subseq not in prev_seqs:
                continue
        # for each tuple element, try removing each single item
        for ind, el in enumerate(c):
            if type(el) is tuple:
                for i in range(len(el)):
                    con_subseq = c[:ind] + (el[:i] + el[i + 1:]) + c[ind + 1:]
                    if con_subseq not in prev_seqs:
                        continue_flag = True
                        break
            if continue_flag:
                break
        if continue_flag:
            continue

        candidates_pruned.append(tuple(c))

    return candidates_pruned


def seq_c_sup(seq_tids: dict, c: tuple) -> bool:
    """ Check if given sequence (represented by a transaction time list for all elements) supports a candidate
    """

    def add_times(el_i: int):
        """Add list of tuples (start_time, end_time) to el_times for given el_i corresponding to given candidate c
        """
        el = c[el_i]
        # if not tuple, start_time = end_time
        if type(el) is not tuple:
            if el in seq_tids:
                for t in seq_tids[el]:
                    el_times[el_i].append((t, t))
        # for tuple elements, append all possible times
        elif all(it in seq_tids for it in el):
            # check all possible time combinations
            for times in itertools.product(*itemgetter(*el)(seq_tids)):
                # add those with appropriate wSize and are not duplicates
                if max(times) - min(times) <= config.wSize and (min(times), max(times)) not in el_times[el_i]:
                    el_times[el_i].append((min(times), max(times)))
        el_times_done[el_i] = True

    # for tests in development
    # c = ((1, 2, 3), 1, 2, 4)
    # seq_tids = {1: [10, 20, 21], 2: [15, 25, 30], 3: [5, 10, 15], 4: [45]}

    # list of start and end times of all elements of candidates
    el_times = [[] for _ in c]
    # to check if element times already appended
    el_times_done = [False for _ in c]

    # check first element
    add_times(0)
    if not el_times[0]:
        return False

    # start checking next elements in forward and backward phases
    el_ind = 1
    move_forward = True
    t1, t2 = el_times[0][0]
    # check while all elements checked, if not supported False returned immediately
    while el_ind < len(c):
        # generate start and end times if not already done
        if not el_times_done[el_ind]:
            add_times(el_ind)
        # if element not present (no times), c not supported
        if not el_times[el_ind]:
            return False
        # check next elements
        if move_forward:
            t1_new, t2_new = el_times[el_ind][0]
            if t1_new - t2 <= config.minGap:
                del el_times[el_ind][0]
                continue
            if t2_new - t1 > config.maxGap:
                del el_times[el_ind - 1][0]
                t1_next, t2_next = t1_new, t2_new
                move_forward = False
                el_ind -= 1
                continue
            t1, t2 = t1_new, t2_new
            el_ind += 1
        # backward phase - when maxGap exceeded
        else:
            t1, t2 = el_times[el_ind][0]
            if t2_next - t1 > config.maxGap:
                del el_times[el_ind][0]
                continue
            if 0 == el_ind:
                move_forward = True
                el_ind += 1
                continue
            t1_prev, t2_prev = el_times[el_ind - 1][0]
            # only check for maxGap, minGap checked in forward phase
            if t2 - t1_prev > config.maxGap:
                del el_times[el_ind - 1][0]
                el_ind -= 1
                continue
            move_forward = True
            el_ind += 1

    return True


def find_freq_seqs(candidates: list, tid_list: dict, minsup: int) -> dict:
    """ Find generalized frequent sequences among candidates based on wSize, minGap and maxGap

        :param candidates: candidates to find freq_seqs among
        :type candidates: list

        :param tid_list: dictionary of time occurrences of all items in all sequences. Used to speed up search
        :type tid_list: dict

        :param minsup: minimum number of occurrences in sequences to be considered frequent
        :type minsup: int

        :return freq_seqs: frequent sequences found in seqs
        :rtype: dict
    """
    if len(candidates) == 0:
        return {}
    new_freq_seqs = {}
    # create supports list
    sups = [0] * len(candidates)

    if config.print_info or config.print_info_progress:
        print(f'Sequences done of {len(tid_list)}: ', end='')
    # go through each transaction
    for sid in tid_list:
        if (config.print_info or config.print_info_progress) and sid % 1000 == 1:
            print(sid - 1, end=' ')
        for cid, c in enumerate(candidates):
            if seq_c_sup(tid_list[sid], c):
                sups[cid] += 1
    if config.print_info_progress and not config.print_info:
        print()
    if config.print_info:
        print(f'\nMaximum support found: {max(sups)}')
    for i in range(len(sups)):
        if sups[i] > minsup:
            new_freq_seqs[candidates[i]] = sups[i]
    return new_freq_seqs


def gsp_algorithm(seqs: list, seqs_times: list, dname, minsup: float = 0.2) -> dict:
    """ Run the gsp algorithm

        :param seqs: sequences to apply the algorithm to
        :type seqs: list

        :param seqs_times: sequences' transactions times
        :type seqs_times: list

        :param dname: dataset name for reading and saving help files with supports and tid_list
        :type minsup: str

        :param minsup: minimum sequence support defined as a fraction of sequence number
        :type minsup: float

        :return freq_seqs: frequent sequences found in seqs
        :rtype: dict
        """
    if config.print_info:
        print(f'\nMin support: {minsup} = ', end='')
    minsup *= len(seqs)
    if config.print_info:
        print(f'{minsup}')
        print(f'wSize: {config.wSize}')
        print(f'minGap: {config.minGap}')
        print(f'maxGap: {config.maxGap}')
    c_len = 1
    freq_seqs = {c_len: {}}

    candidates = list({el for seq in seqs for tran in seq for el in tran})
    candidates.sort()

    tid_list = {}
    # load from files and if no files, read raw data and create files
    try:
        with open(f'data/help/{dname}_sups.txt', 'r') as f:
            sups = [int(x) for x in f.read()[1:-1].split(', ')]
        with open(f'data/help/{dname}_tid_list.pkl', 'rb') as f:
            tid_list = pickle.load(f)
    except FileNotFoundError:
        print('Cannot find files with supports or tid_list, loading raw data files and saving help files.')
        sups = [0] * len(candidates)
        # count support
        for sid, seq in enumerate(seqs):
            tid_list[sid] = {}
            for tid, tran in enumerate(seq):
                for cid, c in enumerate(candidates):
                    if c in tran:
                        if c not in tid_list[sid]:
                            tid_list[sid][c] = [seqs_times[sid][tid]]
                            sups[cid] += 1
                        else:
                            tid_list[sid][c].append(seqs_times[sid][tid])

    # find frequent sequences of length 1
    finally:
        for cid, c in enumerate(candidates):
            if sups[cid] > minsup:
                freq_seqs[c_len][tuple([c])] = sups[cid]
        # save files
        with open(f'data/help/{dname}_sups.txt', 'w') as f:
            f.write(str(sups))
        with open(f'data/help/{dname}_tid_list.pkl', 'wb') as f:
            pickle.dump(tid_list, f)
    # for performance experiment
    # minsup *= 4
    # keys = list(tid_list.keys())
    # keys_len = len(keys)
    # for key in keys:
    #     tid_list[key + keys_len] = tid_list[key]
    # for key in keys:
    #     tid_list[key + 2 * keys_len] = tid_list[key]
    # for key in keys:
    #     tid_list[key + 3 * keys_len] = tid_list[key]

    # repeat as long as frequent sequences generated
    while len(freq_seqs[c_len]) > 0:
        c_len += 1
        if config.print_info:
            print(f'\nCurrent candidate len: {c_len}')
        candidates = generate_candidates(list(freq_seqs[c_len - 1].keys()))
        if config.print_info:
            print(f'Candidates: {len(candidates)}')
            # print(candidates)
        candidates_pruned = prune_candidates(candidates, list(freq_seqs[c_len - 1].keys()))
        if config.print_info:
            print(f'Candidates pruned: {len(candidates_pruned)}')
        freq_seqs[c_len] = find_freq_seqs(candidates_pruned, tid_list, minsup)
        if config.print_info:
            print(f'Found freq_seqs: {len(freq_seqs[c_len])}', end='')
            # print(freq_seqs[c_len])
            print()

    freq_seqs.popitem()
    return freq_seqs
