#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" config
Author: Jakub Ciemięga
"""

# logging flags
print_info_progress = 1
print_info = 0

# action flags
conduct_unit_tests = 0
conduct_experiments = 1

# datasets parameters
data_ind = 2
data_sources = ['bike', 'leviathan', 'retail']
dname = data_sources[data_ind]

# gsp parameters
minsup = 0.2
wSize = 0
minGap = 0
maxGap = float('inf')
