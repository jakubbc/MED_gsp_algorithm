#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" config
Author: Jakub Ciemięga
"""

# logging flags
print_info_progress = 1
print_info = 1

# action flags
conduct_unit_tests = 0
conduct_experiments = 0
conduct_performance_experiment = 1

# datasets parameters
data_ind = 1
data_sources = ['bike', 'leviathan', 'retail']
dname = data_sources[data_ind]

# gsp parameters
minsup = 0.15
wSize = 2
minGap = 2
maxGap = 10
