#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 10:10:00 2019

@author: uli
"""

import ramics
from sys import argv

numpoints = ramics.parse_argv(argv)
if numpoints is None:
    print('\
Computes number of non-isomorphic posets with n points.\n\
Takes n as argument.\n\
Timing: n=3 : 0.05s\n\
        n=4 : 0.14s\n\
        n=5 : 62s\n\
        n=6 : ?? > 80min\n')
else:
    print(len(ramics.all_posets(numpoints)))
