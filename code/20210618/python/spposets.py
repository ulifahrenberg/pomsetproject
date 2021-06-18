#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 10:24:30 2019

@author: uli
"""

import ramics
from sys import argv

numpoints = ramics.parse_argv(argv)
if numpoints is None:
    print('\
Computes number of non-isomorphic sp-posets with n points.\n\
Takes n as argument.\n\
Timing: n=3 : 0.06s\n\
        n=4 : 0.06s\n\
        n=5 : 0.09s\n\
        n=6 : 2.1s\n\
        n=7 : 150s\n')
else:
    print(len(ramics.sp_posets_memo(numpoints)))
