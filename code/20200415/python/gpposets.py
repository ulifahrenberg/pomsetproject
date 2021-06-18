#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 11:42:57 2019

@author: uli
"""

import ramics
from sys import argv

numpoints = ramics.parse_argv(argv)
if numpoints is None:
    print('\
Computes number of non-isomorphic gp-posets with n points.\n\
Takes n as argument.\n\
Timing: n=3 : 0.06s\n\
        n=4 : 0.4s\n\
        n=5 : 46s / 31s\n\
        n=6 : >2h\n\
	n=7 : ? / >9:40h\n')
else:
    print(len(ramics.gp_posets(numpoints)))
