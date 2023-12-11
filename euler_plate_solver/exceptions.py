#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 05:22:27 2023

@author: fred
"""

class CouldNotSolveError(Exception):
    "Raised if we cannot at least semi-robustly determine where the source of interest is in the field"
    pass