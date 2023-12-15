#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 05:22:27 2023

@author: fred
"""

class CouldNotSolveError(Exception):
    "Raised if we cannot at least semi-robustly determine where the source of interest is in the field"
    pass

class PlateSolvedButNoFluxAtObject(Exception):
    "We could plate solve the image, but the acquisition image does not show a detection at the coordinates of the object"
    pass

class TooFewStarsForPlateSolving(Exception):
    "We have several stars, none of which seems to be much brighter than the others. Manual acquisition required"
    pass

class SuspiciousUniqueDetectionError(Exception):
    "We do have a detection, we could acquire, but it is very suspicious. Manual acquisition required."

class EmptyFieldError(Exception):
    "We could not detect stars in the field. Automatic acquisition not possible."
    
class NotWithinImageBoundError(Exception):
    "Thec catalogue coordinates do not fall in the footprint of the acquisition image."