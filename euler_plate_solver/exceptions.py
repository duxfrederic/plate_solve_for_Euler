#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 05:22:27 2023

@author: fred
"""


class CustomException(Exception):
    def __init__(self, message):
        super().__init__(message)
        self.message = message

class CouldNotSolveError(CustomException):
    def __init__(self, message="Raised if we cannot at least semi-robustly determine where the source of interest is in the field"):
        super().__init__(message)

class PlateSolvedButNoFluxAtObject(CustomException):
    def __init__(self, message="We could plate solve the image, but the acquisition image does not show a detection at the coordinates of the object"):
        super().__init__(message)

class TooFewStarsForPlateSolving(CustomException):
    def __init__(self, message="We have several stars, none of which seems to be much brighter than the others. Manual acquisition required"):
        super().__init__(message)

class SuspiciousUniqueDetectionError(CustomException):
    def __init__(self, message="We do have a detection, we could acquire, but it is very suspicious. Manual acquisition required."):
        super().__init__(message)

class EmptyFieldError(CustomException):
    def __init__(self, message="We could not detect stars in the field. Automatic acquisition not possible."):
        super().__init__(message)

class NotWithinImageBoundError(CustomException):
    def __init__(self, message="The catalogue coordinates do not fall in the footprint of the acquisition image."):
        super().__init__(message)


