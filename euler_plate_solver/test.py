#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 07:09:33 2023

@author: fred
"""

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from astropy.visualization import ZScaleInterval

import matplotlib.pyplot as plt

fits_file_path = "../example_data/plate_solve_try3.fits"

image = fits.getdata(fits_file_path).astype(float)
mean, median, stddev = sigma_clipped_stats(image)
image -= median

# Create a DAOStarFinder object to detect sources in the image
daofind = DAOStarFinder(fwhm=5.0, threshold=15*stddev)


# Detect sources in the image
sources = daofind(image)

# Plot the image with detected sources marked
plt.figure(figsize=(15, 15))
interval = ZScaleInterval(contrast=0.1)
vmin, vmax = interval.get_limits(image)
plt.imshow(image, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
plt.colorbar()
plt.plot(sources['xcentroid'], sources['ycentroid'], 'o',
            color='red', label='Detected Sources',
            mfc='None')
plt.legend()
plt.show()
