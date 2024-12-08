# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 08:04:20 2023

@author: fred
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from pathlib import Path

from coralie_plate_solver.star_finder import extract_stars, plot_stars



base_path = Path('example_data')

file_names = [
    # '9h02m_-70d56m.fits',
    # '10h31m58.0_-58:14:42.fits',
    '6h45m53.0_10:03:48.fits',
    '8h44m12.0_-35:28:24.fits',
    '3h07m18.6_-13:45:42.fits',
    '7h40m32.8_2:05:55.fits',
    '07h14m51.5_-29:25:50.fits',
    '9h55m26.2_-58:25:47.fits',
    '9h02m23.3_-70:56:51.fits',
    '10h31m58.0_-58:14:42_BIS.fits',
    '11h12m10.1_-61:45:18.fits',
    '9h40m41.9_-66:39:16.fits',
    '2h54m35.0_-49:28:04.fits',
    '5h48m17.8_-15:11:58.fits',
    '6h45m53.0_10:03:48.fits',
    '8h32m33.2_-23:23:05.fits'
]

fits_file_paths = [base_path / file_name for file_name in file_names]

savedir = Path('_output_test_extract_stars')
savedir.mkdir(exist_ok=True)


for fits_file_path in fits_file_paths:
    
    fits_file_path = Path(fits_file_path)
    
    sources, image_sub = extract_stars(fits_file_path)
    
    savepath = savedir / f"{fits_file_path.name}.jpeg"
    plot_stars(sources, image_sub, savepath)        