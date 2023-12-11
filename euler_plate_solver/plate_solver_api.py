#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 06:17:28 2023

plate solves an image given a list of detections, using the astrometry.net API.

@author: fred
"""

from astroquery.astrometry_net import AstrometryNet
import requests
import json
from pathlib import Path
import os
import sys
from astropy.io import fits

from star_finder import extract_stars



sex = SExtractorRunner('sex', '.')
# fits_file_path = "../example_data/10h31m58s_-58d14m42s.fits"
# fits_file_path = "../example_data/plate_solve_try3.fits"
fits_file_path = '../example_data/9h02m_-70d56m.fits'
detections = sex.read_catalog('output.fits')
detections = sorted(detections, key=lambda x: -x['FLUX_APER'][0])
tupledetections = [(star['X_IMAGE'], star['Y_IMAGE']) for star in detections if star['FLAGS'] < 4][:10]

x, y = list(zip(*tupledetections))


import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
data = fits.getdata(fits_file_path)
plt.imshow(np.log10(data), origin='lower')
plt.plot(x, y, 'o', mfc='None', mew=2.5, ms=8, color='red')
plt.waitforbuttonpress()
plt.title('click on plot to continue')
#%%

# get an astrometry.net session:
R = requests.post('http://nova.astrometry.net/api/login', 
                  data={'request-json': json.dumps({"apikey": os.environ['astrometry_net_api_key']})})

# create an object to interact with our session:
ast = AstrometryNet()
ast._session_id = R.json()['session']

# load our ref image:
hh = fits.getheader(fits_file_path)
nx, ny = hh['NAXIS1'], hh['NAXIS2']

# aaand get the WCS. If not internet connection .......
# well you can still do it by hand.
wcs = ast.solve_from_source_list(x, y, nx, ny, 
                                 scale_est=0.255,
                                 scale_err=10,
                                 scale_units='arcsecperpix')
                                 # center_ra=157.991,
                                 # center_dec=-58.245,
                                 # radius=0.01)
if len(wcs) == 0:
    print('astrometry.net failed!!!!! wcs empty. try with different stars or a different image?')
    sys.exit()
    


#%%
with fits.open(fits_file_path, mode="update") as hdul:
    hdul[0].header.update(wcs)
    hdul.flush()  
