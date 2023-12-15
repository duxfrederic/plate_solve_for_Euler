#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 08:04:20 2023

@author: fred
"""

from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
import requests
import json
import os
from astroquery.astrometry_net import AstrometryNet
import logging
logger = logging.getLogger(__name__)


from euler_plate_solver.star_finder import extract_stars, plot_stars
from euler_plate_solver.exceptions import CouldNotSolveError


N_BRIGHTEST_TO_USE = 15
REDO = False

def plate_solve_with_API(fits_file_path, sources, ra_approx=None, dec_approx=None):
    """
    Calculate the WCS using the nova.astrometry.net API.
    In this case, we first extract the sources and only send those over
    the internet.
    We get back an astrometric solution contained in a fits file.

    Parameters:
    fits_file_path (str): Path to the FITS file.
    sources (astropy.table.Table): Table of detected sources.
    ra_approx: float in degrees, approximate center of the field coord
    dec_approx: float in degrees, approximate center of the field coord

    Returns:
    WCS header if successful, None otherwise.

    """
    logger.info(f"plate_solve_with_API on {fits_file_path}")
    # check, maybe we already solved it
    header = fits.getheader(fits_file_path)
    if 'PL-SLVED' in header and not REDO:
        logger.info(f"{fits_file_path} was already plate solved, no redoing.")
        # we already treated it.
        # the header contains the WCS info.
        return header

    # unpack the positions
    tupledetections = [(star['xcentroid'], star['ycentroid']) for star in sources[:N_BRIGHTEST_TO_USE]]
    x, y = list(zip(*tupledetections))

    # create a session
    R = requests.post('http://nova.astrometry.net/api/login',
                      data={'request-json': json.dumps({"apikey": os.environ['astrometry_net_api_key']})})

    ast = AstrometryNet()
    ast._session_id = R.json()['session']

    # the dimensions of the image...
    nx, ny = header['NAXIS1'], header['NAXIS2']

    # can we make the solver find the solution faster?
    # let's give it the pointing of the telescope.
    morekwargs = {}
    if ra_approx is not None and dec_approx is not None:
        morekwargs['center_ra'] = ra_approx
        morekwargs['center_dec'] = dec_approx
        morekwargs['radius'] = 1. # 1 degree ...if it's worse than that, tchao.
    try:
        # TODO FACTOR OUT PIXEL SCALE ESTIMATE
        # HERE I HARDCODED THE PIXEL SCALE OF CORALIE'S BIGEYE
        wcs = ast.solve_from_source_list(x, y, nx, ny,
                                         scale_est=0.26,
                                         scale_err=10,
                                         scale_units='arcsecperpix',
                                         publicly_visible='n',
                                         **morekwargs)
    except Exception as e:
        # anything ...
        logger.info(f"plate_solve_with_API: something went wrong with API when trying to solve {fits_file_path}")
        raise CouldNotSolveError(f'Exception when trying to plate solve the image: {e}')

    # it can also have actually failed.
    # happens if the selected stars are too faint, or if our sources
    # contain too many artifacts.
    # or if we don't have enough stars
    # or ...many other reasons.
    if len(wcs) == 0:
        logger.info(f"plate_solve_with_API: could not plate solve {fits_file_path}")

        raise CouldNotSolveError('Astrometry.net failed! WCS empty. Try with different stars or a different image?')

    # else, we're probably fine. Like, 99.9% confidence from my
    # experience with astrometry.net's plate solver.
    logger.info(f"plate_solve_with_API: {fits_file_path} solved, writing the WCS")
    with fits.open(fits_file_path, mode="update") as hdul:
        # little flag to indicate that we plate solved this file:
        wcs['PL-SLVED'] = 'done'
        # add all this info to the file:
        hdul[0].header.update(wcs)
        hdul.flush()
    return wcs




def plot_stars_with_wcs(fits_file_path, sources):
    """
    Plot the updated image with WCS projection and marked detected sources.
    This is to be used for diagnostics.

    Parameters:
    fits_file_path (str): Path to the updated FITS file with WCS.
    sources (astropy.table.Table): Table of detected sources.
    """
    with fits.open(fits_file_path) as hdul:
        wcs = WCS(hdul[0].header)
        image = hdul[0].data

    fig = plt.figure(figsize=(15, 15))
    ax = plt.subplot(projection=wcs)
    interval = ZScaleInterval(contrast=0.1)
    vmin, vmax = interval.get_limits(image)

    ax.imshow(image, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
    ax.grid(color='white', ls='solid')
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')

    ra, dec = wcs.all_pix2world(sources['xcentroid'], sources['ycentroid'], 0)
    ax.plot(ra, dec, 'o', color='red', label='Detected Sources', mfc='None',
            transform=ax.get_transform('world'))

    plt.legend()
    plt.show()



# if __name__ == '__main__':
#     # fits_file_path = "../example_data/plate_solve_try3.fits"
#     # fits_file_path = "../example_data/10h31m58s_-58d14m42s.fits"
#     # fits_file_path = "../example_data/8h44m12-35d28m24.fits"
#     # fits_file_path = "../example_data/8h44m12.0_-35:28:24.fits"

#     # fits_file_path = "../example_data/plate_solve_try2.fits"
#     # fits_file_path = "/home/fred/Documents/plate_solve_for_Euler/example_data/11h16m06-30d10m40.fits"
#     # fits_file_path = "../example_data/ECAM.2023-12-11T00:42:12.000.fits"
#     # fits_file_path = "../example_data/ECAM.2023-12-11T01:14:25.000.fits"
#     # fits_file_path = "../example_data/5h48m17_-15:11:58.fits"
#     # fits_file_path = "../example_data/9h02m_-70d56m.fits"
#     fits_file_path = '../example_data/10h31m58.0_-58:14:42.fits'
#     # fits_file_path = '../example_data/6h45m53.0_10:03:48.fits'
#     # fits_file_path = '../example_data/8h44m12.0_-35:28:24.fits'
#     # fits_file_path = '../example_data/3h07m18.6_-13:45:42.fits'
#     # fits_file_path = '../example_data/7h40m32.8_2:05:55.fits'
#     fits_file_path = '../example_data/07h14m51.5_-29:25:50.fits'

#     # fits_file_path = '../example_data/9h55m26.2_-58:25:47.fits'
#     # fits_file_path = '../example_data/9h02m23.3_-70:56:51.fits'
#     # fits_file_path = '../example_data/10h31m58.0_-58:14:42_BIS.fits'

#     fits_file_path = '../tests/example_data/11h12m10.1_-61:45:18.fits'
#     # fits_file_path = '../example_data/9h40m41.9_-66:39:16.fits'
#     # fits_file_path = '../example_data/2h54m35.0_-49:28:04.fits'
#     # fits_file_path = '../example_data/5h48m17.8_-15:11:58.fits'
#     # fits_file_path = '../example_data/6h45m53.0_10:03:48.fits'
#     # fits_file_path = '../example_data/8h32m33.2_-23:23:05.fits'

#     def get_coords(fitspath):
#         from astropy.coordinates import SkyCoord
#         from astropy import units as u
#         from pathlib import Path
#         try:
#             name = Path(fitspath).stem
#             ra, dec = name.split('_')[:2]
#             c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
#             return c.ra.degree, c.dec.degree
#         except:
#             return None
#     cc = get_coords(fits_file_path)

#     if cc is not None:

#         extra_args = {'ra_approx':cc[0], 'dec_approx': cc[1]}
#         print(extra_args)


#     sources, imageskysub = extract_stars(fits_file_path)
#     sources = sources[sources['flux']>2]
#     image = fits.getdata(fits_file_path)
#     from astropy.stats import sigma_clipped_stats
#     _, mm, _ = sigma_clipped_stats(image)
#     image = image - mm
#     # plot_stars(sources, image)        #%%
#     #import matplotlib.pyplot as plt
#     #plt.show(block=True)
#     # #%%
#     wcs = plate_solve_with_API(fits_file_path,  sources, **extra_args)
#     print((wcs['CD1_1']**2 + wcs['CD1_2']**2)**0.5 * 3600)
#     plot_stars_with_wcs(fits_file_path, sources)
