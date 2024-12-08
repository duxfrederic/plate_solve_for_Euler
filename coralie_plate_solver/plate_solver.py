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
import subprocess
from astroquery.astrometry_net import AstrometryNet
import logging
from pathlib import Path
import tempfile
from astropy.table import Table

from coralie_plate_solver.exceptions import CouldNotSolveError
from coralie_plate_solver.config import Config

config = Config()
logger = logging.getLogger(__name__)

REDO = False


def plate_solve_locally(fits_file_path, sources,
                        ra_approx=None, dec_approx=None,
                        scale_min=None, scale_max=None, use_n_brightest_only=None):
    """
    Calculate the WCS using a local installation of astrometry.net.
    The solve-field binary must be in the path, preferably in /usr/bin.
    You need to have it configured to use the appropriate index files.
    For wide field, a mixture of the 4100 series and 5200 series works well.

    Parameters:
    fits_file_path (Path or str): Path to the FITS file.
    sources (astropy.table.Table): Table of detected sources.
    use_existing_wcs_as_guess (bool): if a wcs information is already present in the fits file, use it to complete
                                      the rest of the arguments of the function. (ra_approx, dec_approx, etc.)
    ra_approx (float): Approximate RA in degrees.
    dec_approx (float): Approximate DEC in degrees.
    scale_min (float): lowest pixel scale to consider in arcsec/pixel
    scale_max (float): largest pixel scale to consider in arcsec/pixel
    use_n_brightest_only (int): number of sources to consider. If None using all.
    redo_if_done (bool): redo even if our solved keyword is already in the header?

    Returns:
    WCS header if successful, None otherwise.
    Also, updates the given fits file with the same WCS if successful.
    """
    fits_file_path = Path(fits_file_path)
    logger.info(f"plate_solve_locally on {fits_file_path}")
    solve_field_path = config.get('solve_field_path')
    n_brightest_stars_to_use = config.get('n_brightest_stars_to_use')

    if use_n_brightest_only is None:
        use_n_brightest_only = len(sources)

    # a work dir for the plate solving arguments
    with tempfile.TemporaryDirectory() as tmpdirname:
        # xylist.fits from sources
        xylist_path = Path(tmpdirname) / 'xylist.fits'
        t = Table([sources['xcentroid'], sources['ycentroid'], sources['flux']],
                  names=('X', 'Y', 'FLUX'))
        t.sort('FLUX', reverse=True)
        t = t[:use_n_brightest_only]
        hdu = fits.BinTableHDU(data=t)
        hdu.writeto(xylist_path, overwrite=True)

        # we'll need this command to use the system python interpreter, not the one running this script.
        # (unless solve-field was installed within this environment ...but likely not the case)
        new_env = os.environ.copy()
        #new_env["PATH"] = "/usr/bin:" + '/usr/local/astrometry/bin:' + new_env["PATH"]  # assuming python is in /usr/bin
        #new_env["PATH"] = '/home/euler/astrometrynet/bin:' + new_env["PATH"]  # assuming python is in /usr/bin

        # build solve-field command
        command = [solve_field_path, str(xylist_path), '--no-plots', '--x-column', 'X', '--y-column', 'Y',
                   '--sort-column', 'FLUX']  # by default solve-field sort by largest first, so ok to give flux this way
        # we also need the dimensions of the image:
        header = fits.getheader(fits_file_path)
        command += ['--width', str(header['NAXIS1']), '--height', str(header['NAXIS2'])]
        if ra_approx is not None and dec_approx is not None:
            command += ['--ra', str(ra_approx), '--dec', str(dec_approx), '--radius', '1']
        # add the scale estimate as well
        command += ['--scale-low', str(scale_min), '--scale-high', str(scale_max), '--scale-units', 'arcsecperpix']
        # limit the number of sources we consider
        command += ['--depth', str(n_brightest_stars_to_use)]

        # now call the local astrometry.net binary
        try:
            subprocess.run(command, check=True, cwd=tmpdirname, env=new_env)
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running solve-field: {e}")
            raise CouldNotSolveError('Error running solve-field.')

        # check for solution and update FITS file
        solved_path = next(Path(tmpdirname).glob('*.wcs'), None)
        if solved_path is None:
            logger.error("Astrometry failed: No solution found.")
            raise CouldNotSolveError('failed to solve astrometry')
        else:
            wcs = WCS(fits.getheader(solved_path)).to_header()
            # else, we're probably fine. Like, 99.99% confidence from my
            # experience with astrometry.net's plate solver.
            logger.info(f"plate_solve_with_API: {fits_file_path} solved, writing the WCS")
            with fits.open(fits_file_path, mode="update") as hdul:
                # little flag to indicate that we plate solved this file:
                wcs['PL-SLVED'] = 'done'
                # add all this info to the file:
                hdul[0].header.update(wcs)
                hdul.flush()
    return wcs


def plate_solve_with_API(fits_file_path, sources, ra_approx=None, dec_approx=None):
    """
    FALL BACK ONLY -- do not use in production if possible.
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
    n_brightest_stars_to_use = config.get('n_brightest_stars_to_use')
    # check, maybe we already solved it
    header = fits.getheader(fits_file_path)
    if 'PL-SLVED' in header and not REDO:
        logger.info(f"{fits_file_path} was already plate solved, no redoing.")
        # we already treated it.
        # the header contains the WCS info.
        return header

    # unpack the positions
    tupledetections = [(star['xcentroid'], star['ycentroid']) for star in sources[:n_brightest_stars_to_use]]
    x, y = list(zip(*tupledetections))

    # create a session
    R = requests.post('http://nova.astrometry.net/api/login',
                      #data={'request-json': json.dumps({"apikey": os.environ['astrometry_net_api_key']})})
                      data={'request-json': json.dumps({"apikey": "fyikoajaybczbugv"})})

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

