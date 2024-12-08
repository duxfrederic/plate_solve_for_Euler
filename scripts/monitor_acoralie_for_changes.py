#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 01:37:19 2023

@author: fred

this script is a prototype running in parallel to observation.

it will add plate solved fields in a desktop directory of remote,
so the observer can use them for reference


after a few days of testing, we can incorporate it into coralie-autoguider-backend.
"""

import time
from pathlib import Path
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS
import logging

from coralie_plate_solver.telescope_reader import (get_telescope_position_skycoord_from_etcs,
                                                   get_current_catalogue_skycoord_from_tcs)
from coralie_plate_solver.star_finder import extract_stars
from coralie_plate_solver.process_acquisition_image import process_acquisition_image


WORKDIR = "/home/remote/Desktop/helper_astrometry_pointing"
logger = logging.getLogger(__name__)
logging.basicConfig(filename=WORKDIR + '/' + 'my_log.log', level=logging.INFO)


def diagnostic_plot(fits_file_path, sources, object_position, catalogue_skycoord):
    """
    Makes a diagnostic plot of our target localization.

    Args:
    - fits_file_path (str): Path to the acquisition FITS file.
    - sources (Table): Astropy Table of extracted sources.
    - object_position (tuple): Estimated position of the object.
    """

    # Load the FITS file
    with fits.open(fits_file_path) as hdu_list:
        image_data = hdu_list[0].data

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    z_scale = ZScaleInterval()
    v_min, v_max = z_scale.get_limits(image_data)
    
    # If plate solving worked, show the coordinates
    if 'PL-SLVED' in hdu_list[0].header:
        # make a new ax with the WCS projection
        logger.info('monitor_acoralie_for_changes.diagnostic_plot: PL-SLVED in header -> using WCS')
        wcs = WCS(hdu_list[0].header)
        ax = plt.subplot(projection=wcs)
        ax.coords.grid(True, color='white', ls='solid')
        ax.coords[0].set_axislabel('Right Ascension')
        ax.coords[1].set_axislabel('Declination')
        ax.plot([catalogue_skycoord.ra.deg], [catalogue_skycoord.dec.deg], 
                'o', mfc='None', label='Catalogue coordinates', ms=15,
                color='green', transform=ax.get_transform('world'))
    else:
        logger.info('script.monitor_acoralie_for_changes.diagnostic_plot: no wcs available for plot.')

    # image
    ax.imshow(image_data, origin='lower', cmap='gray', vmin=v_min, vmax=v_max)

    # Plot the extracted sources
    ax.scatter(sources['xcentroid'], sources['ycentroid'], s=30, edgecolor='red', 
               facecolor='none', label='Extracted Sources')

    # Plot the estimated position of the source
    if object_position is not None:
        ax.plot(object_position[0], object_position[1], 'x', color='blue', 
                markersize=10, label='Estimated Position')

    ax.legend()
    plt.tight_layout()

    # Save the plot
    plot_path = Path(fits_file_path).with_suffix('.jpg')
    plt.savefig(plot_path)


def get_date_obs(file_path):
    with fits.open(file_path) as hdu_list:
        date_obs = hdu_list[0].header['DATE-OBS']
    return date_obs.replace(':', '-').replace(' ', 'T')


def process_file_changes(source_file, destination_folder, last_modified):
    current_modified = source_file.stat().st_mtime
    if last_modified is None or current_modified > last_modified:
        logger.info("""\n\n CHANGE IN ACORALIE.FITS DETECTED, PROCESSING \n""")
        last_modified = current_modified
        # this means we just finished slewing, and a new
        # acquisition image was just written.
        # So, we need to
        # - get the telescope pointing for faster astrometry
        # - extract the sources in the image
        # - call our image processing routine
        # - write the plot to the workdir
        # 
        time.sleep(1) 
        # get the coords where telescope is pointing
        # not that useful anymore now that we have the catalogue 
        # coordinates. 
        # will write them down, so we can have diagnostic data for the pointing model.
        coord = get_telescope_position_skycoord_from_etcs()
        # get the coords of the current target
        target = get_current_catalogue_skycoord_from_tcs()

        date_obs = get_date_obs(source_file)
        destination_file = destination_folder / f"{date_obs}.fits"
        if not destination_file.exists():
            with fits.open(source_file, mode='readonly') as hdul_source:
                logger.info(f'writing file at {destination_file}')
                hdul_source.writeto(destination_file)
        
        # ok, we wrote our file ... we can now process it.
        # this process function expects the true coordinates of the object, so it can
        # provide its pixel position (to be used to calculate the offset)
        try:
            logger.info(f'running process_ac image on {destination_file} with {target.ra.deg}, {target.dec.deg}')
            object_position = process_acquisition_image(destination_file, target.ra.deg, target.dec.deg)
            logger.info(f'process_acquisition_image determined target position at {object_position}')
            
            # here we do redundant things (e.g. extracting sources 2 times, once
            # in extract_stars below and another in process_acquisition_image)
            # it's just for test
            sources, image_skysub = extract_stars(destination_file)
            logger.info(f'about to call diagnostic_plot with args: destination_file: {destination_file},'
                        f'{len(sources)} sources, object_position: {object_position}, target: {target}')
            diagnostic_plot(destination_file, 
                            sources, 
                            object_position,
                            target)
        except Exception as e:
            logger.info(f'Exception raised: {e}')

    return last_modified


def main():
    import socket
    hostname = socket.gethostname()
    if hostname == 'lenovolaptop':
        source_file = Path("/tmp/acoralie.fits")
        destination_folder = Path('/tmp/')

    else:
        source_file = Path("/raw/COR-GUIDING/acoralie.fits")
        destination_folder = Path(WORKDIR)
        
    last_modified = None

    while True:
        if source_file.exists():
            last_modified = process_file_changes(source_file, destination_folder, last_modified)
        time.sleep(2)


if __name__ == "__main__":
    main()

