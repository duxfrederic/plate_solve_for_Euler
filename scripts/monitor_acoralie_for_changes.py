#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 01:37:19 2023

@author: fred

this script is a prototype running in parallel to observation.

it will add plate solved fields in a desktop directory of remote,
so the observer can use them for reference


after a a few days of testing, we can incorporate it into coralie-autoguider-backend.


WE WILL NEED A WAY TO GET THE CATALOGUE COORDINATES OF THE OBJECT
WE WANT TO SLEW TO.
"""

import time
from pathlib import Path
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS


from euler_plate_solver.telescope_reader import (get_telescope_position_skycoord_from_ETCS,
                                                 get_current_catalogue_skycoord_from_TCS)
from euler_plate_solver.star_finder import extract_stars
from euler_plate_solver.process_acquisition_image import process_acquisition_image
from euler_plate_solver.exceptions import PlateSolvedButNoFluxAtObject


WORKDIR = "/home/remote/Desktop/helper_astrometry_pointing"



def diagnostic_plot(fits_file_path, sources, object_position, catalogue_skycoord):
    """
    Makes a diagnostic plot of our target localization.

    Args:
    - fits_file_path (str): Path to the acquisition FITS file.
    - sources (Table): Astropy Table of extracted sources.
    - object_position (tuple): Estimated position of the object.
    """

    # Load the FITS file
    with fits.open(fits_file_path) as hdulist:
        image_data = hdulist[0].data

    # Set up the plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(image_data)

    # If plate solving worked, show the coordinates
    try:
        wcs = WCS(hdulist[0].header)
        ax = plt.subplot(projection=wcs)
        ax.coords.grid(True, color='white', ls='solid')
        ax.coords[0].set_axislabel('Right Ascension')
        ax.coords[1].set_axislabel('Declination')
        # ax.coords[0].set_ticks(number=20)
        # ax.coords[1].set_ticks(number=20)
        ax.plot([catalogue_skycoord.ra.deg], [catalogue_skycoord.dec.deg], 
                'o', mfc='None', label='Catalogue coordinates', ms=15,
                color='green', transform=ax.get_transform('world'))
        
    except:
        pass
    
    # Plot the image
    ax.imshow(image_data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)

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
    with fits.open(file_path) as hdul:
        date_obs = hdul[0].header['DATE-OBS']
    return date_obs.replace(':', '-').replace(' ', 'T')  

def main():
    source_file = Path("/raw/COR-GUIDING/acoralie.fits")
    destination_folder = Path(WORKDIR)
    last_modified = None

    while True:
        if source_file.exists():
            current_modified = source_file.stat().st_mtime
            if last_modified is None or current_modified > last_modified:
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
                # will write them down so we can have diagnostic data
                # for the pointing model.
                coord = get_telescope_position_skycoord_from_ETCS()
                # get the coords of the current target
                target = get_current_catalogue_skycoord_from_TCS()
                
                
                date_obs = get_date_obs(source_file)
                destination_file = destination_folder / f"{date_obs}.fits"
                if not destination_file.exists():
                    destination_file.write_bytes(source_file.read_bytes())
                
                # ok, we wrote out file ... we can now process it.
                # this process function expects the true coordinates of the object, so it can
                # provide its pixel position (to be used to calculate the offset)
                object_position = None
                try:
                    object_position = process_acquisition_image(destination_file, target.ra.deg, target.dec.deg)
                except Exception as e:
                    print(e)
                    raise
                    # something else,just a print for debug for now, 
                    # I don't know where the best place would be for logging.
                # true_position = None  # You need to determine this based on WCS solving
                
                # here we do redundant things (e.g. extracting sources 2 times, once
                # in extract_stars below and another in process_acquisition_image)
                # it's just for test
                sources, imageskysub = extract_stars(destination_file)
                diagnostic_plot(destination_file, sources, object_position,
                                target)

        time.sleep(5)


if __name__ == "__main__":
    main()

