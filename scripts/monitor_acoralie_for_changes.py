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

from euler_plate_solver.telescope_reader import get_skycoord_from_ETCS
from euler_plate_solver.star_finder import extract_stars
from euler_plate_solver.process_acquisition_image import process_acquisition_image, diagnostic_plot
from euler_plate_solver.exceptions import PlateSolvedButNoFluxAtObject


WORKDIR = "/home/remote/Desktop/helper_astrometry_pointing"


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
                coord = get_skycoord_from_ETCS()
                
                
                date_obs = get_date_obs(source_file)
                destination_file = destination_folder / f"{date_obs}.fits"
                destination_file.write_bytes(source_file.read_bytes())
                
                # ok, we wrote out file ... we can now process it.
                # here we do redundant things (e.g. extracting sources 2 times, once
                # in extract_stars below and another in process_acquisition_image)
                # it's just for test
                sources = extract_stars(destination_file)
                
                # WARNING THIS IS NOT HOW THE FUNCTION BELOW IS TO BE CALLED ULTIMATELY
                # this function expects the true coordinates of the object, so it can
                # provide its pixel position (to be used to calculate the offset)
                # here we provide the pointing coordinates
                # thus, the second exception will always be raised unless the telescope
                # pointed exactly on a star.
                object_position = None
                try:
                    object_position = process_acquisition_image(destination_file, coord.ra.deg, coord.dec.deg)
                except PlateSolvedButNoFluxAtObject:
                    # this is normal due to the warning above.
                    # we just keep going.
                    pass
                except Exception as e:
                    print(e)
                    # something else,just a print for debug for now, 
                    # I don't know where the best place would be for logging.
                # true_position = None  # You need to determine this based on WCS solving
                
                
                diagnostic_plot(destination_file, sources, object_position, coord.ra.deg, coord.dec.deg)

        time.sleep(5)


if __name__ == "__main__":
    main()