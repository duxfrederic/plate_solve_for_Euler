#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 05:21:01 2023

@author: fred
"""

import numpy as np
from scipy.ndimage import median_filter
from astropy.io import fits
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval


from euler_plate_solver.star_finder import extract_stars
from euler_plate_solver.plate_solver import plate_solve_with_API
from euler_plate_solver.exceptions import CouldNotSolveError



def verify_object_position(sources_with_wcs, RA_obj, DEC_obj, tolerance=5):
    """
    Verifies if the object of interest is in the list of sources.

    Args:
    - sources_with_wcs (Table): Astropy table of sources with RA and DEC.
    - RA_obj (float): Right Ascension of the object of interest.
    - DEC_obj (float): Declination of the object of interest.
    - tolerance (int): Tolerance in pixels for matching the position.

    Returns:
    - Tuple (x, y): Pixel coordinates of the object if found, else None.
    """

    # Create a SkyCoord object for the object of interest
    obj_coord = SkyCoord(ra=RA_obj*u.degree, dec=DEC_obj*u.degree, frame='icrs')

    for source in sources_with_wcs:
        source_coord = SkyCoord(ra=source['RA']*u.degree, dec=source['DEC']*u.degree, frame='icrs')
        sep = obj_coord.separation(source_coord).arcsec

        if sep < tolerance:
            return (source['xcentroid'], source['ycentroid'])

    return None



def check_flux_at_position(fits_file_path, RA_obj, DEC_obj,
                           aperture_radius=5, annulus_radii=(8, 10), significance=10):
    """
    Performs aperture photometry at a specified position with background subtraction
    and plots a diagnostic showing the aperture and annulus.

    Args:
    - fits_file_path (str): Path to the FITS file.
    - RA_obj (float): Right Ascension of the object of interest.
    - DEC_obj (float): Declination of the object of interest.
    - aperture_radius (int): Radius for the aperture in pixels.
    - annulus_radii (tuple): Radii (inner and outer) for the background annulus in pixels.
    - significance (float): Flux significance (units of sigma) for claiming detection.

    Returns:
    - bool: Is there significant flux at the position or not?
    """

    # Load the FITS file and WCS
    with fits.open(fits_file_path) as hdulist:
        image_data = hdulist[0].data
        wcs = WCS(hdulist[0].header)

    # Calculate sigma-clipped statistics for background estimation
    _, median, std = sigma_clipped_stats(image_data)
    image_data = image_data - median  # Subtract the median background

    # Transform RA and DEC to pixel coordinates
    x_obj, y_obj = wcs.all_world2pix([[RA_obj, DEC_obj]], 1)[0]

    # Define the aperture and annulus for photometry
    aperture = CircularAperture((x_obj, y_obj), r=aperture_radius)
    annulus_aperture = CircularAnnulus((x_obj, y_obj), r_in=annulus_radii[0], r_out=annulus_radii[1])
    
    # Perform photometry to get the total flux within the aperture and annulus
    phot_table = aperture_photometry(image_data, [aperture, annulus_aperture])
    
    # Calculate the mean background per pixel within the annulus
    background_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
    
    # Calculate the total background within the aperture
    background_sum = background_mean * aperture.area
    
    # Subtract the background from the aperture sum
    flux = phot_table['aperture_sum_0'] - background_sum
    dflux = std * aperture.area**0.5  # Error in the flux

    # Plotting section: display the image, aperture, and annulus
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(image_data)
    plt.imshow(image_data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
    aperture.plot(color='blue', lw=1.5, label='Aperture')
    annulus_aperture.plot(color='red', lw=1.5, linestyle='dashed', label='Background Annulus')
    plt.legend()
    plt.xlabel('X Pixel')
    plt.ylabel('Y Pixel')
    
    # Save the plot
    plot_path = Path(fits_file_path).with_suffix('.png')
    plt.savefig(plot_path)
    plt.close()

    # Check for significant flux in the aperture
    return (flux / dflux > significance)


def find_brightest_source(sources):
    """
    Finds the brightest source in the list.

    Args:
    - sources (Table): Astropy table of sources.

    Returns:
    - Row: The row of the brightest source in the table.
    """

    # Assuming 'flux' column exists and represents the brightness
    brightest_index = np.argmax(sources['flux'])
    return sources[brightest_index]



def is_significantly_brighter(brightest_source, sources, brightness_factor=5):
    """
    Is the brightest source is significantly brighter than others?

    Args:
    - brightest_source (Row): The row of the brightest source in the table.
    - sources (Table): Astropy table of sources.
    - brightness_factor (int): Factor by which the brightest source should be brighter.

    Returns:
    - bool: True if the brightest source is significantly brighter, False otherwise.
    """

    brightest_flux = brightest_source['flux']

    for source in sources:
        if source == brightest_source:
            continue
        if brightest_flux < brightness_factor * source['flux']:
            return False

    return True




def confirm_with_peak_detection(fits_file_path, source, window_size=5, peak_tolerance=5):
    """
    Confirms the detection of a single source using peak detection.

    Args:
    - fits_file_path (str): Path to the FITS file.
    - source (Row): The row of the single source detected.
    - window_size (int): Size of the median filter window.
    - peak_tolerance (int): Tolerance in pixels for confirming the peak.

    Returns:
    - Tuple (x, y): Pixel coordinates of the peak if confirmed, else None.
    """

    # Load the FITS image data
    with fits.open(fits_file_path) as hdulist:
        image_data = hdulist[0].data

    # median filter to remove hot pixels or other unholy manifestations of bad luck
    filtered_data = median_filter(image_data, size=window_size)

    # then find the peak
    y_peak, x_peak = np.unravel_index(np.argmax(filtered_data), filtered_data.shape)

    # check that we're close.
    distance2 = (x_peak - source['xcentroid'])**2 + (y_peak - source['ycentroid'])**2
    return distance2 < peak_tolerance**2




def process_acquisition_image(fits_file_path, RA_obj, DEC_obj):
    """
    Process the acquisition image to determine the position of the 
    object of interest in the field.

    Args:
    - fits_file_path (str): Path to the FITS file.
    - RA_obj (float): Right Ascension of the object of interest.
    - DEC_obj (float): Declination of the object of interest.

    Returns:
    - Tuple (x, y): Pixel coordinates of the object in the image.
    - Raises CouldNotSolveError if the object cannot be located or other conditions are not met.
    """

    try:
        # 1: Extract stars from the image
        sources = extract_stars(fits_file_path)
        sources = [] if sources is None else sources
        
        
        # 2: depends on the number of detections
        num_sources = len(sources)

        if num_sources > 6:
            # Try plate solving
            try:
                wcs_header = plate_solve_with_API(fits_file_path, sources)
                wcs = WCS(wcs_header)
        
                # Convert pixel coordinates to RA and DEC and add them to the sources table
                sources['RA'] = None
                sources['DEC'] = None
                for source in sources:
                    ra, dec = wcs.all_pix2world(source['xcentroid'], source['ycentroid'], 1)
                    source['RA'] = ra
                    source['DEC'] = dec                
                # check that we have a detection near our object
                object_position = verify_object_position(sources, RA_obj, DEC_obj)
                if object_position:
                    # if yes, then surely we can't be too wrong!
                    return object_position
                else:
                    # then ...our source extractor might have missed it.
                    # we check the aperture photometry at that location!
                    thereisflux = check_flux_at_position(fits_file_path, RA_obj, DEC_obj)
                    if thereisflux:
                        # yay, we plate solved, and there is flux at the
                        # coordinates where our object should be.
                        # we're probably good.
                        object_position = wcs.world_to_pixel_values(RA_obj, DEC_obj)
                        return object_position
                    else:
                        mm = "Plate solve succeeded, but no flux where our object should be!"
                        raise CouldNotSolveError(mm)
            except CouldNotSolveError:
                raise

        elif num_sources > 1:
            # well not great, only 2-6 sources ...
            # hard to plate solve, but we can probably assume the pointing
            # was good enough such that our target is in the field.
            # There is also a strong probability that our target
            # will be the brightest in the field.
            # So, check if one source is significantly brighter
            brightest_source = find_brightest_source(sources)
            if is_significantly_brighter(brightest_source, sources):
                return (brightest_source['xcentroid'], brightest_source['ycentroid'])
            else:
                # if no, then it's really ambigous, can't say anything with
                # any kind of confidence
                raise CouldNotSolveError("No dominant bright source found")

        elif num_sources == 1:
            # that's probably our target if we're not too unlucky with the pointing.
            # let's confirm that it isn't a fluke by looking for the 
            # source with another technique.
            if confirm_with_peak_detection(fits_file_path, sources):
                uniquesource = sources[0]
                return (uniquesource['xcentroid'], uniquesource['ycentroid'])
            else:
                raise CouldNotSolveError("Unique source in the field, but it's too sus.")

        else:
            raise CouldNotSolveError("Seemingly empty field")

    except Exception as e:
        # place holder for future issues ...
        raise e



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.visualization import ZScaleInterval
    from astropy.wcs import WCS
    from pathlib import Path
    
    def get_coords(fitspath):
        from astropy.coordinates import SkyCoord
        from astropy import units as u
        from pathlib import Path
        try:
            name = Path(fitspath).stem
            ra, dec = name.split('_')[:2]
            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            return c.ra.degree, c.dec.degree
        except:
            return None
        
    def diagnostic_plot(fits_file_path, sources, object_position, RA_obj, DEC_obj):
        """
        Makes a diagnostic plot of our target localization.
    
        Args:
        - fits_file_path (str): Path to the acquisition FITS file.
        - sources (Table): Astropy Table of extracted sources.
        - object_position (tuple): Estimated position of the object.
        - true_position (tuple, optional): True position of the object if WCS succeeded.
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
            
            ax.plot([RA_obj], [DEC_obj], 'o', mfc='None', label='Catalogue coordinates', ms=15,
                    color='green', transform=ax.get_transform('world'))
        except:
            pass
        
        coo = SkyCoord(RA_obj, DEC_obj, unit=u.degree)
        rastr = coo.ra.to_string(unit=u.hourangle, sep=":", precision=2, pad=True)
        decstr = coo.dec.to_string(sep=":", precision=2, alwayssign=True, pad=True)

        ax.set_title(rastr + ' ' + decstr)
        # Plot the image
        ax.imshow(image_data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
    
        # Plot the extracted sources
        ax.scatter(sources['xcentroid'], sources['ycentroid'], s=30, edgecolor='red', facecolor='none', label='Extracted Sources')
    
        # Plot the estimated position of the source
        ax.plot(object_position[0], object_position[1], 'x', color='blue', markersize=10, label='Estimated Position')
    
        # Plot the true position of the source, if available
        # if true_position:
            # ax.plot(true_position[0], true_position[1], marker='o', color='green', markersize=10, label='True Position')
    
        # Add a legend
        ax.legend()
        plt.tight_layout()
    
        # Save the plot
        plot_path = Path(fits_file_path).with_suffix('.jpg')
        plt.savefig(plot_path)
        # plt.close()
    
    # Example usage
    if __name__ == '__main__':
        dd = Path('../example_data')
        for fits_file_path in dd.glob('*.fits'):
            if not '5h48' in fits_file_path.name:
                continue
    
            # Assume these functions and variables are defined as per your previous logic
            sources = extract_stars(fits_file_path)
            RA_obj, DEC_obj = get_coords(fits_file_path)
            object_position = process_acquisition_image(fits_file_path, RA_obj, DEC_obj)
            # true_position = None  # You need to determine this based on WCS solving
        
            diagnostic_plot(fits_file_path, sources, object_position, RA_obj, DEC_obj)
