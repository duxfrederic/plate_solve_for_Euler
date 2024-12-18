#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 05:21:01 2023

@author: fred
"""
from pathlib import Path
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
import logging

from coralie_plate_solver.star_finder import extract_stars
from coralie_plate_solver.plate_solver import plate_solve_locally
from coralie_plate_solver.exceptions import (CouldNotSolveError, PlateSolvedButNoFluxAtObject,
                                             TooFewStarsForPlateSolving, EmptyFieldError,
                                             SuspiciousUniqueDetectionError,
                                             NotWithinImageBoundError)


logger = logging.getLogger(__name__)


def verify_object_position(sources_with_wcs, ra_obj, dec_obj, tolerance=5):
    """
    Verifies if the object of interest is in the list of sources.

    Args:
    - sources_with_wcs (Table): Astropy table of sources with RA and DEC.
    - ra_obj (float): Right Ascension of the object of interest.
    - dec_obj (float): Declination of the object of interest.
    - tolerance (int): Tolerance in pixels for matching the position.

    Returns:
    - Tuple (x, y): Pixel coordinates of the object if found, else None.
    """
    logger.info(f'verify_object_position at {ra_obj}, {dec_obj} on {len(sources_with_wcs)} sources.')

    # Create a SkyCoord object for the object of interest
    obj_coord = SkyCoord(ra=ra_obj * u.degree, dec=dec_obj * u.degree, frame='icrs')

    for source in sources_with_wcs:
        source_coord = SkyCoord(ra=source['RA']*u.degree, dec=source['DEC']*u.degree, frame='icrs')
        sep = obj_coord.separation(source_coord).arcsec

        if sep < tolerance:
            object_position_pixels = (source['xcentroid'], source['ycentroid'])
            return object_position_pixels

    return None


def check_flux_at_position(fits_file_path, ra_obj, dec_obj,
                           aperture_radius=5, annulus_radii=(8, 10), significance=10,
                           do_diagnostic_plot=True):
    """
    Performs aperture photometry at a specified position with background subtraction
    and optionally plots a diagnostic showing the aperture and annulus.

    Args:
    - fits_file_path (str): Path to the FITS file.
    - ra_obj (float): Right Ascension of the object of interest.
    - dec_obj (float): Declination of the object of interest.
    - aperture_radius (int): Radius for the aperture in pixels.
    - annulus_radii (tuple): Radii (inner and outer) for the background annulus in pixels.
    - significance (float): Flux significance (units of sigma) for claiming detection.
    - do_diagnostic_plot (bool): make a diagnostic plot?

    Returns:
    - bool: Is there significant flux at the position or not?
    """
    logger.info(f'check_flux_at_position on {fits_file_path} at {ra_obj}, {dec_obj}')
    with fits.open(fits_file_path) as hdulist:
        image_data = hdulist[0].data
        wcs = WCS(hdulist[0].header)

    _, median, std = sigma_clipped_stats(image_data)
    image_data = image_data - median  # Subtract the median background

    # transform to pixel coordinates
    x_obj, y_obj = wcs.world_to_pixel_values(ra_obj, dec_obj)

    # setup aperture photometry
    aperture = CircularAperture((x_obj, y_obj), r=aperture_radius)
    annulus_aperture = CircularAnnulus((x_obj, y_obj), r_in=annulus_radii[0], r_out=annulus_radii[1])
    
    # do the photometry to get the total flux within the aperture and annulus
    phot_table = aperture_photometry(image_data, [aperture, annulus_aperture])
    
    # get the mean background per pixel, then subtract it times area of aperture
    background_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
    background_sum = background_mean * aperture.area
    flux = phot_table['aperture_sum_0'] - background_sum
    d_flux = std * aperture.area**0.5  # error in the flux, since sigma is per pixel.
    if do_diagnostic_plot:
        z_scale = ZScaleInterval()
        v_min, v_max = z_scale.get_limits(image_data)
        plt.figure()
        plt.imshow(image_data, origin='lower', cmap='gray', vmin=v_min, vmax=v_max)
        aperture.plot(color='blue', lw=1.5, label='Aperture')
        annulus_aperture.plot(color='red', lw=1.5, linestyle='dashed', label='Background Annulus')
        plt.legend()
        plt.xlabel('X Pixel')
        plt.ylabel('Y Pixel')
        plt.tight_layout()
        fp = Path(fits_file_path)
        plot_path = fp.with_stem(f'{fp.stem}_aperture_diagnostic').with_suffix('.jpg')
        logger.info(f'check_flux_at_position on {fits_file_path}, saving diagnostic plot at {plot_path}')
        plt.savefig(plot_path, dpi=300)
        plt.close()

    # Check for significant flux in the aperture
    snr = (flux / d_flux > significance)
    return snr


def find_brightest_source(sources):
    """
    Finds the brightest source in the list.

    Args:
    - sources (Table): Astropy table of sources.

    Returns:
    - Row: The row of the brightest source in the table.
    """
    logger.info(f'find_brightest_source on {len(sources)} sources.')

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
    logger.info(f'is_significantly_brighter on {len(sources)} sources with {brightest_source}')

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
    logger.info(f'confirm_with_peak_detection on {fits_file_path}')

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


def find_star_that_clearly_pops_out(sources, brightness_factor=5):
    """
    given a list of stars, selects the brightest ones if all the others
    are much dimmer.

    Args:
    - sources: astropy Table with at least columns xcentroid, ycentroid, flux

    Returns:
    - Tuple (x, y): Pixel coordinates of the peak if confirmed, else raises
    
    Raises:
        raises TooFewStarsForPlateSolving indicates that we really can't be sure
    """
    logger.info(f'find_star_that_clearly_pops_out on {len(sources)} sources.')

    brightest_source = find_brightest_source(sources)
    if is_significantly_brighter(brightest_source, sources, brightness_factor):
        logger.info('found a very bright star! using this as our target.')
        object_position_pixels = (brightest_source['xcentroid'], brightest_source['ycentroid'])
        return object_position_pixels
    else:
        logger.info('did not find an obviously brighter star, cannot find object, giving up.')
        # if no, then it's really ambiguous, can't say anything with
        # any kind of confidence
        raise TooFewStarsForPlateSolving


def process_acquisition_image(fits_file_path, ra_obj, dec_obj):
    """
    Process the acquisition image to determine the position of the 
    object of interest in the field.

    Args:
    - fits_file_path (str): Path to the FITS file.
    - ra_obj (float): Right Ascension of the object of interest.
    - dec_obj (float): Declination of the object of interest.

    Returns:
    - Tuple (x, y): Pixel coordinates of the object in the image.
    
    Raises:
    - different exceptions if the object cannot be located or other conditions are not met.
    """
    logger.info(f'process_acquisition_image in {fits_file_path} at ({ra_obj:.05f}, {dec_obj:.05f})')

    # 1: Extract stars from the image
    sources, skysub_image = extract_stars(fits_file_path)
    sources = [] if sources is None else sources

    # 2: depends on the number of detections
    num_sources = len(sources)
    logger.info(f'process_acquisition_image on {fits_file_path} at  {ra_obj}, {dec_obj}:'
                f' extracted {num_sources} sources.')
    
    if num_sources > 6:
        # Try plate solving
        logger.info(f'trying plate solving on {fits_file_path}')
        try:
            wcs_header = plate_solve_locally(fits_file_path, sources,
                                             ra_approx=ra_obj,
                                             dec_approx=dec_obj,
                                             scale_min=0.15,
                                             scale_max=0.3)
            wcs = WCS(wcs_header)
            # we might need this as well
            image_header = fits.getheader(fits_file_path)
            # Convert pixel coordinates to RA and DEC and add them to the sources table
            sources['ra'] = None
            sources['dec'] = None
            for source in sources:
                ra, dec = wcs.all_pix2world(source['xcentroid'], source['ycentroid'], 1)
                source['ra'] = ra
                source['dec'] = dec
            # check that we have a detection near our object
            object_position_pixels = verify_object_position(sources, ra_obj, dec_obj)
            if object_position_pixels is not None:
                # if yes, then surely we can't be too wrong!
                logger.info('plate solved and found object at catalogue position.')
                return object_position_pixels
            else:
                # then ...our source extractor might have missed it, or it falls outside of the image.
                # check that it's inside the image:$
                logger.info('plate solved but no source corresponding to our object.')
                px_obj, py_obj = wcs.world_to_pixel_values(ra_obj, dec_obj)
                nx, ny = image_header.get('NAXIS1'), image_header.get('NAXIS2')
                is_within_bounds = (0 <= px_obj < nx) and (0 <= py_obj < ny)
                if not is_within_bounds:
                    logger.info('  object not within bounds of image.')
                    raise NotWithinImageBoundError
                logger.info('  object should be in the range of the image.')
                # ok, if it is in the footprint...
                # we check the aperture photometry at that location!
                there_is_flux = check_flux_at_position(fits_file_path, ra_obj, dec_obj, do_diagnostic_plot=False)
                if there_is_flux:
                    logger.info('  the source extractor probably missed it, there is flux there.')
                    # yay, we plate solved, and there is flux at the
                    # coordinates where our object should be.
                    # we're probably good.
                    object_position_pixels = wcs.world_to_pixel_values(ra_obj, dec_obj)
                    return object_position_pixels
                else:
                    logger.info(' no flux where the object should be. giving up.')
                    raise PlateSolvedButNoFluxAtObject
        except CouldNotSolveError:
            logger.info('could not plate solve, trying to see if there is an obviously brighter star.')
            # Well, here we can still do the old school select the obviously brighter star.
            return find_star_that_clearly_pops_out(sources)
            # this function will return the position of the obviously brighter star,
            # or raise another exception if no such star.

    elif num_sources > 1:
        logger.info(f'found {num_sources}')
        logger.info('no plate solving, checking sources')
        # well not great, only 2-6 sources ...
        # hard to plate solve, but we can probably assume the pointing
        # was good enough such that our target is in the field.
        # There is also a strong probability that our target
        # will be the brightest in the field.
        # So, check if one source is significantly brighter
        return find_star_that_clearly_pops_out(sources)

    elif num_sources == 1:
        logger.info('found a single star! checking that it is not a fluke with a different peak detection.')
        # that's probably our target if we're not too unlucky with the pointing.
        # let's confirm that it isn't a fluke by looking for the 
        # source with another technique.
        if confirm_with_peak_detection(fits_file_path, sources):
            logger.info('single star found with other technique as well, using this as target.')
            unique_source = sources[0]
            object_position_pixels = (unique_source['xcentroid'], unique_source['ycentroid'])
            return object_position_pixels
        else:
            logger.info('other technique found different peak, suspicious, giving up.')
            raise SuspiciousUniqueDetectionError

    else:
        logger.info('no star found in this field. giving up.')
        raise EmptyFieldError


def diagnostic_plot(fits_file_path, sources, object_position_pixels, ra_obj, dec_obj):
    """
    Makes a diagnostic plot of our target localization.

    Args:
    - fits_file_path (str): Path to the acquisition FITS file.
    - sources (Table): Astropy Table of extracted sources.
    - object_position_pixels (tuple): Estimated position of the object.
    - ra_obj (float): right ascension of the target object
    - dec_obj (float): declination of the target object.
    """
    # Load the image
    with fits.open(fits_file_path) as hdu_list:
        image_data = hdu_list[0].data

    # Set up the plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    z_scale = ZScaleInterval()
    v_min, v_max = z_scale.get_limits(image_data)

    # If plate solving worked, show the coordinates
    try:
        wcs = WCS(hdu_list[0].header)
        ax = plt.subplot(projection=wcs)
        ax.coords.grid(True, color='white', ls='solid')
        ax.coords[0].set_axislabel('Right Ascension')
        ax.coords[1].set_axislabel('Declination')
        
        ax.plot([ra_obj], [dec_obj], 'o', mfc='None', label='Catalogue coordinates', ms=15,
                color='green', transform=ax.get_transform('world'))
    except Exception as e:
        logger.info(f"Not showing position of the target on the plot: {e}")
    
    coo = SkyCoord(ra_obj, dec_obj, unit=u.degree)
    ra_str = coo.ra.to_string(unit=u.hourangle, sep=":", precision=2, pad=True)
    dec_str = coo.dec.to_string(sep=":", precision=2, alwayssign=True, pad=True)

    ax.set_title(ra_str + ' ' + dec_str)
    # Plot the image
    ax.imshow(image_data, origin='lower', cmap='gray', vmin=v_min, vmax=v_max)

    # Plot the extracted sources
    ax.scatter(sources['xcentroid'], sources['ycentroid'], s=30, edgecolor='red', 
               facecolor='none', label='Extracted Sources')

    # Plot the estimated position of the source
    if object_position_pixels is not None:
        ax.plot(object_position_pixels[0], object_position_pixels[1], 'x', color='blue',
                markersize=10, label='Estimated Position')

    ax.legend()
    plt.tight_layout()

    plot_path = Path(fits_file_path).with_suffix('.jpg')
    plt.savefig(plot_path)
