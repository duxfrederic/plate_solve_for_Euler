#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 08:51:53 2023

@author: fred
"""
import numpy as np

def pixel_to_azi_ele(px_x, px_y):


    d_xscale = 0.2575#d_xpixel / te_scale
    # x(y)scales      echelle "/pixel
    d_yscale = 0.2575#d_ypixel / te_scale
    # x(y)scales      echelle "/pixel

    d_xscazin = -.17452406e-01 #-sindeg(1)
    # Matrice normalisee de transformation X-Y -> azimut-elevation
    d_yscazin = 0.99984771 #cosdeg(1)
    # Correspond a la rotation de ~1 degree des axes avec
    d_xscelen = -.99984771 #-cosdeg(1)
    # Azimut || Y+ et elevation || X-
    d_yscelen = -.17452406e-01 #-sindeg(1)
    # i.e. un mouvement dx/dy  d'image est obtenu par mouvement

    d_xscazi = d_xscazin * d_xscale
    d_yscazi = d_yscazin * d_yscale
    d_xscele = d_xscelen * d_xscale
    d_yscele = d_yscelen * d_yscale

    ce_argazi = (px_x * d_xscazi + px_y * d_yscazi)
    ce_argele = (px_x * d_xscele + px_y * d_yscele)

    return (ce_argazi,ce_argele)




def pixel_to_azi_ele_dep(px_x, px_y, altitude_of_center_deg):
    # Constants and scales
    d_xscale = 0.2575
    d_yscale = 0.2575

    # Rotation matrix components
    d_xscazin = -0.017452406
    d_yscazin = 0.99984771
    d_xscelen = -0.99984771
    d_yscelen = -0.017452406

    # Apply scales
    d_xscazi = d_xscazin 
    d_yscazi = d_yscazin 
    d_xscele = d_xscelen 
    d_yscele = d_yscelen 

    # Calculate rotated pixel offsets
    ce_argazi_rotated = (px_x * d_xscazi + px_y * d_yscazi)
    ce_argele_rotated = (px_x * d_xscele + px_y * d_yscele)

    from astropy.wcs import WCS

    def setup_wcs(pix_scale_deg, center):
        # Initialize a WCS object with 2 axes (Altitude and Azimuth)
        w = WCS(naxis=2)

        # Set the pixel scale in degrees/pixel
        w.wcs.cdelt = [pix_scale_deg, pix_scale_deg]

        # Define the coordinate system type: 'RA---TAN' and 'DEC--TAN' can be used
        # as proxies for Altitude and Azimuth
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

        # Set the reference pixel to the center of your image
        # For example, if your image is 1000x1000 pixels:
        w.wcs.crpix = [1.0, 1.0]

        # Set the reference coordinate to the center pointing of your telescope
        # For example, if pointing to zenith (Alt=90 degrees, Az=0 degrees):
        w.wcs.crval = center

        return w

    def pixel_to_alt_az(wcs, px_x, px_y):
        # Convert pixel coordinates to world coordinates
        alt_az = wcs.pixel_to_world_values(px_x, px_y)
        return alt_az
    
    
    # Setup WCS with pixel scale (convert arcsec to degrees)
    pix_scale_deg = 0.2575 / 3600
    wcs = setup_wcs(pix_scale_deg, [altitude_of_center_deg, 0.])


    # Get Altitude and Azimuth
    alt, az = pixel_to_alt_az(wcs, ce_argazi_rotated, ce_argele_rotated)
    delta_alt_arcsec = (alt-altitude_of_center_deg)*3600
    delta_az_arcsec = az * 3600
    # # Convert arcseconds to radians
    # ce_argazi_rad = np.radians(ce_argazi_rotated / 3600)
    # ce_argele_rad = np.radians(ce_argele_rotated / 3600)
    
    # altitude_of_center_rad = np.radians(altitude_of_center_deg)
    # # Deprojection to sphere
    # delta_az = np.arctan2(np.sin(ce_argazi_rad), np.cos(ce_argazi_rad) * np.sin(altitude_of_center_rad) + np.tan(ce_argele_rad) * np.cos(altitude_of_center_rad))
    # delta_alt = np.arcsin(np.cos(ce_argele_rad) * np.sin(altitude_of_center_rad) - np.sin(ce_argele_rad) * np.cos(altitude_of_center_rad) * np.cos(ce_argazi_rad))

    # # Convert radians to arcseconds or degrees as needed
    # delta_az_arcsec = np.degrees(delta_az) * 3600
    # delta_alt_arcsec = np.degrees(delta_alt) * 3600

    return delta_az_arcsec, delta_alt_arcsec

# Example usage
# altitude_of_center_rad = np.radians(45) # Example altitude, adjust as needed
# px_x, px_y = 10, 20 # Example pixel offsets
# azimuth, elevation = pixel_to_azi_ele(px_x, px_y, altitude_of_center_rad)


#%%




origalt, origaz = 30., 0.

px_x = 1
px_y = 0


daz, dalt = pixel_to_azi_ele_dep(px_x, px_y, origalt)
print(f"Altitude: {dalt:.02f} arcsec, Azimuth: {daz:.02f} arcsec")
daz, dalt = pixel_to_azi_ele(px_x, px_y)
print(f"Altitude: {dalt:.02f} arcsec, Azimuth: {daz:.02f} arcsec")
