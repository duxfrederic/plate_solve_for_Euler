#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 04:23:16 2023

@author: fred
"""
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
from pathlib import Path


def find_stars(fitspath):

    try:


        # Resolve star
        data = fits.getdata(fitspath)
        #base64_encoded = self.convert_array_to_png(data)


        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        #threshold = 1000000*(random.randint(0, 1)) + median + (5.0 * std)
        threshold = 10.0 * std
        daofind = DAOStarFinder(threshold=threshold, fwhm=5.0)
        sources = daofind(data - median)
        if sources is None:
            sorted_stars = sorted_stars_filtered = []
        else:

            sorted_stars = []
            for source in sources:
                x = source['xcentroid']
                y = source['ycentroid']
                flux = source['flux']
                sorted_stars.append({
                    "flux": int(flux),
                    "x": int(x),
                    "y": int(y)})

            sorted_stars.sort(key=lambda x: -x["flux"])
            # threshold: pretty high since exoplanets can usually only be studied around bright stars
            # so if another star in the field has less than a third of the flux of the brightest,
            # it is probably safe to assume it is not the one we are looking for.
            # so we discard anything that has less than the flux of the brightest star.
            # fred note: the stars are still saturated, we cannot really compare their flux.
            threshold = 0.5 * sorted_stars[0]["flux"]
            sorted_stars_filtered = [star for star in sorted_stars if star["flux"] > threshold]

        return sorted_stars, sorted_stars_filtered

    except Exception as E:
        print(E)

        return [], []
    
    

def process_fits_files(directory):
    fits_files = Path(directory).glob('*.fits')

    for fits_file in fits_files:
        data = fits.getdata(fits_file)

        # Use Astropy to apply a stretch for better visualization
        norm = simple_norm(data, 'linear', percent=99.5)
        plt.imshow(data, cmap='gray', norm=norm)

        # Get sources
        original_sources, filtered_sources = find_stars(fits_file)

        # Plot original sources
        for source in original_sources:
            plt.plot(source['x'], source['y'], 'o', mfc='None', markersize=10, markeredgecolor='red')

        # Plot filtered sources
        for source in filtered_sources:
            plt.plot(source['x'], source['y'], 'o', mfc='None', markersize=5, markeredgecolor='green')

        # Save the plot
        plt.colorbar()
        plt.title(fits_file.name)
        output_filename = fits_file.with_suffix('.png')
        plt.savefig(output_filename)
        plt.close()

process_fits_files('example_data/')



