import subprocess
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt


class SExtractorRunner:
    def __init__(self, sextractor_path, config_path):
        self.sextractor_path = sextractor_path
        self.config_path = config_path

    def run_sextractor(self, imagepath, pixsize, saturlevel, catfilename):
        cmd = f"{self.sextractor_path} {imagepath} -c {self.config_path}/template.sex"
        cmd += f" -PIXEL_SCALE {pixsize:.3f} "
        cmd += f"-SATUR_LEVEL {saturlevel:.3f} "
        cmd += f"-CATALOG_NAME {catfilename}"
        subprocess.call(cmd, shell=True)

    def read_catalog(self, catfilename):
        with fits.open(catfilename, ignore_missing_simple=True) as hdul:
            data = hdul[1].data  # Assuming the data is in the first extension
            return data


def extract_stars(fits_file_path):
    """
    Extract star positions from an image using DAOStarFinder.
    We are quite strict here, requiring 40 sigmas detections.
    The fwhm is tuned on the average seeing at Euler, ~1.2-1.4 arcsec.

    Parameters:
    fits_file_path (str): Path to the FITS file.

    Returns:
    astropy.table.Table: Table of detected sources.
    """
    image = fits.getdata(fits_file_path).astype(float)
    mean, median, stddev = sigma_clipped_stats(image)
    image -= median

    daofind = DAOStarFinder(fwhm=7.0, threshold=40*stddev)
    sources = daofind(image)
    return sources


def plot_stars(sources, image):
    """
    Plot the image with detected sources marked (for debugging).

    Parameters:
    sources (astropy.table.Table): Table of detected sources.
    image (numpy.ndarray): Image data.
    """
    plt.figure(figsize=(15, 15))
    interval = ZScaleInterval(contrast=0.1)
    vmin, vmax = interval.get_limits(image)
    plt.imshow(image, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.plot(sources['xcentroid'], sources['ycentroid'], 'o',
             color='red', label='Detected Sources',
             mfc='None')
    plt.legend()
    plt.show()
    
    
if __name__ == '__main__':

    im = "../example_data/10h31m58s_-58d14m42s.fits"
    im = "../example_data/plate_solve_try3.fits"
    im = '../example_data/9h02m_-70d56m.fits'
    # im = '../example_data/ECAM.2023-12-11T00:32:18.000.fits'
    sextractor = SExtractorRunner("sex", ".")
    sextractor.run_sextractor(im, 0.25, 4096, "output.fits")
    catalog_data = sextractor.read_catalog("output.fits")
    
    
    
    def plot_stars_on_image(image_file, catalog_data):
        # Load the image data
        with fits.open(image_file) as hdul:
            image_data = hdul[0].data
    
        # Prepare the plot
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
        # Adjust the image contrast using ZScaleInterval
        interval = ZScaleInterval(contrast=0.1)
        vmin, vmax = interval.get_limits(image_data)
        ax.imshow(image_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
    
        # Plot the stars and annotate with flags
        for star in catalog_data:
            x, y, flag = star['X_IMAGE'], star['Y_IMAGE'], star['FLAGS']
            flux = star['FLUX_APER'][0]/100
            if flag > 3:
                continue
            ax.plot(x, y, 'ro', mfc='None', ms=int(flux**0.5))  # red circle for each star
    
        ax.set_xlabel('X Pixel')
        ax.set_ylabel('Y Pixel')
        plt.title('Stars with Flag Numbers')
        plt.show()
    
    plot_stars_on_image(im, catalog_data)
    
