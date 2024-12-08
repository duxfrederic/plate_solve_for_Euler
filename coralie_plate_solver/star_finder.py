from scipy.ndimage import median_filter
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
import sep
import pathlib

    
def extract_stars(fits_file_path_or_2darray):
    """
    Extract star positions from an image using SEP (Source Extractor as a Python library).
    
    Parameters:
    fits_file_path (str): Path to the FITS file.

    Returns:
    astropy.table.Table: Table of detected sources.
    numpy 2D array: background subtracted image
    """
    if type(fits_file_path_or_2darray) is str or type(fits_file_path_or_2darray) is pathlib.PosixPath:
        image = fits.getdata(fits_file_path_or_2darray).astype(float)
    else:
        image = fits_file_path_or_2darray
    
    bkg = sep.Background(image, bw=64, bh=64, fw=3, fh=3)
    image_filtered = median_filter(image, size=2)

    image_sub = image_filtered - bkg
    objects = sep.extract(image_sub, thresh=4,  err=bkg.globalrms,
                          minarea=10)

    sources = Table()
    for col in objects.dtype.names:
        sources[col] = objects[col]
        
    # just to stick to the daostarfinder way
    sources['xcentroid'] = sources['x']
    sources['ycentroid'] = sources['y']

    # Sorting the sources by flux to have the brightest first
    sources.sort('flux', reverse=True)

    return sources, image_sub


def plot_stars(sources, image, savepath=None):
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
             color='red', label='Detected Sources', alpha=0.5,
             ms=10,
             mfc='None')
    plt.legend()
    plt.tight_layout()
    if savepath is None:
        plt.show()
    else:
        plt.savefig(savepath)
