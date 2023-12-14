from scipy.ndimage import median_filter
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
import sep

    
def extract_stars(fits_file_path):
    """
    Extract star positions from an image using SEP (Source Extractor as a Python library).
    
    Parameters:
    fits_file_path (str): Path to the FITS file.

    Returns:
    astropy.table.Table: Table of detected sources.
    numpy 2D array: background subtracted image
    """
    image = fits.getdata(fits_file_path).astype(float)
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
    
    
if __name__ == '__main__':


    
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
    
