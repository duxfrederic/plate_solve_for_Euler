import unittest
from pathlib import Path
from shutil import rmtree
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS

from euler_plate_solver.star_finder import extract_stars
from euler_plate_solver.plate_solver import plate_solve_locally


def parse_ra_dec_from_filename(filename):
    """
    Extracts RA and Dec from the given filename formatted as 'HHhMMmSS.S_DD:MM:SS.fits'

    Parameters:
    filename (str): Filename to parse.

    Returns:
    tuple: (RA in degrees, Dec in degrees)
    """
    ra_str, dec_str = filename.replace('.fits', '').split('_')
    ra_str = ra_str.replace('h', 'h ').replace('m', 'm ')
    dec_str = dec_str.replace(':', ' ')

    coord = SkyCoord(ra_str + ' ' + dec_str, unit=(u.hourangle, u.deg))
    return coord.ra.degree, coord.dec.degree


class TestPlateSolveLocally(unittest.TestCase):

    def test_plate_solve_locally(self):
        filename = '6h45m53.0_10:03:48.fits'
        fits_file_path = Path('example_data') / filename  # Update path as needed
        tempdir = Path('/tmp/tmp_astro')
        if tempdir.exists():
            rmtree(tempdir)
        tempdir.mkdir(exist_ok=True)
        # write the array without wcs elsewhere:
        data = fits.getdata(fits_file_path)
        newpath = tempdir / f"{fits_file_path.name}"
        fits.writeto(tempdir / newpath, data, overwrite=True)

        sources, image_sub = extract_stars(newpath)
        ra_approx, dec_approx = parse_ra_dec_from_filename(newpath.name)

        wcs_header = plate_solve_locally(newpath, sources, ra_approx, dec_approx)

        wcs = WCS(fits.getheader(newpath))

        self.assertIsNotNone(wcs_header, "WCS header should not be None.")
        self.assertIn('PL-SLVED', wcs_header, "WCS header should contain 'PL-SLVED' key.")


# Running the test
if __name__ == '__main__':
    unittest.main()
