# -*- coding: utf-8 -*-

import subprocess
from astropy.coordinates import SkyCoord
import astropy.units as u
from pathlib import Path

glscriptpath = Path('/opt/t4/beta/scripts/')


def get_skycoord_from_ETCS():
    position_script = str(glscriptpath / 'ETCS_axis_positions')
    try:
        result = subprocess.run([position_script], capture_output=True, text=True, check=True)
        output = result.stdout

        # parse the string output
        parts = output.split('#')
        alpha = float(parts[parts.index('alpha') + 1])
        delta = float(parts[parts.index('delta') + 1])

        sky_coord = SkyCoord(ra=alpha * u.degree, dec=delta * u.degree)

        return sky_coord

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
    except ValueError as e:
        print(f"Could not parse the output: {e}")

