# -*- coding: utf-8 -*-

import subprocess
from astropy.coordinates import SkyCoord
import astropy.units as u
from pathlib import Path

glscriptpath = Path('/opt/t4/beta/scripts/')


def get_skycoord_from_ETCS():
    command_name =  'ETCS_axis_positions'
    position_script = str(glscriptpath / command_name)
    try:
        result = subprocess.run([position_script], capture_output=True, text=True, check=True)
        output = result.stdout

        # parse the string output
        separator = output[0]
        parts = output.split(separator)
        alpha = float(parts[parts.index('alpha') + 1])
        delta = float(parts[parts.index('delta') + 1])

        sky_coord = SkyCoord(ra=alpha * u.degree, dec=delta * u.degree)

        return sky_coord

    except subprocess.CalledProcessError as e:
        print(f"An error occurred in the subprocess when executing {command_name}: {e}")
    except ValueError as e: # probably from parts.index, alpha and delta not in string
        print(f"ValueError: {e}, when parsing the output of {command_name}, it wasn't what we expected:", output)
    except Exception as e:
        print(f"Undefined error when running {command_name}: {e}")
        raise # the above print just for more info, but can't let it go through.

