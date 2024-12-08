# -*- coding: utf-8 -*-
import subprocess
from astropy.coordinates import SkyCoord
import astropy.units as u

from coralie_plate_solver.config import Config

config = Config()


def get_telescope_position_skycoord_from_etcs():
    command_name = 'ETCS_axis_positions'
    position_script = config.get('etcs_telescope_position_script')
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
        raise  # can't let it go through
    except ValueError as e:  # probably from parts.index, alpha and delta not in string
        print(f"ValueError: {e}, when parsing the output of {command_name}, it wasn't what we expected:", output)
        raise  # again, print above for context, but can't let it through
    except Exception as e:
        print(f"Undefined error when running {command_name}: {e}")
        raise  # same


def get_current_catalogue_skycoord_from_tcs():
    command_name = 'TCS_GetNode'
    position_script = config.get('tcs_script_node')
    arg_ra = config.get('tcs_ra_argument')
    arg_dec= config.get('tcs_dec_argument')
    try:
        ra = subprocess.run([position_script, arg_ra], capture_output=True, text=True, check=True)
        ra = float(ra.stdout)
        dec = subprocess.run([position_script, arg_dec], capture_output=True, text=True, check=True)
        dec = float(dec.stdout)

        sky_coord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)

        return sky_coord

    except subprocess.CalledProcessError as e:
        print(f"An error occurred in the subprocess when executing {command_name}: {e}")
        raise  # can't let it go through
    except Exception as e:
        print(f"Undefined error when running {command_name}: {e}")
        raise  # same
