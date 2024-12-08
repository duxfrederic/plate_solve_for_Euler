# A plate solver for Coralie's `bigeye` acquisition camera
Used during acquisition of Coralie targets: determines the position of the target on the
acquisition image; which can then be converted into a telescope offset to center the
target under the fiber. 

The logic is the following:
- acquisition image contains a lot of stars: plate solve acquisition image to find target
- acquisition image contains a few (< 6) stars with similar fluxes: cannot automatically handle this case for now, fall back to manual
- acquisition image contains a few (< 6) stars with one much brighter than the others: given that expolanet stars targets are typically bright, go for it. **The PI should specify a "manual" acquisition should their target star be next to a much brighter one.**
- acquisition image contains only one bright star: go for that one.


## Setup
This code is meant to be run on the `glsguiding` machine.

### Local installation of `astrometry.net`
I compiled `astrometry.net` at `/home/euler/compile_astrometrynet/astrometry.net-0.94`.
The `solve-field` executable, used by this package to determine astrometric solutions, is found at
`/home/euler/astrometrynet/bin`.
The index files (need 4000 and 5000 series for this to work) 
were downloaded at `/home/euler/astrometry_net_index_files`.

### YAML config file
The package expects to find an environment variable, `CORALIE_PLATE_SOLVER_CONFIG_PATH`, to point to the YAML 
file containing the configuration. Copy `sample_config.yaml` to the production directory, adapt it, and
set the environment variable above.

### Installation of the package
Your environment should contain the following packages:

- numpy
- scipy
- astropy -- for handling coordinates, cutouts ...
- photutils -- for aperture photometry, confirming presence of sources
- PyYAML -- for parsing the config file
- requests
- sep -- for source extraction
- astroquery -- for using the API service of astrometry.net at nova.astrometry.net. Not used by default.

I am not making the installation script of this package take care of its own dependencies to avoid conflicting versions.

Then, `pip install -e /path/to/location/of/this/repository`.

## Old: `scripts`
The `scripts.monitor_acoralie_for_changes.py` applies the logic above in parallel to observation,
monitoring the `acoralie.fits` file containing the acquisition image. 
Every time it is overwritten, a plot showing the position of the target is made at
`/home/remote/desktop/helper_astrometry_pointing`. 