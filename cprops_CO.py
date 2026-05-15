import pycprops
import astropy.units as u
import os
from spectral_cube import SpectralCube


cubefile = 'NMolR_aca_12CO21_K_cube.fits' #Your cube
mask = 'NMolR_aca_image_v3_12CO21_mask_rms_fix.fits' # Mask defining where to find emission
d = 50.0 * u.kpc           # Distance (with units)


pycprops.fits2props(cubefile,
                    mask_file=mask,
                    distance=d, 
                    noise_file='noise_cube12.fits',
                    asgnname='12CO21.asgn.fits',
                    propsname='12CO21.props.fits')