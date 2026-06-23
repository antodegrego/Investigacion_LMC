import pycprops
import astropy.units as u
import os
from spectral_cube import SpectralCube

# Permitir operaciones grandes
SpectralCube.allow_huge_operations = True

cubefile = 'NMolR_aca_12CO21_K_cube.fits' #Your cube
mask = 'NMolR_aca_13CO21_mask_3sigma.fits' # Mask defining where to find emission
d = 50.0 * u.kpc           # Distance (with units)


pycprops.fits2props(cubefile,
                    mask_file=mask,
                    distance=d, 
                    #noise_file='noise_cube12.fits',
                    asgnname='12CO21_mask13CO.asgn.fits',
                    propsname='12CO21_mask13CO.props.fits')