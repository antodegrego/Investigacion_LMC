from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
import numpy as np

hdu12 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_moment0_12CO21_Wm2sr.fits")
hdu13 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_moment0_13CO21_Wm2sr.fits")


data12 = hdu12[0].data
data13 = hdu13[0].data

# Header, cambia 12 o 13 dependiendo
header = hdu13[0].header.copy()

fwhm_initial = 7.0      # arcsec
fwhm_target  = 15.3     # arcsec

fwhm_kernel = np.sqrt(
    fwhm_target**2 - fwhm_initial**2
)

pixel_scale = 1.2       # arcsec/pixel

sigma_pix = (
    fwhm_kernel / 2.355
) / pixel_scale

kernel = Gaussian2DKernel(sigma_pix)

data_conv = convolve(
    data13,
    kernel,
    boundary='extend',
    preserve_nan=True
)

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\13CO_15p3arcsec.fits",
    data_conv,
    header,
    overwrite=True
)