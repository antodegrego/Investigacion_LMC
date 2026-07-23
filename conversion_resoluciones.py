from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
import numpy as np

hdu12 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_12CO21.fits")
hdu13 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_13CO21.fits")


data12 = hdu12[0].data
data13 = hdu13[0].data

# Header, cambia 12 o 13 dependiendo
header = hdu12[0].header.copy()

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
    data12,
    kernel,
    boundary='extend',
    preserve_nan=True
)

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\12CO_15p3arcsec_Kkm_s.fits",
    data_conv,
    header,
    overwrite=True
)