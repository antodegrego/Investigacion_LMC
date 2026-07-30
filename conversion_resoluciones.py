from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
import numpy as np

hdu12 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_12CO21.fits")
hdu13 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_13CO21.fits")
hdu0 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_ratio_13CO21_12CO21_13mask_only.fits")
hdu12_8 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment8_12CO21.fits")
hdu13_8 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment8_13CO21.fits")
hdu8 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_ratio_M8_13CO21_12CO21.fits")

data12 = hdu12[0].data
data13 = hdu13[0].data

# Header, cambiar 12 o 13 dependiendo
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

#12CO_15p3arcsec_Kkm_s.fits

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\12CO_15p3arcsec_Kkm_s.fits",
    data_conv,
    header,
    overwrite=True
)

# ==========================================================
# Comparar flujo antes y después de convolucionar
# ==========================================================

# Elige una región suficientemente grande
# alrededor de la fuente
y0, y1 = 100, 500
x0, x1 = 100, 700

data_original_region = data12[y0:y1, x0:x1]
data_conv_region = data_conv[y0:y1, x0:x1]

# Flujo integrado
flux_original = np.nansum(data_original_region)
flux_conv = np.nansum(data_conv_region)

print("Flujo integrado pre-convolución =", flux_original)
print("Flujo integrado post-convolución =", flux_conv)

print(
    "Fracción conservada =",
    flux_conv / flux_original
)