from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from reproject import reproject_interp
from astropy.convolution import convolve, Gaussian2DKernel

#########################################
# CARGAR DATOS
#########################################
hdu_co = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_12CO21.fits")
hdu_cii = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\line_flux_CII.fits")

data_ciii = np.squeeze(hdu_cii[0].data)
data_co   = np.squeeze(hdu_co[0].data)

header_cii = hdu_cii[0].header

#########################################
# ESTIMAR RESOLUCIONES
#########################################

# tamaño de píxel en grados → arcsec
beam_co  = hdu_co[0].header["BMAJ"] * 3600  # en arcsec
beam_cii = abs(header_cii["CDELT1"]) * 3600  # aproximación

print("Beam CO:", beam_co)
print("Beam CII (approx):", beam_cii)

#########################################
# convertir FWHM → sigma
#########################################

sigma_co  = beam_co  / 2.355
sigma_cii = beam_cii / 2.355

#########################################
# kernel
#########################################

sigma_kernel = np.sqrt(max(0, sigma_cii**2 - sigma_co**2))

print("Sigma kernel:", sigma_kernel)

#########################################
# convertir a píxeles del mapa CO
#########################################

pix_co = abs(hdu_co[0].header["CDELT1"]) * 3600  # arcsec/pixel

sigma_pix = sigma_kernel / pix_co

kernel = Gaussian2DKernel(sigma_pix)

if sigma_kernel > 0:
    sigma_pix = sigma_kernel / pix_co
    kernel = Gaussian2DKernel(sigma_pix)
    data_co_smooth = convolve(data_co, kernel)
else:
    print("No se necesita convolución")
    data_co_smooth = data_co

#########################################
# SUAVIZAR CO
#########################################

data_co_smooth = convolve(data_co, kernel)

#########################################
# REPROYECTAR CO → GRILLA DE CIII
#########################################

co_reproj, footprint = reproject_interp(
    (data_co_smooth, hdu_co[0].header),
    header_cii
)

#########################################
# FUNCIÓN DE CORRELACIÓN LOCAL
#########################################

def local_correlation(map1, map2, size=7):
    pad = size // 2
    corr_map = np.full(map1.shape, np.nan)

    for i in range(pad, map1.shape[0] - pad):
        for j in range(pad, map1.shape[1] - pad):

            sub1 = map1[i-pad:i+pad+1, j-pad:j+pad+1]
            sub2 = map2[i-pad:i+pad+1, j-pad:j+pad+1]

            mask = (
                np.isfinite(sub1) &
                np.isfinite(sub2)
            )

            if np.sum(mask) > 5:
                x = sub1[mask]
                y = sub2[mask]

                if np.std(x) > 0 and np.std(y) > 0:
                    corr_map[i, j] = np.corrcoef(x, y)[0, 1]

    return corr_map

#########################################
# CALCULAR MAPA DE CORRELACIÓN
#########################################

corr_map = local_correlation(data_ciii, co_reproj, size=7)

#########################################
# GUARDAR RESULTADO
#########################################

fits.writeto(
    "correlation_map_CII_12CO.fits",
    corr_map,
    header_cii,
    overwrite=True
)

#########################################
# VISUALIZACIÓN
#########################################

plt.figure(figsize=(6,8))

im = plt.imshow(
    corr_map,
    origin="lower",
    cmap="coolwarm",
    vmin=-1,
    vmax=1
)

plt.colorbar(im, label="Pearson correlation")
plt.title("Spatial Correlation Map: CII vs 12CO")

plt.tight_layout()
plt.show()