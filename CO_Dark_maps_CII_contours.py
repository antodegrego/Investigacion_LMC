from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt


####################################
# Cargar mapas
####################################

cii_hdu = fits.open(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\line_flux_CII.fits"
)

co12_hdu = fits.open(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_12CO_15p3arcsec_Wm2sr_reprojected.fits"
)

co13_hdu = fits.open(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_13CO_15p3arcsec_Wm2sr_reprojected.fits"
)


cii = np.squeeze(cii_hdu[0].data)
co12 = np.squeeze(co12_hdu[0].data)
co13 = np.squeeze(co13_hdu[0].data)

header_cii = cii_hdu[0].header
header_12 = co12_hdu[0].header
header_13 = co13_hdu[0].header


####################################
# Revisar dimensiones
####################################

print("================================")
print("DIMENSIONES")
print("================================")

print("CII  :", cii.shape)
print("12CO :", co12.shape)
print("13CO :", co13.shape)


####################################
# RMS de CII
####################################

# IMPORTANTE:
# Esta región debe estar libre de emisión de CII.

noise_region_cii = cii[330:380, 22:77]

valid_cii = noise_region_cii[np.isfinite(noise_region_cii)]

if len(valid_cii) == 0:
    raise ValueError(
        "La región de ruido del CII contiene solamente NaN."
    )

rms_cii = np.std(valid_cii)


print()
print("================================")
print("RMS CII")
print("================================")

print(
    f"RMS CII = {rms_cii:.4e} "
    r"W/m$^2$/sr"
)


####################################
# Niveles de CII
####################################

# Múltiplos de RMS que quieres usar

levels_xrms_cii = np.array([
    3, 6, 9, 12
])

# Convertir a unidades del mapa

levels_cii = levels_xrms_cii * rms_cii


####################################
# Mostrar niveles
####################################

print()
print("================================")
print("NIVELES CII")
print("================================")

for n, level in zip(levels_xrms_cii, levels_cii):

    print(
        f"{n:3.0f} RMS = {level:.4e} W/m2/sr"
    )


####################################
# Comprobar niveles disponibles
####################################

print()
print("Máximo CII =", np.nanmax(cii))

levels_cii_valid = levels_cii[
    levels_cii <= np.nanmax(cii)
]

print("Niveles utilizados:")
print(levels_cii_valid)


####################################
# WCS
####################################

wcs_12 = WCS(header_12)
wcs_13 = WCS(header_13)


####################################
# MAPA 1
# 12CO + contornos CII
####################################

fig = plt.figure(figsize=(9, 7))

ax = fig.add_subplot(
    111,
    projection=wcs_12
)


####################################
# 12CO como fondo
####################################

im = ax.imshow(
    co12,
    origin="lower",
    cmap="inferno"
)


####################################
# Colorbar
####################################

cbar = plt.colorbar(
    im,
    ax=ax
)

cbar.set_label(
    r"$^{12}$CO(2-1) [W m$^{-2}$ sr$^{-1}$]"
)


####################################
# Contornos CII
####################################

contours = ax.contour(
    cii,
    levels=levels_cii_valid,
    colors="cyan",
    linewidths=1.2
)


####################################
# Etiquetas
####################################

ax.clabel(
    contours,
    inline=True,
    fontsize=9,
    fmt=lambda x: f"{x/rms_cii:.0f}$\\sigma$"
)


####################################
# Ejes
####################################

ax.set_xlabel("RA")
ax.set_ylabel("Dec")

ax.set_title(
    r"$^{12}$CO(2-1) + [CII] contours"
)

plt.tight_layout()
plt.show()


####################################
# MAPA 2
# 13CO + contornos CII
####################################

fig = plt.figure(figsize=(9, 7))

ax = fig.add_subplot(
    111,
    projection=wcs_13
)


####################################
# 13CO como fondo
####################################

im = ax.imshow(
    co13,
    origin="lower",
    cmap="inferno"
)


####################################
# Colorbar
####################################

cbar = plt.colorbar(
    im,
    ax=ax
)

cbar.set_label(
    r"$^{13}$CO(2-1) [W m$^{-2}$ sr$^{-1}$]"
)


####################################
# Contornos CII
####################################

contours = ax.contour(
    cii,
    levels=levels_cii_valid,
    colors="cyan",
    linewidths=1.2
)


####################################
# Etiquetas
####################################

ax.clabel(
    contours,
    inline=True,
    fontsize=9,
    fmt=lambda x: f"{x/rms_cii:.0f}$\\sigma$"
)



####################################
# Ejes
####################################

ax.set_xlabel("RA")
ax.set_ylabel("Dec")

ax.set_title(
    r"$^{13}$CO(2-1) + [CII] contours"
)

plt.tight_layout()
plt.show()


####################################
# Cerrar archivos
####################################

cii_hdu.close()
co12_hdu.close()
co13_hdu.close()