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

header = cii_hdu[0].header

wcs = WCS(header)


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
# Revisar valores
####################################

print()
print("================================")
print("VALORES DE LOS MAPAS")
print("================================")

print("Máximo CII  =", np.nanmax(cii))
print("Máximo 12CO =", np.nanmax(co12))
print("Máximo 13CO =", np.nanmax(co13))


####################################
# RMS
####################################

# IMPORTANTE:
# Estos RMS se calculan directamente sobre los mapas
# que están en W/m2/sr.

# Región de ruido
noise_region_12 = co12[160:180, 90:120]
noise_region_13 = co13[160:180, 90:120]

# RMS
rms_12 = np.nanstd(noise_region_12)
rms_13 = np.nanstd(noise_region_13)


print()
print("================================")
print("RMS")
print("================================")

print(
    f"RMS 12CO = {rms_12:.4e} W/m2/sr"
)

print(
    f"RMS 13CO = {rms_13:.4e} W/m2/sr"
)


####################################
# Niveles de contorno
####################################

# Multiplicadores de RMS

levels_xrms_12 = np.array([
    3, 6, 9, 12
])

levels_xrms_13 = np.array([
    3, 6, 9, 12
])


# Convertir a unidades físicas

levels_12 = levels_xrms_12 * rms_12
levels_13 = levels_xrms_13 * rms_13


####################################
# Mostrar niveles
####################################

print()
print("================================")
print("NIVELES 12CO")
print("================================")

for n, level in zip(levels_xrms_12, levels_12):

    print(
        f"{n:3.0f} RMS = {level:.4e} W/m2/sr"
    )


print()
print("================================")
print("NIVELES 13CO")
print("================================")

for n, level in zip(levels_xrms_13, levels_13):

    print(
        f"{n:3.0f} RMS = {level:.4e} W/m2/sr"
    )


####################################
# Comprobar qué niveles existen
# realmente en los mapas
####################################

levels_12_valid = levels_12[
    levels_12 <= np.nanmax(co12)
]

levels_13_valid = levels_13[
    levels_13 <= np.nanmax(co13)
]


print()
print("================================")
print("NIVELES UTILIZABLES")
print("================================")

print("12CO:")
print(levels_12_valid)

print()
print("13CO:")
print(levels_13_valid)


####################################
# MAPA 1
# CII + 12CO
####################################

fig = plt.figure(figsize=(9, 7))

ax = fig.add_subplot(
    111,
    projection=wcs
)


####################################
# CII como fondo
####################################

im = ax.imshow(
    cii,
    origin="lower",
    cmap="inferno"
)


####################################
# Barra de color
####################################

cbar = plt.colorbar(
    im,
    ax=ax
)

cbar.set_label(
    r"[CII] [W m$^{-2}$ sr$^{-1}$]"
)


####################################
# Contornos 12CO
####################################

contours = ax.contour(
    co12,
    levels=levels_12_valid,
    colors="cyan",
    linewidths=1.2
)


####################################
# Etiquetas de los contornos
####################################

ax.clabel(
    contours,
    inline=True,
    fontsize=9,
    fmt=lambda x: f"{x/rms_12:.0f}$\\sigma$"
)


####################################
# Ejes
####################################

ax.set_xlabel("RA")
ax.set_ylabel("Dec")


####################################
# Título
####################################

ax.set_title(
    r"[CII] moment 0 + $^{12}$CO(2-1)"
)


plt.tight_layout()
plt.show()


####################################
# MAPA 2
# CII + 13CO
####################################

fig = plt.figure(figsize=(9, 7))

ax = fig.add_subplot(
    111,
    projection=wcs
)


####################################
# CII como fondo
####################################

im = ax.imshow(
    cii,
    origin="lower",
    cmap="inferno"
)


####################################
# Barra de color
####################################

cbar = plt.colorbar(
    im,
    ax=ax
)

cbar.set_label(
    r"[CII] [W m$^{-2}$ sr$^{-1}$]"
)


####################################
# Contornos 13CO
####################################

contours = ax.contour(
    co13,
    levels=levels_13_valid,
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
    fmt=lambda x: f"{x/rms_13:.0f}$\\sigma$"
)


####################################
# Ejes
####################################

ax.set_xlabel("RA")
ax.set_ylabel("Dec")


####################################
# Título
####################################

ax.set_title(
    r"[CII] moment 0 + $^{13}$CO(2-1)"
)


plt.tight_layout()
plt.show()


####################################
# Cerrar archivos
####################################

cii_hdu.close()
co12_hdu.close()
co13_hdu.close()