from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

####################################
# Cargar mapas
####################################

cii_hdu = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\line_flux_CII.fits")
co12_hdu = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_12CO_reprojected.fits")
co13_hdu = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_13CO_reprojected.fits")

cii = np.squeeze(cii_hdu[0].data)
co12 = np.squeeze(co12_hdu[0].data)
co13 = np.squeeze(co13_hdu[0].data)

header = cii_hdu[0].header.copy()

####################################
# Región de ruido para 13CO
####################################

# Ajustar índices según tu mapa
noise_region = co13[50:100,50:100]

rms_13 = np.nanstd(noise_region)

print("RMS 13CO =", rms_13)

####################################
# Máscara común
####################################

mask = (
    (co13 > 3*rms_13)
    & np.isfinite(co13)
    & np.isfinite(co12)
    & np.isfinite(cii)
    & (co12 > 0)
    & (cii > 0)
)

print("Pixeles válidos:", np.sum(mask))

####################################
# Mapas de razón
####################################

valor_negativo = -1e6

ratio_co12_cii = np.where(
    mask,
    co12/cii,
    valor_negativo
)

ratio_co13_cii = np.where(
    mask,
    co13/cii,
    valor_negativo
)

####################################
# Guardar mapas
####################################

header12 = header.copy()
header12["BUNIT"] = "12CO/CII"

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_ratio_12CO_CII.fits",
    ratio_co12_cii,
    header12,
    overwrite=True
)

header13 = header.copy()
header13["BUNIT"] = "13CO/CII"

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_ratio_13CO_CII.fits",
    ratio_co13_cii,
    header13,
    overwrite=True
)

####################################
# Mostrar mapas
####################################

#plt.figure(figsize=(8,6))
#plt.imshow(
#    ratio_cii_co12,
#    origin="lower",
#    cmap="inferno"
#)
#plt.colorbar(label="[CII]/12CO")
#plt.title("[CII]/12CO")
#plt.show()

#plt.figure(figsize=(8,6))
#plt.imshow(
#    ratio_cii_co13,
#    origin="lower",
#    cmap="inferno"
#)
#plt.colorbar(label="[CII]/13CO")
#plt.title("[CII]/13CO")
#plt.show()

####################################
# Scatter [CII] vs 12CO
####################################

co12_valid = co12[mask]
co13_valid = co13[mask]
cii_valid = cii[mask]

plt.figure(figsize=(6,6))

plt.scatter(
    cii_valid,
    co12_valid,
    s=2,
    alpha=0.3
)

plt.xlabel(r"$^{12}$CO [W m$^{-2}$ sr$^{-1}$]")
plt.ylabel(r"[CII] [W m$^{-2}$ sr$^{-1}$]")

plt.xscale("log")
plt.yscale("log")

plt.grid(alpha=0.5)

####################################
# Fiteo
####################################

coef12 = np.polyfit(
    np.log10(cii_valid),
    np.log10(co12_valid),
    1
)

m12, b12 = coef12

xfit = np.linspace(
    np.min(np.log10(cii_valid)),
    np.max(np.log10(cii_valid)),
    1000
)

yfit = m12*xfit + b12

plt.plot(
    10**xfit,
    10**yfit,
    'r',
    lw=2,
    label=fr'$\log(\mathrm{{CII}})={m12:.3f}\log(\mathrm{{CO}})+{b12:.3f}$'
)

plt.legend()
plt.show()

####################################
# Scatter [CII] vs 13CO
####################################

plt.figure(figsize=(6,6))

plt.scatter(
    cii_valid,
    co13_valid,
    s=2,
    alpha=0.3
)

plt.xlabel(r"$^{13}$CO [W m$^{-2}$ sr$^{-1}$]")
plt.ylabel(r"[CII] [W m$^{-2}$ sr$^{-1}$]")

plt.xscale("log")
plt.yscale("log")

plt.grid(alpha=0.5)

coef13 = np.polyfit(
    np.log10(cii_valid),
    np.log10(co13_valid),
    1
)

m13, b13 = coef13

xfit = np.linspace(
    np.min(np.log10(cii_valid)),
    np.max(np.log10(cii_valid)),
    1000
)

yfit = m13*xfit + b13

plt.plot(
    10**xfit,
    10**yfit,
    'r',
    lw=2,
    label=fr'$\log(\mathrm{{CII}})={m13:.3f}\log(\mathrm{{CO}})+{b13:.3f}$'
)

plt.legend()
plt.show()

####################################
# Razones globales
####################################

R12 = np.nansum(co12_valid) / np.nansum(cii_valid)
R13 = np.nansum(co13_valid) / np.nansum(cii_valid)

print()
print("12CO/[CII] =", R12)
print("13CO/[CII] =", R13)