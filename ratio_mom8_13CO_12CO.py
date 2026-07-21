from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom

# ==========================================================
# Leer mapas momento 8
# ==========================================================

hdu12 = fits.open(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment8_12CO21.fits"
)

hdu13 = fits.open(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment8_13CO21.fits"
)

data12 = np.squeeze(hdu12[0].data)
data13 = np.squeeze(hdu13[0].data)

header = hdu12[0].header.copy()

valor_negativo = -1e6

# ==========================================================
# Rebinning (7 pix/beam -> 2 pix/beam)
# ==========================================================

factor = 2/7

data12 = zoom(data12, factor, order=1)
data13 = zoom(data13, factor, order=1)

# ==========================================================
# Calcular RMS del mapa rebineado
# ==========================================================

y0, y1 = int(846*factor), int(902*factor)
x0, x1 = int(317*factor), int(440*factor)

noise_region = data13[y0:y1, x0:x1]

rms = np.nanstd(noise_region)

print(f"RMS = {rms:.5f}")

# ==========================================================
# Máscara
# ==========================================================

mask = (
    (data13 > 3*rms) &
    np.isfinite(data12) &
    np.isfinite(data13) &
    (data12 != 0)
)

# ==========================================================
# Ratio
# ==========================================================

ratio = np.full(data12.shape, valor_negativo, dtype=float)

ratio[mask] = data13[mask]/data12[mask]

# ==========================================================
# Actualizar header
# ==========================================================

header["NAXIS"] = 2
header["NAXIS1"] = ratio.shape[1]
header["NAXIS2"] = ratio.shape[0]

header["CDELT1"] /= factor
header["CDELT2"] /= factor

for key in [
    "NAXIS3","CTYPE3","CRVAL3",
    "CRPIX3","CDELT3","CUNIT3"
]:
    if key in header:
        del header[key]

# ==========================================================
# Guardar FITS
# ==========================================================

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_ratio_M8_13CO21_12CO21_rebinned.fits",
    ratio,
    header,
    overwrite=True
)

# ==========================================================
# Mostrar mapa
# ==========================================================

plt.figure(figsize=(8,6))

plt.imshow(
    ratio,
    origin="lower",
    cmap="inferno",
    vmin=0,
    vmax=0.1
)

plt.colorbar(label="Moment 8 Ratio (13CO / 12CO)")

plt.title("Moment 8 Ratio Map")

plt.tight_layout()

plt.show()

# ==========================================================
# Scatter
# ==========================================================

co12 = data12[mask]
co13 = data13[mask]

plt.figure(figsize=(6,6))

plt.scatter(
    co12,
    co13,
    s=1,
    alpha=0.3
)

plt.xlabel("Moment 8 12CO(2-1)")
plt.ylabel("Moment 8 13CO(2-1)")

plt.grid(alpha=0.4)

plt.tight_layout()

plt.show()

# ==========================================================
# Ajuste lineal
# ==========================================================

coef = np.polyfit(co12, co13, 1)

m, b = coef

xfit = np.linspace(
    np.min(co12),
    np.max(co12),
    1000
)

yfit = np.polyval(coef, xfit)

plt.figure(figsize=(6,6))

plt.scatter(
    co12,
    co13,
    s=1,
    alpha=0.3
)

plt.plot(
    xfit,
    yfit,
    color="red",
    lw=2,
    label=fr"$y={m:.5f}x+{b:.5f}$"
)

plt.legend()

plt.xlabel("Moment 8 12CO(2-1)")
plt.ylabel("Moment 8 13CO(2-1)")

plt.grid(alpha=0.4)

plt.tight_layout()

plt.show()

print(f"Pendiente = {m:.5f}")
print(f"Intercepto = {b:.5f}")