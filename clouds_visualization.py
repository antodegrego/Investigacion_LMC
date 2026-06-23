import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.patches import Ellipse

# -----------------------
# Leer mapa momento 0
# -----------------------
fitsfile = r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_12CO21.fits"

hdu = fits.open(fitsfile)[0]
data = hdu.data.squeeze()
header = hdu.header
wcs = WCS(header)

# -----------------------
# Leer catálogo
# -----------------------
cat = pd.read_csv(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_12CO21_mask13CO_corrected.props.csv"
)

# -----------------------
# Eliminar radios NaN
# -----------------------
cat = cat[np.isfinite(cat["RAD_PC"])]

print(f"Número final de nubes: {len(cat)}")

# -----------------------
# Escala angular del mapa
# -----------------------
pixscale_deg = np.abs(header["CDELT1"])
pixscale_rad = np.deg2rad(pixscale_deg)

# -----------------------
# Figura
# -----------------------
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(projection=wcs.celestial)

vmax = np.nanpercentile(data,99)

ax.imshow(
    data,
    origin="lower",
    cmap="inferno",
    vmax=vmax
)

# -----------------------
# Dibujar nubes
# -----------------------
for _, row in cat.iterrows():

    x = row["XCTR_PIX"]
    y = row["YCTR_PIX"]

    R_pc = row["RAD_PC"]
    D_pc = row["DISTANCE_PC"]

    theta_rad = R_pc / D_pc

    radius_pix = theta_rad / pixscale_rad

    ellipse = Ellipse(
        (x, y),
        width=2*radius_pix,
        height=2*radius_pix,
        angle=0,
        fill=False,
        lw=1,
        color="cyan"
    )

    ax.add_patch(ellipse)

ax.set_xlabel("RA")
ax.set_ylabel("Dec")

plt.tight_layout()
plt.show()