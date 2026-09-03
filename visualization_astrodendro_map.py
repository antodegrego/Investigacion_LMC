import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.patches import Ellipse


# --------------------------------------------------
# Leer mapa momento 0
# --------------------------------------------------

fitsfile = 'moment0_12CO21.fits'

hdu = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_13CO21.fits")[0]

data = hdu.data.squeeze()
header = hdu.header

wcs = WCS(header)


# --------------------------------------------------
# Leer catálogo astrodendro
# --------------------------------------------------

cat = pd.read_csv(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\catalog_13CO_resolved.csv")

print("Número de estructuras:", len(cat))


# --------------------------------------------------
# Escala angular del mapa
# --------------------------------------------------

pixscale_deg = np.abs(header['CDELT1'])
pixscale_arcsec = pixscale_deg * 3600.0

print("Escala de píxel:", pixscale_arcsec, "arcsec/pixel")


# --------------------------------------------------
# Figura
# --------------------------------------------------

fig = plt.figure(figsize=(10, 10))

ax = plt.subplot(projection=wcs.celestial)

vmax = np.nanpercentile(data, 99)

ax.imshow(
    data,
    origin='lower',
    cmap='inferno',
    vmax=vmax
)


# --------------------------------------------------
# Dibujar estructuras
# --------------------------------------------------

for _, row in cat.iterrows():

    # Centroide en coordenadas celestes
    ra = row['x_cen']
    dec = row['y_cen']

    # Radio de la estructura en arcsec
    radius_arcsec = row['radius']

    # Convertir RA/Dec -> píxeles
    x, y = wcs.celestial.world_to_pixel_values(ra, dec)

    # Radio de la estructura en píxeles
    radius_pix = radius_arcsec / pixscale_arcsec

    ellipse = Ellipse(
        (x, y),
        width=2 * radius_pix,
        height=2 * radius_pix,
        angle=0,
        fill=False,
        lw=1,
        color='cyan'
    )

    ax.add_patch(ellipse)


# --------------------------------------------------
# Etiquetas
# --------------------------------------------------

ax.set_xlabel('RA')
ax.set_ylabel('Dec')

plt.tight_layout()
plt.show()


