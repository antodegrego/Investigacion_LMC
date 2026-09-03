import numpy as np
from astropy.table import Table

cat12 = Table.read(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\catalog_12CO.fits")
cat13 = Table.read(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\catalog_13CO.fits")

print(cat12.colnames)
print(cat13.colnames)

print('12CO:', len(cat12))
print('13CO:', len(cat13))

cat12[
    '_idx',
    'x_cen',
    'y_cen',
    'v_cen',
    'flux',
    'v_rms',
    'major_sigma',
    'minor_sigma',
    'position_angle'
].pprint(max_lines=30)

cat13[
    '_idx',
    'x_cen',
    'y_cen',
    'v_cen',
    'flux',
    'v_rms',
    'major_sigma',
    'minor_sigma',
    'position_angle'
].pprint(max_lines=30)

import pandas as pd

df12 = cat12.to_pandas()
df13 = cat13.to_pandas()

df12['v_rms']
df12['flux']

cat12.write(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\catalog_12CO.csv",
    format='csv',
    overwrite=True
)

cat13.write(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\catalog_13CO.csv",
    format='csv',
    overwrite=True
)


# -----------------------
# Leer catálogo
# -----------------------

cat = Table.read(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\catalog_13CO.fits")

# -----------------------
# Número inicial
# -----------------------

n_inicial = len(cat)

print("Número inicial de estructuras:", n_inicial)

# -----------------------
# Resolución angular
# -----------------------

resolution = 7.0   # arcsec

# -----------------------
# Eliminar NaN
# -----------------------

mask_finite = np.isfinite(cat['radius'])

cat = cat[mask_finite]

print("Después de eliminar radius NaN:", len(cat))

# -----------------------
# Eliminar estructuras
# más pequeñas que la resolución
# -----------------------

cat_resolved = cat[cat['radius'] >= resolution]

print("Número final de estructuras:", len(cat_resolved))
print("Estructuras eliminadas:", n_inicial - len(cat_resolved))

cat_resolved.write(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\catalog_13CO_resolved.csv",
    format='csv',
    overwrite=True
)