from astropy.io import fits
import numpy as np

archivo = r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\12CO21.asgn.fits"

d = fits.getdata(archivo)

print("Mínimo:", np.min(d))
print("Máximo:", np.max(d))
print("Primeros valores únicos:")
print(np.unique(d)[:20])

print("Número de estructuras:", len(np.unique(d))-1)

import numpy as np

canales = np.where(np.sum(d > 0, axis=(1,2)) > 0)[0]

print("Canales con estructuras:")
print(canales)
print("Primer canal útil:", canales[0])

import matplotlib.pyplot as plt

k = canales[0]

plt.figure(figsize=(8,8))
plt.imshow(d[k,:,:] > 0, origin='lower')
plt.colorbar()
plt.title(f'Canal {k}')
plt.show()

mascara = np.max(d > 0, axis=0)

plt.figure(figsize=(8,8))
plt.imshow(mascara, origin='lower')
plt.colorbar()
plt.title("Todas las estructuras detectadas")
plt.show()

proy = np.max(d, axis=0)

plt.figure(figsize=(8,8))
plt.imshow(proy, origin='lower')
plt.colorbar(label='ID estructura')
plt.show()

###########################################################################


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

mom0 = fits.getdata(
r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_12CO21.fits"
)

asgn = fits.getdata(
r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\12CO21.asgn.fits"
)

# Proyección 2D de estructuras
mask2D = np.max(asgn > 0, axis=0)
print(mask2D.shape)
print(np.sum(mask2D))

plt.figure(figsize=(10,8))

# momento cero con contraste más razonable
plt.imshow(
    mom0,
    origin='lower',
    cmap='gray',
    vmin=np.nanpercentile(mom0,5),
    vmax=np.nanpercentile(mom0,99)
)

# contornos visibles
plt.contour(
    mask2D.astype(float),
    levels=[0.5],
    colors='red',
    linewidths=1.0
)

plt.show()



from astropy.io import fits

hdul = fits.open(
r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\12CO21.props.fits"
)

hdul.info()

from astropy.table import Table

tabla = Table.read(
r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\12CO21.props.fits",
hdu=1
)

print(tabla.colnames)

print(tabla[:10])


import pandas as pd

# Leer archivo
df = pd.read_csv(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_12CO21.props.csv")

# Seleccionar columnas de interés
tabla = df[["CLOUDNUM", "RAD_PC", "SIGV_KMS"]]

print(tabla)