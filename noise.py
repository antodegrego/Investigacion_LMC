from astropy.io import fits
import numpy as np

cubefile = 'NMolR_aca_13CO21_K_cube.fits'

# abrir cubo original
hdul = fits.open(cubefile)

data = hdul[0].data
header = hdul[0].header

# crear cubo con ruido constante
noise_data = np.full(data.shape, 0.2, dtype=np.float32)

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\noise_cube13.fits",
    noise_data,
    header,
    overwrite=True
)

print("noise_cube.fits creado")