from astropy.io import fits
import numpy as np

# Leer cubo
hdul = fits.open(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMolR_aca_13CO21_K_cube.fits",
    memmap=True
)

data = hdul[0].data
header = fits.getheader(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMolR_aca_13CO21_K_cube.fits")

# Calcular RMS usando canales sin emisión
rms = np.nanstd(data[:20])

print(f"RMS = {rms:.4f}")

# Máscara 3 sigma
mask = (data > 3*rms).astype(np.uint8)

# Guardar
fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMolR_aca_13CO21_mask_3sigma.fits",
    mask,
    header,
    overwrite=True
)

print("Máscara guardada.")