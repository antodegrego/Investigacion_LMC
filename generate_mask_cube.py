from astropy.io import fits
import numpy as np
from scipy.ndimage import label

# ==========================
# Leer cubo
# ==========================

cubefile = r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMolR_aca_13CO21_K_cube.fits"

data = fits.getdata(cubefile, memmap=True)
header = fits.getheader(cubefile)

print("Shape:", data.shape)

# ==========================
# RMS local
# ==========================

# Usa canales sin emisión
noise_map = np.nanstd(data[:20], axis=0)

median_rms = np.nanmedian(noise_map)

print("RMS mediano =", median_rms)

# ==========================
# Eliminar bordes muy ruidosos
# ==========================

# Ajustar entre 1.5 y 3 según resultados
spatial_mask = noise_map < 2.0 * median_rms

# ==========================
# Detección 3 sigma local
# ==========================

above = data > (3.0 * noise_map[np.newaxis, :, :])

# Aplicar máscara espacial
above &= spatial_mask[np.newaxis, :, :]

# ==========================
# Exigir 3 canales consecutivos
# ==========================

mask = (
    above[:-2] &
    above[1:-1] &
    above[2:]
)

mask_full = np.zeros_like(data, dtype=np.uint8)

mask_full[1:-1] = mask

# ==========================
# Eliminar estructuras pequeñas
# ==========================

labels, nlab = label(mask_full)

clean_mask = np.zeros_like(mask_full, dtype=np.uint8)

for i in range(1, nlab + 1):

    region = labels == i

    # mínimo número de vóxeles
    if region.sum() >= 50:
        clean_mask[region] = 1

print("Número de regiones:", nlab)

# ==========================
# Guardar
# ==========================

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMolR_aca_13CO21_mask_clean.fits",
    clean_mask,
    header,
    overwrite=True
)

print("Máscara guardada.")