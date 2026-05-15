import numpy as np
from astropy.io import fits
from scipy import ndimage
import os

def crear_mascara_co(
    input_fits,
    output_dir,
    output_name="mask_CO.fits",
    threshold=None,
    sigma=None,
    snr=None,
    min_pixels=10,
    smooth_sigma=None
):
    """
    Genera una máscara binaria de emisión para un mapa de CO
    y la guarda en un directorio específico.

    Parámetros:
    -----------
    input_fits : str
        Archivo FITS de entrada.
    
    output_dir : str
        Carpeta donde se guardará el archivo.
    
    output_name : str
        Nombre del archivo de salida.
    
    threshold / sigma / snr :
        Definición del umbral.
    
    min_pixels : int
        Tamaño mínimo de estructuras.
    
    smooth_sigma : float
        Suavizado gaussiano opcional.
    """

    # --- crear carpeta si no existe ---
    os.makedirs(output_dir, exist_ok=True)

    # --- path completo ---
    output_path = os.path.join(output_dir, output_name)

    # --- leer datos ---
    data, header = fits.getdata(input_fits, header=True)
    data = np.nan_to_num(data, nan=0.0)

    # --- suavizado opcional ---
    if smooth_sigma is not None:
        data = ndimage.gaussian_filter(data, smooth_sigma)

    # --- definir threshold ---
    if snr is not None:
        if sigma is None:
            raise ValueError("Debes proporcionar sigma si usas snr.")
        threshold = snr * sigma

    if threshold is None:
        raise ValueError("Debes definir threshold o snr.")

    # --- máscara ---
    mask = data > threshold

    # --- limpiar regiones pequeñas ---
    labeled, num_features = ndimage.label(mask)
    sizes = ndimage.sum(mask, labeled, range(1, num_features + 1))

    clean_mask = np.zeros_like(mask)

    for i, size in enumerate(sizes):
        if size >= min_pixels:
            clean_mask[labeled == i + 1] = 1

    # --- guardar ---
    fits.writeto(output_path, clean_mask.astype(np.uint8), header, overwrite=True)

    print(f"Máscara guardada en: {output_path}")

crear_mascara_co(
    input_fits= r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_12CO21.fits",
    output_dir= r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC",
    output_name="12CO_mask.fits",
    snr=3,
    sigma=0.2,
    min_pixels=20
)

crear_mascara_co(
    input_fits= r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_13CO21.fits",
    output_dir=  r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC",
    output_name="13CO_mask.fits",
    snr=3,
    sigma=0.2,
    min_pixels=20
)