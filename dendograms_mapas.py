from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astrodendro import Dendrogram
from astropy.stats import mad_std
from astropy.table import Table



def analizar_mapa(path_fits, label, output_fits):
    
    # --- cargar datos ---
    hdu = fits.open(path_fits)
    data = np.squeeze(hdu[0].data)
    header = hdu[0].header

    print(f"{label} shape:", data.shape)

    # --- máscara válida (sin tocar los datos) ---
    mask = np.isfinite(data)

    # --- estimar ruido SOLO en fondo ---
    noise_region = data[mask & (data < np.percentile(data[mask], 30))]
    sigma = mad_std(noise_region)

    print(f"{label} sigma:", sigma)

    # --- filtrar fondo ---
    #data_masked = np.where(data > 3*sigma, data, np.nan)

    # --- dendrograma ---
    d = Dendrogram.compute(
        data,
        min_value = 3 * sigma,
        min_delta = 0.2 * sigma,
        min_npix = 50
    )

    # --- crear catálogo ---

    catalog = []

    for i, leaf in enumerate(d.leaves):

        mask = leaf.get_mask()
        values = data[mask]

        # --- propiedades básicas ---
        npix = np.sum(mask)
        flux = np.nansum(values)
        peak = np.nanmax(values)

        # centroide (en pixeles)
        y, x = np.where(mask)
        x_mean = np.mean(x)
        y_mean = np.mean(y)

        # radio efectivo (asumiendo área circular)
        r_eff = np.sqrt(npix / np.pi)

        catalog.append([
            i+1, x_mean, y_mean, npix, flux, peak, r_eff
        ])

    #--- crear tabla ---
    table = Table(
        rows=catalog,
        names=("id", "x", "y", "npix", "flux", "peak", "r_eff")
        )

    # --- labels ---
    id_map = np.zeros_like(data, dtype=int)

    for i, leaf in enumerate(d.leaves):
        leaf_mask = leaf.get_mask()
        id_map[leaf_mask] = i + 1

    fits.writeto(
        output_fits.replace(".fits", "_labels.fits"),
        id_map,
        header,
        overwrite=True
    )


    # --- guardar catalogo ---
    catalog_name = f"catalog_{label}.fits"

    table.write(catalog_name, overwrite=True)

    # --- plot ---
    fig, ax = plt.subplots()

    im = ax.imshow(data, origin="lower", cmap="inferno", vmin=0, vmax=10)
    plt.colorbar(im, ax=ax, label=label)

    # contornos jerárquicos (esto es lo importante)
    levels = np.linspace(3*sigma, np.nanmax(data), 40)

    ax.contour(
        data,
        levels=levels,
        colors='cyan',
        linewidths=0.5
    )

    ax.set_title(f"Estructura en {label}")

    plt.show()


    # --- guardar mapa ---
    #fits.writeto(
    #    output_fits,
    #    data.astype("float32"),
    #    header,
    #    overwrite=True
    #)

    plt.show()

    return d

d12 = analizar_mapa(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_12CO21.fits",
    "12CO(2-1)",
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\clumps_12CO21.fits"
)

d13 = analizar_mapa(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_13CO21.fits",
    "13CO(2-1)",
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\clumps_13CO21.fits"
)