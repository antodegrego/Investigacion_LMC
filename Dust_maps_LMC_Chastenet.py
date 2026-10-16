import tarfile
import os
from astropy.io import fits

# ============================================================
# ARCHIVO .tar.gz
# ============================================================

tar_file = r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\lmc_themis_singleu.tar.gz"

# Carpeta donde vamos a extraer los archivos
outdir = r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\lmc_themis_singleu"

os.makedirs(outdir, exist_ok=True)


# ============================================================
# EXTRAER EL .tar.gz
# ============================================================

with tarfile.open(tar_file, "r:gz") as tar:

    for member in tar.getmembers():
        print(member.name)

    tar.extractall(outdir)


# ============================================================
# ARCHIVO QUE QUEREMOS INSPECCIONAR
# ============================================================

filename = os.path.join(
    outdir,
    "lmc_aj13_single_13aug15_dust_param_max.fits"
)


# ============================================================
# INSPECCIONAR LAS EXTENSIONES DEL FITS
# ============================================================

with fits.open(filename) as hdul:

    print("\nNúmero de extensiones:", len(hdul))

    for i, hdu in enumerate(hdul):

        print("\n-----------------------------")
        print("HDU:", i)

        print("EXTNAME:", hdu.header.get("EXTNAME"))
        print("BUNIT:", hdu.header.get("BUNIT"))
        print(
            "Shape:",
            hdu.data.shape if hdu.data is not None else None
        )