import numpy as np
from astropy.io import fits
from astropy.constants import c


# ============================================================
# ARCHIVOS
# ============================================================

archivo_hi = "HI4PI.fits"
archivo_co = "NMolR_aca_12CO21_K_cube.fits"
archivo_cii = "F0842_FI_IFS_0900361_RED_WXY_100010-300230.fits"


# ============================================================
# RANGO DE VELOCIDAD QUE NOS INTERESA
# ============================================================

vmin = 190.0       # km/s
vmax = 290.0       # km/s


# ============================================================
# FUNCIONES
# ============================================================

def velocity_axis_from_header(header):
    """
    Construye el eje espectral para un cubo cuyo eje 3
    está expresado directamente en velocidad.
    """

    nchan = header["NAXIS3"]

    pix = np.arange(1, nchan + 1)

    velocity = (
        header["CRVAL3"]
        + (pix - header["CRPIX3"]) * header["CDELT3"]
    )

    # Convertir m/s -> km/s si corresponde
    if header.get("CUNIT3", "").strip().lower() in ["m/s", "m s-1"]:
        velocity /= 1000.0

    return velocity


def cii_velocity_axis(header):
    """
    Convierte el eje de longitud de onda del cubo [CII]
    a velocidad usando la longitud de onda de reposo.
    """

    nchan = header["NAXIS3"]

    pix = np.arange(1, nchan + 1)

    wavelength = (
        header["CRVAL3"]
        + (pix - header["CRPIX3"]) * header["CDELT3"]
    )

    # Longitud de onda de reposo de [CII] 158 um
    lambda0 = 157.7409   # um

    # Velocidad Doppler aproximada
    velocity = (
        (wavelength - lambda0) / lambda0
        * c.to("km/s").value
    )

    return wavelength, velocity


def channel_limits(velocity, vmin, vmax):

    mask = (velocity >= vmin) & (velocity <= vmax)

    indices = np.where(mask)[0]

    if len(indices) == 0:
        raise ValueError(
            f"No hay canales entre {vmin} y {vmax} km/s"
        )

    return indices[0], indices[-1], mask


# ============================================================
# 1. H I
# ============================================================

print("\n========================================")
print("H I")
print("========================================")

with fits.open(archivo_hi, memmap=True) as hdul:

    header_hi = hdul[0].header
    data_hi = hdul[0].data

    velocity_hi = velocity_axis_from_header(header_hi)

    i0_hi, i1_hi, mask_hi = channel_limits(
        velocity_hi,
        vmin,
        vmax
    )

    print("Shape original:", data_hi.shape)

    print(
        f"Velocidad original: "
        f"{velocity_hi.min():.2f} - "
        f"{velocity_hi.max():.2f} km/s"
    )

    print(
        f"Canales seleccionados: "
        f"{i0_hi} - {i1_hi}"
    )

    print(
        f"Velocidad seleccionada: "
        f"{velocity_hi[i0_hi]:.2f} - "
        f"{velocity_hi[i1_hi]:.2f} km/s"
    )

    # Copiamos solamente los canales necesarios
    data_hi_cut = np.asarray(
        data_hi[i0_hi:i1_hi + 1],
        dtype=np.float32
    )

    header_hi_cut = header_hi.copy()

    # Actualizar referencia del eje espectral
    header_hi_cut["CRPIX3"] -= i0_hi
    header_hi_cut["NAXIS3"] = data_hi_cut.shape[0]

    fits.PrimaryHDU(
        data=data_hi_cut,
        header=header_hi_cut
    ).writeto(
        "HI4PI_190_290kms.fits",
        overwrite=True
    )


# ============================================================
# 2. 12CO
# ============================================================

print("\n========================================")
print("12CO")
print("========================================")

with fits.open(archivo_co, memmap=True) as hdul:

    header_co = hdul[0].header
    data_co = hdul[0].data

    velocity_co = velocity_axis_from_header(header_co)

    i0_co, i1_co, mask_co = channel_limits(
        velocity_co,
        vmin,
        vmax
    )

    print("Shape original:", data_co.shape)

    print(
        f"Velocidad original: "
        f"{velocity_co.min():.2f} - "
        f"{velocity_co.max():.2f} km/s"
    )

    print(
        f"Canales seleccionados: "
        f"{i0_co} - {i1_co}"
    )

    print(
        f"Velocidad seleccionada: "
        f"{velocity_co[i0_co]:.2f} - "
        f"{velocity_co[i1_co]:.2f} km/s"
    )

    # IMPORTANTE:
    # no cargamos todo el cubo a RAM.
    data_co_cut = np.asarray(
        data_co[i0_co:i1_co + 1],
        dtype=np.float32
    )

    header_co_cut = header_co.copy()

    header_co_cut["CRPIX3"] -= i0_co
    header_co_cut["NAXIS3"] = data_co_cut.shape[0]

    fits.PrimaryHDU(
        data=data_co_cut,
        header=header_co_cut
    ).writeto(
        "12CO21_190_290kms.fits",
        overwrite=True
    )


# ============================================================
# 3. [CII]
# ============================================================

print("\n========================================")
print("[CII]")
print("========================================")

with fits.open(archivo_cii, memmap=True) as hdul:

    # Buscar la extensión FLUX
    hdu_cii = None

    for hdu in hdul:
        if hdu.header.get("EXTNAME", "").strip() == "FLUX":
            hdu_cii = hdu
            break

    if hdu_cii is None:
        raise ValueError(
            "No se encontró la extensión FLUX del cubo [CII]."
        )

    header_cii = hdu_cii.header
    data_cii = hdu_cii.data

    wavelength_cii, velocity_cii = cii_velocity_axis(
        header_cii
    )

    i0_cii, i1_cii, mask_cii = channel_limits(
        velocity_cii,
        vmin,
        vmax
    )

    print("Shape original:", data_cii.shape)

    print(
        f"Longitud de onda: "
        f"{wavelength_cii.min():.5f} - "
        f"{wavelength_cii.max():.5f} um"
    )

    print(
        f"Velocidad equivalente: "
        f"{velocity_cii.min():.2f} - "
        f"{velocity_cii.max():.2f} km/s"
    )

    print(
        f"Canales seleccionados: "
        f"{i0_cii} - {i1_cii}"
    )

    print(
        f"Velocidad seleccionada: "
        f"{velocity_cii[i0_cii]:.2f} - "
        f"{velocity_cii[i1_cii]:.2f} km/s"
    )

    data_cii_cut = np.asarray(
        data_cii[i0_cii:i1_cii + 1],
        dtype=np.float32
    )

    header_cii_cut = header_cii.copy()

    header_cii_cut["CRPIX3"] -= i0_cii
    header_cii_cut["NAXIS3"] = data_cii_cut.shape[0]

    fits.PrimaryHDU(
        data=data_cii_cut,
        header=header_cii_cut
    ).writeto(
        "CII_190_290kms.fits",
        overwrite=True
    )


print("\n========================================")
print("PROCESO TERMINADO")
print("========================================")

print("\nArchivos creados:")
print("  HI4PI_190_290kms.fits")
print("  12CO21_190_290kms.fits")
print("  CII_190_290kms.fits")