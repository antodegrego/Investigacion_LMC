from astropy.io import fits
from astropy import units as u
from astropy.constants import k_B, c

# Archivo de entrada
hdu12 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_12CO21.fits")
hdu13 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_13CO21.fits")
input_fits = "mapa_original.fits"

# Factor por el que se multiplicará el mapa
nu_12 = 230.538 * u.GHz   # 12CO(2-1)
nu_13 = 220.400 * u.GHz   # 13CO(2-1)

factor = (
    2 * k_B * nu_12**3 / c**3
    * (1*u.km/u.s)
).to(u.W/u.m**2/u.K)


# Abrir el FITS
data = hdu12[0].data
header = hdu12[0].header.copy()

# Multiplicar todos los píxeles
data_new = data * factor

# Actualizar unidades
header['BUNIT'] = 'W/m2/Sr'

# Guardar en un nuevo FITS
fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_moment0_12CO21_Wm2sr.fits",
    data_new,
    header,
    overwrite=True
)
print("Mapa guardado")