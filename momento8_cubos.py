from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

# Leer cubo
cube = SpectralCube.read(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMolR_aca_13CO21_K_cube.fits"
)

cube = cube.with_spectral_unit(u.km/u.s)

print(cube.spectral_axis.min())
print(cube.spectral_axis.max())


# Región donde está la línea
subcube = cube.spectral_slab(201*u.km/u.s, 279*u.km/u.s)

# Región de ruido
noise_cube = cube.spectral_slab(1*u.km/u.s, 55*u.km/u.s)

print(noise_cube.shape)
print(subcube.shape)

# RMS global
rms = np.nanstd(noise_cube.unmasked_data[:].value)

print(f"RMS = {rms:.3f} K")

# Máscara de 3 sigma
masked_cube = subcube.with_mask(subcube > 3*rms*u.K)


# Tpeak (momento 8)
Tpeak = cube.max(axis=0, how='slice')

# Guardar FITS
Tpeak.write(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment8_13CO21.fits",
    overwrite=True
)

# Mostrar
plt.figure(figsize=(8,6))
plt.imshow(Tpeak.value, origin="lower")
plt.colorbar(label=str(Tpeak.unit))
plt.title("13CO(2-1) Tpeak (Moment 8)")
plt.show()