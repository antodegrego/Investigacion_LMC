from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u

cube = SpectralCube.read(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMolR_aca_12CO21_K_cube.fits")

cube = cube.with_spectral_unit(u.km/u.s)


print(cube.shape)
print(cube.spectral_axis.min())
print(cube.spectral_axis.max())
subcube = cube.spectral_slab(201*u.km/u.s, 279*u.km/u.s)
moment0 = subcube.moment(order=0)
moment0.write(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_12CO21.fits", overwrite=True)

plt.imshow(moment0.value, origin="lower")
plt.colorbar(label=str(moment0.unit))
plt.title("12CO21 Moment 0")
plt.show()

#spectrum = cube.mean(axis=(1,2), how="slice")

#plt.plot(cube.spectral_axis.value, spectrum.value)
#plt.xlabel("Velocity (km/s)")
#plt.ylabel("Intensity")
#plt.show()

cube2 = SpectralCube.read(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMolR_aca_13CO21_K_cube.fits")

cube2 = cube2.with_spectral_unit(u.km/u.s)

print(cube2.shape)
print(cube2.spectral_axis.min())
print(cube2.spectral_axis.max())
subcube2 = cube2.spectral_slab(201*u.km/u.s, 279*u.km/u.s)
moment0_13CO = subcube2.moment(order=0)
moment0_13CO.write(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_13CO21.fits", overwrite=True)

plt.imshow(moment0_13CO.value, origin="lower")
plt.colorbar(label=str(moment0_13CO.unit))
plt.title("13CO21 Moment 0")
plt.show()