from astropy.io import fits
import numpy as np

hdu = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMolR_aca_13CO21_K_cube.fits", memmap=True)

cube = hdu[0].data

print("shape =", cube.shape)

print("pixeles validos =", np.sum(np.isfinite(cube)))

noise_cube = cube[:, 580:680, 230:320]

print(noise_cube.shape)

print(np.sum(np.isfinite(noise_cube)))

rms_chan = np.nanstd(noise_cube)

print(rms_chan)

sigma_chan = rms_chan
dv = 0.25
N = 320

sigma_mom0 = sigma_chan * dv * np.sqrt(N)

print(sigma_mom0)