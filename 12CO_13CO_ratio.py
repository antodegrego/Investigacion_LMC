from astropy.io import fits
import numpy as np

hdu12 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_12CO21.fits")
hdu13 = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_13CO21.fits")

data12 = hdu12[0].data
data13 = hdu13[0].data
header = hdu12[0].header

# Valor muy negativo para los píxeles que no cumplen
valor_negativo = -1e6

# definir una región de ruido
noise_region = data13[846:902,317:440]
rms_mom0 = np.nanstd(noise_region)

#Aplicamos una mascara
mask =  np.isfinite(data13) & (data13 > 0) #& (data12 > 0) & np.isfinite(data12)
mask = (data13 > 3*rms_mom0) & np.isfinite(data13) & np.isfinite(data12) 
ratio = np.where(mask, data13/data12, valor_negativo)
ratio[data12==0] = valor_negativo #np.nan

# asegurar que sea 2D
ratio = np.squeeze(ratio)

# corregir header
header["NAXIS"] = 2
header["NAXIS1"] = ratio.shape[1]
header["NAXIS2"] = ratio.shape[0]

# eliminar keywords de ejes extra si existen
for key in ["NAXIS3","CTYPE3","CRVAL3","CRPIX3","CDELT3","CUNIT3"]:
    if key in header:
        del header[key]

print(ratio.shape)
print(hdu12[0].header["NAXIS"])
print(np.nanmin(ratio))
print(np.nanmax(ratio))

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_ratio_13CO21_12CO21_13mask_only.fits",
    ratio,
    header,
    overwrite=True)


import matplotlib.pyplot as plt

plt.imshow(ratio, origin="lower", cmap="inferno", vmin=0, vmax=0.1)
plt.colorbar(label="13CO(2-1) / 12CO(2-1)")
plt.title("Line Ratio Map (13CO(2-1)/12CO(2-1))")
plt.show()

#################################
#Graficar cociente entre 13C
# O y 12CO
###################################

import matplotlib.pyplot as plt
import numpy as np

# asegurar misma forma
data12 = np.squeeze(data12)
data13 = np.squeeze(data13)



# máscara más estricta (recomendado)
mask = (data12 > 0) & (data13 > 0) & np.isfinite(data12) & np.isfinite(data13)
mask = (data13 > 0) & np.isfinite(data13)
mask = (data13 > 3*rms_mom0) & np.isfinite(data13) & np.isfinite(data12) 

# extraer valores válidos
co12 = data12[mask]
co13 = data13[mask]

# scatter plot
plt.figure(figsize=(6,6))
plt.scatter(co12, co13, s=1, alpha=0.3)

plt.xlabel("12CO(2-1)")
plt.ylabel("13CO(2-1)")
plt.title("13CO vs 12CO")

# límites de ejes
plt.xlim(-0.5, np.max(co12))
plt.ylim(-0.5, np.max(co13))

# escala log
#plt.xscale("log")
#plt.yscale("log")

plt.grid(True, which="both", ls="--", alpha=0.5)
plt.show()