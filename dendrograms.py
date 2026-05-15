from astropy.io import fits
import numpy as np

hdu = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\ratio_12CO21_13CO21.fits")
data = hdu[0].data
header = hdu[0].header
print(data.shape)

data = np.nan_to_num(data, nan=0)

sigma = np.nanstd(data[0:100,0:100])

from astrodendro import Dendrogram

d = Dendrogram.compute(
    data,
    min_value = 3*sigma,
    min_delta = 2*sigma,
    min_npix = 20
)

import matplotlib.pyplot as plt

fig, ax = plt.subplots()

im = ax.imshow(data, origin="lower", cmap="inferno")
plt.colorbar(im, ax=ax, label="12CO/13CO")

p = d.plotter()
p.plot_contour(ax, color="cyan", linewidth=0.7)

ax.set_title("Clump Identification")

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\ratio_clumps_map.fits",
    data.astype("float32"),
    header,
    overwrite=True)

plt.show()

#for leaf in d.leaves:
#    mask = leaf.get_mask()
#    values = data[mask]
#    
#    print("Clump ID:", leaf.idx)
#    print("Pixels:", np.sum(mask))
#    print("Mean ratio:", np.mean(values))
#    print("Max ratio:", np.max(values))
#    print()