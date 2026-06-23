from astropy.io import fits

archivo = r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\12CO21_mask13CO_corrected.props.fits"
hdul = fits.open(archivo)
hdul.info()

from astropy.table import Table

t = Table.read(archivo)

t.write("NMoIR_12CO21_mask13CO_corrected.props.csv", format="csv", overwrite=True)