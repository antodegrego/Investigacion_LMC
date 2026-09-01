import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astrodendro import Dendrogram
from astropy.stats import mad_std
from astrodendro.analysis import ppv_catalog


filename = 'NMolR_aca_13CO21_K_cube.fits'

with fits.open(filename) as hdul:
    data = hdul[0].data
    header = hdul[0].header

# --------------------------------------------------
# Beam y escala espacial
# --------------------------------------------------

bmaj = header['BMAJ'] * u.deg
bmin = header['BMIN'] * u.deg

pix_x = abs(header['CDELT1']) * u.deg
pix_y = abs(header['CDELT2']) * u.deg

beam_area = 1.1331 * bmaj * bmin
pixel_area = pix_x * pix_y

npix_beam = (beam_area / pixel_area).decompose().value

min_npix = int(np.ceil(2 * npix_beam))

print('BMAJ =', bmaj.to(u.arcsec))
print('BMIN =', bmin.to(u.arcsec))
print('Pixels por beam =', npix_beam)
print('2 beams =', min_npix, 'pixels')

noise_cube = data[:, 160:180, 90:120]

rms = np.nanstd(noise_cube)

print('RMS =', rms, 'K')


rms = mad_std(noise_cube, ignore_nan=True)

print('RMS =', rms, 'K')

min_value = 3.0 * rms
min_delta = 2.5 * rms

d = Dendrogram.compute(
    data,
    min_value=3.0 * rms,
    min_delta=2.5 * rms,
    min_npix=min_npix
)

print('Número total de estructuras:', len(d))
print('Número de trunks:', len(d.trunk))
print('Número de leaves:', len(d.leaves))

def structure_type(s):

    if s.is_leaf:
        return 'leaf'

    elif s.parent is None:
        return 'trunk'

    else:
        return 'branch'

print("NAXIS =", header['NAXIS'])

for i in range(1, header['NAXIS'] + 1):
    print(
        i,
        header.get(f'NAXIS{i}'),
        header.get(f'CTYPE{i}'),
        header.get(f'CUNIT{i}')
    )

metadata = {}

metadata['data_unit'] = u.K

metadata['spatial_scale'] = pix_x.to(u.arcsec)

metadata['beam_major'] = bmaj.to(u.arcsec)
metadata['beam_minor'] = bmin.to(u.arcsec)

frec13CO = 220.3986842
frec12CO = 230.538

metadata['wavelength'] = (
    frec13CO * u.GHz
).to(u.m, equivalencies=u.spectral())

# WCS original
wcs = WCS(header)

# Eliminar eje STOKES espurio
wcs = wcs.dropaxis(3)

metadata['wcs'] = wcs

print(metadata['wcs'])

metadata['vaxis'] = 0

dv = abs(header['CDELT3']) * u.Unit(header['CUNIT3'])

dv = dv.to(u.km / u.s)

metadata['velocity_scale'] = dv

print(header['CTYPE3'])
print(header['CUNIT3'])
print(header['CDELT3'])

cat = ppv_catalog(d, metadata)

print(cat.colnames)

cat[
    '_idx',
    'x_cen',
    'y_cen',
    'v_cen',
    'flux',
    'v_rms',
    'major_sigma',
    'minor_sigma',
    'position_angle'
].pprint(max_lines=30)

# Guardar dendrograma
d.save_to('dendrogram_13CO.fits')

# Guardar catálogo
cat.write(
    'catalog_13CO.fits',
    overwrite=True
)

print("Dendrograma guardado")
print("Catálogo guardado")