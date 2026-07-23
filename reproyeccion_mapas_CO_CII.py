from astropy.io import fits
from reproject import reproject_interp

# mapa CII (grilla objetivo)
hdu_cii = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\line_flux_CII.fits")[0]

# mapa CO ya convolucionado a 15.3"
hdu_co = fits.open(r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_moment0_12CO21_15p3arcsec_Wm2sr.fits")[0]

co_reproj, footprint = reproject_interp(
    hdu_co,
    hdu_cii.header
)

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_12CO_15p3arcsec_Wm2sr_reprojected.fits",
    co_reproj,
    hdu_cii.header,
    overwrite=True
)