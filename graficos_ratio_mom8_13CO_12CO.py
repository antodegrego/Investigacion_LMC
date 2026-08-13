from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
from scipy.stats import linregress


# ==========================================================
# Leer mapas momento 8
# ==========================================================

hdu12 = fits.open(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment8_12CO21.fits"
)

hdu13 = fits.open(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment8_13CO21.fits"
)

hdu12_m0 = fits.open(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_12CO21.fits"
)

hdu13_m0 = fits.open(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\moment0_13CO21.fits"
)

mom0_12 = np.squeeze(hdu12_m0[0].data)
mom0_13 = np.squeeze(hdu13_m0[0].data)

# Guardar datos originales
data12_orig = np.squeeze(hdu12[0].data)
data13_orig = np.squeeze(hdu13[0].data)

header = hdu12[0].header.copy()

valor_negativo = -1e6


# ==========================================================
# Factor de rebinning
# ==========================================================

factor = 2/7


# ==========================================================
# Definir regiones de las nubes
# EN COORDENADAS ORIGINALES
# ==========================================================

y_orig, x_orig = np.indices(data12_orig.shape)


# N159
cloud3_orig = (
    (x_orig >= 90) & (x_orig <= 680) &
    (y_orig >= 20) & (y_orig <= 410)
)


# N160
cloud2_orig = (
    (x_orig >= 250) & (x_orig <= 700) &
    (y_orig >= 444) & (y_orig <= 740)
)


# N158
cloud1_orig = (
    (x_orig >= 360) & (x_orig <= 860) &
    (y_orig >= 750) & (y_orig <= 1300)
)


# ==========================================================
# Rebinning de los mapas
# ==========================================================

data12 = zoom(
    data12_orig,
    factor,
    order=1
)

data13 = zoom(
    data13_orig,
    factor,
    order=1
)

mom0_12 = zoom(mom0_12, factor, order=1)
mom0_13 = zoom(mom0_13, factor, order=1)

# ==========================================================
# Rebinning de las máscaras de las nubes
# order=0 para no interpolar las máscaras
# ==========================================================

cloud1_rebin = zoom(
    cloud1_orig.astype(float),
    factor,
    order=0
).astype(bool)


cloud2_rebin = zoom(
    cloud2_orig.astype(float),
    factor,
    order=0
).astype(bool)


cloud3_rebin = zoom(
    cloud3_orig.astype(float),
    factor,
    order=0
).astype(bool)


# ==========================================================
# Calcular RMS del mapa rebineado
# ==========================================================

y0, y1 = int(846*factor), int(902*factor)
x0, x1 = int(317*factor), int(440*factor)

noise_region = data13[y0:y1, x0:x1]

rms = np.nanstd(noise_region)

print(f"RMS = {rms:.5f}")

noise_region_12 = mom0_12[y0:y1, x0:x1]
noise_region_13 = mom0_13[y0:y1, x0:x1]

rms_mom0_12 = np.nanstd(noise_region_12)
rms_mom0_13 = np.nanstd(noise_region_13)


# ==========================================================
# Máscara base de emisión
# ==========================================================

mask_emission = (
    (mom0_12 > 3*rms_mom0_12) &
    (mom0_13 > 3*rms_mom0_13) &
    np.isfinite(mom0_12) &
    np.isfinite(mom0_13)
)

mask_base = (
    (data13 > 3*rms) &
    (data12 > 3*rms) &
    np.isfinite(data12) &
    np.isfinite(data13) &
    (data12 != 0) &
    mask_emission
)


# ==========================================================
# Ratio
# ==========================================================

ratio = np.full(
    data12.shape,
    valor_negativo,
    dtype=float
)

ratio[mask_base] = (
    data13[mask_base] /
    data12[mask_base]
)


# Eliminar ratios no físicos muy altos
ratio[ratio > 7] = valor_negativo


# ==========================================================
# Actualizar header
# ==========================================================

header["NAXIS"] = 2
header["NAXIS1"] = ratio.shape[1]
header["NAXIS2"] = ratio.shape[0]

header["CDELT1"] /= factor
header["CDELT2"] /= factor


for key in [
    "NAXIS3",
    "CTYPE3",
    "CRVAL3",
    "CRPIX3",
    "CDELT3",
    "CUNIT3"
]:

    if key in header:
        del header[key]


# ==========================================================
# Guardar FITS
# ==========================================================

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_ratio_M8_13CO21_12CO21_rebinned.fits",
    ratio,
    header,
    overwrite=True
)


# ==========================================================
# Mostrar mapa de ratio
# ==========================================================

plt.figure(figsize=(8,6))

plt.imshow(
    ratio,
    origin="lower",
    cmap="inferno",
    vmin=0,
    vmax=0.1
)

plt.colorbar(
    label="Peak Temperature: Ratio (13CO / 12CO)"
)

plt.title(
    "Peak Temperature Ratio Map"
)

plt.tight_layout()

plt.show()


# ==========================================================
# Crear máscaras individuales para cada nube
# ==========================================================

mask_cloud1 = (
    mask_base &
    cloud1_rebin
)

mask_cloud2 = (
    mask_base &
    cloud2_rebin
)

mask_cloud3 = (
    mask_base &
    cloud3_rebin
)


# ==========================================================
# Extraer datos de cada nube
# ==========================================================

# N158

co12_cloud1 = data12[mask_cloud1]
co13_cloud1 = data13[mask_cloud1]


# N160

co12_cloud2 = data12[mask_cloud2]
co13_cloud2 = data13[mask_cloud2]


# N159

co12_cloud3 = data12[mask_cloud3]
co13_cloud3 = data13[mask_cloud3]


# ==========================================================
# Combinar todas las nubes
# ==========================================================

co12_all = np.concatenate([
    co12_cloud1,
    co12_cloud2,
    co12_cloud3
])

co13_all = np.concatenate([
    co13_cloud1,
    co13_cloud2,
    co13_cloud3
])


# ==========================================================
# Scatter plot separado por nubes
# ==========================================================

plt.figure(figsize=(6,6))


# N159

plt.scatter(
    co12_cloud3,
    co13_cloud3,
    s=3,
    alpha=0.15,
    color="green",
    label="N159"
)


# N160

plt.scatter(
    co12_cloud2,
    co13_cloud2,
    s=3,
    alpha=0.20,
    color="blue",
    label="N160"
)


# N158

plt.scatter(
    co12_cloud1,
    co13_cloud1,
    s=3,
    alpha=0.30,
    color="red",
    label="N158"
)


plt.xlabel(
    "Peak Temperature 12CO(2-1) [K]"
)

plt.ylabel(
    "Peak Temperature 13CO(2-1) [K]"
)

plt.title(
    "Peak temperature: 13CO(2-1) vs 12CO(2-1)"
)

plt.legend()

plt.grid(
    alpha=0.4
)

plt.tight_layout()

plt.show()


# ==========================================================
# FUNCIÓN BOOTSTRAP
# ==========================================================

def bootstrap_linear_fit(
    x,
    y,
    n_boot=10000,
    seed=42
):

    # Convertir a numpy
    x = np.asarray(x)
    y = np.asarray(y)


    # Eliminar valores inválidos
    valid = (
        np.isfinite(x) &
        np.isfinite(y)
    )

    x = x[valid]
    y = y[valid]


    # Número de puntos
    N = len(x)

    print(
        f"Número de puntos = {N}"
    )


    # Generador aleatorio
    rng = np.random.default_rng(seed)


    # Arrays para guardar resultados
    slopes = np.zeros(n_boot)

    intercepts = np.zeros(n_boot)


    # ======================================================
    # Bootstrap
    # ======================================================

    for i in range(n_boot):

        # Muestreo con reemplazo
        indices = rng.integers(
            0,
            N,
            size=N
        )


        # Datos bootstrap
        x_boot = x[indices]
        y_boot = y[indices]


        # Ajuste lineal
        fit = linregress(
            x_boot,
            y_boot
        )


        # Guardar pendiente
        slopes[i] = fit.slope


        # Guardar intercepto
        intercepts[i] = fit.intercept


    # ======================================================
    # Valores centrales
    # ======================================================

    m = np.median(slopes)

    b = np.median(intercepts)


    # ======================================================
    # Errores
    # ======================================================

    m_err = np.std(slopes)

    b_err = np.std(intercepts)


    # Intervalos de confianza 68%
    m_low, m_high = np.percentile(
        slopes,
        [16,84]
    )

    b_low, b_high = np.percentile(
        intercepts,
        [16,84]
    )


    return {

        "m": m,

        "b": b,

        "m_err": m_err,

        "b_err": b_err,

        "m_low": m_low,

        "m_high": m_high,

        "b_low": b_low,

        "b_high": b_high,

        "slopes": slopes,

        "intercepts": intercepts,

        "x": x,

        "y": y

    }


# ==========================================================
# BOOTSTRAP GLOBAL
# ==========================================================

fit_global = bootstrap_linear_fit(
    co12_all,
    co13_all,
    n_boot=10000,
    seed=42
)


# ==========================================================
# BOOTSTRAP N158
# ==========================================================

fit_cloud1 = bootstrap_linear_fit(
    co12_cloud1,
    co13_cloud1,
    n_boot=10000,
    seed=43
)


# ==========================================================
# BOOTSTRAP N160
# ==========================================================

fit_cloud2 = bootstrap_linear_fit(
    co12_cloud2,
    co13_cloud2,
    n_boot=10000,
    seed=44
)


# ==========================================================
# BOOTSTRAP N159
# ==========================================================

fit_cloud3 = bootstrap_linear_fit(
    co12_cloud3,
    co13_cloud3,
    n_boot=10000,
    seed=45
)


# ==========================================================
# IMPRIMIR RESULTADOS
# ==========================================================

print("\n==========================================")
print("RESULTADOS DEL BOOTSTRAP - MOMENTO 8")
print("==========================================")


print("\nAJUSTE GLOBAL")
print("------------------------------------------")

print(
    f"Pendiente = "
    f"{fit_global['m']:.5f} ± "
    f"{fit_global['m_err']:.5f}"
)

print(
    f"Intercepto = "
    f"{fit_global['b']:.5f} ± "
    f"{fit_global['b_err']:.5f}"
)


print("\nN158")
print("------------------------------------------")

print(
    f"Pendiente = "
    f"{fit_cloud1['m']:.5f} ± "
    f"{fit_cloud1['m_err']:.5f}"
)

print(
    f"Intercepto = "
    f"{fit_cloud1['b']:.5f} ± "
    f"{fit_cloud1['b_err']:.5f}"
)


print("\nN160")
print("------------------------------------------")

print(
    f"Pendiente = "
    f"{fit_cloud2['m']:.5f} ± "
    f"{fit_cloud2['m_err']:.5f}"
)

print(
    f"Intercepto = "
    f"{fit_cloud2['b']:.5f} ± "
    f"{fit_cloud2['b_err']:.5f}"
)


print("\nN159")
print("------------------------------------------")

print(
    f"Pendiente = "
    f"{fit_cloud3['m']:.5f} ± "
    f"{fit_cloud3['m_err']:.5f}"
)

print(
    f"Intercepto = "
    f"{fit_cloud3['b']:.5f} ± "
    f"{fit_cloud3['b_err']:.5f}"
)


# ==========================================================
# Preparar rectas
# ==========================================================

x_min = np.min(co12_all)

x_max = np.max(co12_all)

xfit = np.linspace(
    x_min,
    x_max,
    1000
)


# Recta global

yfit_global = (
    fit_global["m"] * xfit
    + fit_global["b"]
)


# Recta N158

yfit_cloud1 = (
    fit_cloud1["m"] * xfit
    + fit_cloud1["b"]
)


# Recta N160

yfit_cloud2 = (
    fit_cloud2["m"] * xfit
    + fit_cloud2["b"]
)


# Recta N159

yfit_cloud3 = (
    fit_cloud3["m"] * xfit
    + fit_cloud3["b"]
)


# ==========================================================
# GRÁFICO FINAL
# ==========================================================

plt.figure(figsize=(7,7))


# ==========================================================
# PUNTOS N159
# ==========================================================

plt.scatter(
    co12_cloud3,
    co13_cloud3,
    s=3,
    alpha=0.15,
    color="green",
    label="N159"
)


# ==========================================================
# PUNTOS N160
# ==========================================================

plt.scatter(
    co12_cloud2,
    co13_cloud2,
    s=3,
    alpha=0.20,
    color="blue",
    label="N160"
)


# ==========================================================
# PUNTOS N158
# ==========================================================

plt.scatter(
    co12_cloud1,
    co13_cloud1,
    s=3,
    alpha=0.30,
    color="red",
    label="N158"
)


# ==========================================================
# AJUSTE GLOBAL
# ==========================================================

plt.plot(
    xfit,
    yfit_global,
    color="black",
    lw=2.5,
    label=(
        "Global: "
        fr"$y=({fit_global['m']:.2g}"
        fr"\pm{fit_global['m_err']:.2g})x"
        fr"+({fit_global['b']:.2g}"
        fr"\pm{fit_global['b_err']:.2g})$"
    )
)


# ==========================================================
# AJUSTE N158
# ==========================================================

plt.plot(
    xfit,
    yfit_cloud1,
    color="red",
    lw=2,
    linestyle="--",
    label=(
        "N158: "
        fr"$y=({fit_cloud1['m']:.2g}"
        fr"\pm{fit_cloud1['m_err']:.2g})x"
        fr"+({fit_cloud1['b']:.2g}"
        fr"\pm{fit_cloud1['b_err']:.2g})$"
    )
)


# ==========================================================
# AJUSTE N160
# ==========================================================

plt.plot(
    xfit,
    yfit_cloud2,
    color="blue",
    lw=2,
    linestyle="--",
    label=(
        "N160: "
        fr"$y=({fit_cloud2['m']:.2g}"
        fr"\pm{fit_cloud2['m_err']:.2g})x"
        fr"+({fit_cloud2['b']:.2g}"
        fr"\pm{fit_cloud2['b_err']:.2g})$"
    )
)


# ==========================================================
# AJUSTE N159
# ==========================================================

plt.plot(
    xfit,
    yfit_cloud3,
    color="green",
    lw=2,
    linestyle="--",
    label=(
        "N159: "
        fr"$y=({fit_cloud3['m']:.2g}"
        fr"\pm{fit_cloud3['m_err']:.2g})x"
        fr"+({fit_cloud3['b']:.2g}"
        fr"\pm{fit_cloud3['b_err']:.2g})$"
    )
)


# ==========================================================
# ETIQUETAS
# ==========================================================

plt.xlabel(
    "Peak Temperature 12CO(2-1) [K]"
)

plt.ylabel(
    "Peak Temperature 13CO(2-1) [K]"
)

plt.title(
    "T_peak[13CO(2-1)] vs T_peak[12CO(2-1)]"
)

plt.legend(
    fontsize=8
)

plt.grid(
    True,
    alpha=0.4
)

plt.tight_layout()

plt.show()

print("================================")
print("DIAGNÓSTICO")
print("================================")

print("RMS =", rms)

print("Píxeles totales =", data12.size)

print("Píxeles que cumplen máscara =", np.sum(mask_base))

print("Fracción válida =", np.sum(mask_base)/data12.size)

print("12CO mínimo válido =", np.min(data12[mask_base]))
print("12CO máximo válido =", np.max(data12[mask_base]))

print("13CO mínimo válido =", np.min(data13[mask_base]))
print("13CO máximo válido =", np.max(data13[mask_base]))

print("Ratio mínimo =", np.min(ratio[mask_base]))
print("Ratio máximo =", np.max(ratio[mask_base]))

print("Ratio mediana =", np.median(ratio[mask_base]))