from astropy.io import fits
import numpy as np
from scipy.ndimage import zoom

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
mask = (data13 > 3*rms_mom0) & np.isfinite(data13) & np.isfinite(data12) & (data12 > 3*rms_mom0)  & (data12 > 0)
ratio = np.where(mask, data13/data12, valor_negativo)
ratio[data12==0] = valor_negativo #np.nan
# Eliminar ratios mayores que 7 no fisicos causados por errores
ratio[ratio > 7] = valor_negativo

# Rebinning para evitar oversalpling
factor = 2/7
#factor = 1

data13_rebin = zoom(data13, factor, order=1)
data12_rebin = zoom(data12, factor, order=1)

# La misma región de ruido, escalada
y0, y1 = int(846*factor), int(902*factor)
x0, x1 = int(317*factor), int(440*factor)

noise_region = data13_rebin[y0:y1, x0:x1]

rms_rebin = np.nanstd(noise_region)

print(rms_rebin)

mask = (data13_rebin > 3*rms_rebin) & np.isfinite(data13_rebin) & np.isfinite(data12_rebin) & (data12_rebin > 3*rms_rebin)
ratio = np.where(mask, data13_rebin/data12_rebin, valor_negativo)
ratio[data12_rebin==0] = valor_negativo #np.nan
# Eliminar ratios mayores que 7 no fisicos causados por errores
ratio[ratio > 7] = valor_negativo

# asegurar que sea 2D
ratio = np.squeeze(ratio)

# corregir header
header["NAXIS"] = 2
header["NAXIS1"] = ratio.shape[1]
header["NAXIS2"] = ratio.shape[0]
header["CDELT1"] /= factor
header["CDELT2"] /= factor

# eliminar keywords de ejes extra si existen
for key in ["NAXIS3","CTYPE3","CRVAL3","CRPIX3","CDELT3","CUNIT3"]:
    if key in header:
        del header[key]


print(ratio.shape)
print(hdu12[0].header["NAXIS"])
print(np.nanmin(ratio))
print(np.nanmax(ratio))

fits.writeto(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_ratio_13CO21_12CO21_13mask_only_test2.fits",
    ratio,
    header,
    overwrite=True)


import matplotlib.pyplot as plt

plt.imshow(ratio, origin="lower", cmap="inferno", vmin=0, vmax=1)
plt.colorbar(label="13CO(2-1) / 12CO(2-1)")
plt.title("Line Ratio Map (13CO(2-1)/12CO(2-1))")
plt.show()

#################################
#Graficar cociente entre 13C
# O y 12CO
###################################

import matplotlib.pyplot as plt
import numpy as np
data12 = hdu12[0].data
data13 = hdu13[0].data

# Crear coordenadas de cada píxel
y, x = np.indices(data12.shape)

# asegurar misma forma
data12 = np.squeeze(data12)
data13 = np.squeeze(data13)

# ==========================================
# Definir regiones de las nubes
# ==========================================

y_orig, x_orig = np.indices(data12.shape)
cloud1_orig = (
    (x_orig >= 360) & (x_orig <= 860) &
    (y_orig >= 750) & (y_orig <= 1300)
)

cloud2_orig = (
    (x_orig >= 250) & (x_orig <= 700) &
    (y_orig >= 444) & (y_orig <= 740)
)

cloud3_orig = (
    (x_orig >= 90) & (x_orig <= 680) &
    (y_orig >= 20) & (y_orig <= 410)
)

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

data12 = zoom(data12, factor, order=1)
data13 = zoom(data13, factor, order=1)


# máscara
#mask = (data12 > 0) & (data13 > 0) & np.isfinite(data12) & np.isfinite(data13)
#mask = (data13 > 0) & np.isfinite(data13)
mask = (data13 > 3*rms_mom0) & np.isfinite(data13) & np.isfinite(data12) & (data12 > 3*rms_mom0)

mask = (
    (data13 > 3*rms_rebin) &
    np.isfinite(data12) &
    np.isfinite(data13) & 
    (data12 > 3*rms_mom0)
)


# extraer valores válidos
co12 = data12[mask]
co13 = data13[mask]

# scatter plot
plt.figure(figsize=(6,6))

plt.scatter(co12, co13, s=1, alpha=0.3)

plt.xlabel("Intensidad 12CO(2-1) [K km s-1]")
plt.ylabel("Intensidad 13CO(2-1) [K km s-1]")
plt.title("Intensidad 13CO(2-1) vs  Intensidad 12CO(2-1)")

# límites de ejes
plt.xlim(-0.5, np.max(co12))
plt.ylim(-0.5, np.max(co13))

# escala log
#plt.xscale("log")
#plt.yscale("log")

plt.grid(True, which="both", ls="--", alpha=0.5)
plt.show()

#Pesos para ajuste ponderado para resultado lineal
weights = co12**2

# Fiteo a los datos del scatter plot
coef = np.polyfit(co12, co13, 1)

m, b = coef

xfit = np.linspace(np.min(co12), np.max(co12), 1000)
yfit = np.polyval(coef, xfit)

plt.figure(figsize=(6,6))
plt.scatter(co12, co13, s=1, alpha=0.3)

plt.plot(
    xfit,
    yfit,
    'r',
    lw=2,
    label = fr"$y={m:.5f}x+{b:.5f}$" #(fr"$y={a:.5f}x^2+{b:.5f}x+{c:.3f}$")#fr"$y={a:.3e}x^2+{b:.3e}x+{c:.3e}$"
)

plt.legend()
plt.xlabel("Intensidad 12CO(2-1) [K km s$^{-1}$]")
plt.ylabel("Intensidad 13CO(2-1) [K km s$^{-1}$]")
plt.grid(alpha=0.5)
plt.show()

#################################
# Graficar 13CO vs 12CO pero distinguiendo las nubes
#################################

import matplotlib.pyplot as plt
import numpy as np


# ==========================================
# Coordenadas espaciales de los píxeles
# ==========================================

y, x = np.indices(data12.shape)

# ==========================================
# Máscara base de emisión
# ==========================================

mask_base = (
    (data13 > 3*rms_rebin) &
    (data12 > 3*rms_rebin) &
    np.isfinite(data12) &
    np.isfinite(data13)
)

# ==========================================
# Combinar región espacial + emisión válida
# ==========================================

mask_cloud1 = mask_base & cloud1_rebin
mask_cloud2 = mask_base & cloud2_rebin
mask_cloud3 = mask_base & cloud3_rebin

# ==========================================
# Extraer valores para cada nube
# ==========================================

co12_cloud1 = data12[mask_cloud1]
co13_cloud1 = data13[mask_cloud1]

co12_cloud2 = data12[mask_cloud2]
co13_cloud2 = data13[mask_cloud2]

co12_cloud3 = data12[mask_cloud3]
co13_cloud3 = data13[mask_cloud3]

# ==========================================
# Scatter plot
# ==========================================

plt.figure(figsize=(6,6))


plt.scatter(
    co12_cloud3,
    co13_cloud3,
    s=2,
    alpha=0.3,
    color='green',
    label='N159'
)

plt.scatter(
    co12_cloud2,
    co13_cloud2,
    s=2,
    alpha=0.3,
    color='blue',
    label='N160'
)

plt.scatter(
    co12_cloud1,
    co13_cloud1,
    s=2,
    alpha=0.3,
    color='red',
    label='N158'
)

plt.xlabel("Intensidad 12CO(2-1) [K km s$^{-1}$]")
plt.ylabel("Intensidad 13CO(2-1) [K km s$^{-1}$]")

plt.title(
    "Intensidad 13CO(2-1) vs 12CO(2-1)"
)

plt.legend()

plt.grid(
    True,
    which="both",
    ls="--",
    alpha=0.5
)

plt.show()

#############################################
# BOOTSTRAP LINEAR FIT
# Ajuste global + ajuste individual por nube
#############################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


# ============================================================
# FUNCIÓN PARA REALIZAR BOOTSTRAP
# ============================================================

def bootstrap_linear_fit(x, y, n_boot=10000, seed=42):

    # Asegurar arrays numpy
    x = np.asarray(x)
    y = np.asarray(y)

    # Eliminar valores no válidos
    valid = (
        np.isfinite(x) &
        np.isfinite(y)
    )

    x = x[valid]
    y = y[valid]

    # Número de puntos
    N = len(x)

    print("Número de puntos:", N)

    # Generador aleatorio
    rng = np.random.default_rng(seed)

    # Arrays para guardar resultados
    slopes = np.zeros(n_boot)
    intercepts = np.zeros(n_boot)

    # ========================================================
    # Bootstrap
    # ========================================================

    for i in range(n_boot):

        # Seleccionar índices con reemplazo
        indices = rng.integers(
            0,
            N,
            size=N
        )

        # Muestra bootstrap
        x_boot = x[indices]
        y_boot = y[indices]

        # Ajuste lineal
        fit = linregress(
            x_boot,
            y_boot
        )

        # Guardar resultados
        slopes[i] = fit.slope
        intercepts[i] = fit.intercept

    # ========================================================
    # Valores centrales
    # ========================================================

    # Usamos la mediana de las distribuciones bootstrap
    m = np.median(slopes)
    b = np.median(intercepts)

    # ========================================================
    # Intervalos de confianza de 68%
    # ========================================================

    m_low, m_high = np.percentile(
        slopes,
        [16, 84]
    )

    b_low, b_high = np.percentile(
        intercepts,
        [16, 84]
    )

    # Errores simétricos aproximados
    m_err = np.std(slopes)
    b_err = np.std(intercepts)

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


# ============================================================
# AJUSTE GLOBAL
# ============================================================

# Combinar todas las nubes

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


# Bootstrap global

fit_global = bootstrap_linear_fit(
    co12_all,
    co13_all,
    n_boot=10000,
    seed=42
)


# ============================================================
# AJUSTE NUBE 1
# ============================================================

fit_cloud1 = bootstrap_linear_fit(
    co12_cloud1,
    co13_cloud1,
    n_boot=10000,
    seed=43
)


# ============================================================
# AJUSTE NUBE 2
# ============================================================

fit_cloud2 = bootstrap_linear_fit(
    co12_cloud2,
    co13_cloud2,
    n_boot=10000,
    seed=44
)


# ============================================================
# AJUSTE NUBE 3
# ============================================================

fit_cloud3 = bootstrap_linear_fit(
    co12_cloud3,
    co13_cloud3,
    n_boot=10000,
    seed=45
)


# ============================================================
# IMPRIMIR RESULTADOS
# ============================================================

print("\n==========================================")
print("RESULTADOS DEL BOOTSTRAP")
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


print("\nNUBE 1")
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


print("\nNUBE 2")
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


print("\nNUBE 3")
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


# ============================================================
# PREPARAR RECTAS
# ============================================================

# Rango de X común

x_min = np.min(co12_all)
x_max = np.max(co12_all)

xfit = np.linspace(
    x_min,
    x_max,
    1000
)


# Ajuste global

yfit_global = (
    fit_global["m"] * xfit
    + fit_global["b"]
)


# Ajuste nube 1

yfit_cloud1 = (
    fit_cloud1["m"] * xfit
    + fit_cloud1["b"]
)


# Ajuste nube 2

yfit_cloud2 = (
    fit_cloud2["m"] * xfit
    + fit_cloud2["b"]
)


# Ajuste nube 3

yfit_cloud3 = (
    fit_cloud3["m"] * xfit
    + fit_cloud3["b"]
)


# ============================================================
# GRAFICAR
# ============================================================

plt.figure(figsize=(6, 6))


# ------------------------------------------------------------
# Puntos Nube 3
# ------------------------------------------------------------

plt.scatter(
    co12_cloud3,
    co13_cloud3,
    s=3,
    alpha=0.15,
    color="green",
    label="N159"
)


# ------------------------------------------------------------
# Puntos Nube 2
# ------------------------------------------------------------

plt.scatter(
    co12_cloud2,
    co13_cloud2,
    s=3,
    alpha=0.20,
    color="blue",
    label="N160"
)


# ------------------------------------------------------------
# Puntos Nube 1
# ------------------------------------------------------------

plt.scatter(
    co12_cloud1,
    co13_cloud1,
    s=3,
    alpha=0.30,
    color="red",
    label="N158"
)


# ============================================================
# AJUSTE GLOBAL
# ============================================================

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


# ============================================================
# AJUSTE NUBE 1
# ============================================================

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


# ============================================================
# AJUSTE NUBE 2
# ============================================================

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


# ============================================================
# AJUSTE NUBE 3
# ============================================================

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


# ============================================================
# ETIQUETAS
# ============================================================

plt.xlabel(
    "Intensidad 12CO(2-1) [K km s$^{-1}$]"
)

plt.ylabel(
    "Intensidad 13CO(2-1) [K km s$^{-1}$]"
)

plt.title(
    "Intensidad 13CO(2-1) vs 12CO(2-1)"
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