import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import linregress


# ============================================================
# PARÁMETROS
# ============================================================

# Distancia a la LMC
D_pc = 50000.0

# Resolución de la imagen
# BMAJ = BMIN = 7 arcsec
resolution_arcsec = 7.0


# ============================================================
# CONVERTIR RESOLUCIÓN A pc
# ============================================================

resolution_pc = (
    resolution_arcsec / 206265.0
) * D_pc

print("Resolución angular:", resolution_arcsec, "arcsec")
print("Resolución espacial:", resolution_pc, "pc")


# ============================================================
# LEER CATÁLOGOS
# ============================================================

# CPROPS
cprops = pd.read_csv(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_12CO21_mask13CO_corrected.props.csv"
)

# Astrodendro
astro = pd.read_csv(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\catalog_12CO_resolved.csv"
)


print("\nNúmero inicial de estructuras:")
print("CPROPS:", len(cprops))
print("Astrodendro:", len(astro))


# ============================================================
# CONVERTIR RADIO ASTRODENDRO
# arcsec -> pc
# ============================================================

astro['radius_pc'] = (
    astro['radius'] / 206265.0
) * D_pc


# ============================================================
# ELIMINAR VALORES NO FINITOS
# ============================================================

cprops = cprops[
    np.isfinite(cprops['RAD_PC']) &
    np.isfinite(cprops['SIGV_KMS'])
].copy()


astro = astro[
    np.isfinite(astro['radius_pc']) &
    np.isfinite(astro['v_rms'])
].copy()


# ============================================================
# ELIMINAR ESTRUCTURAS NO RESUELTAS
# ============================================================

# CPROPS

n_cprops_before = len(cprops)

cprops = cprops[
    cprops['RAD_PC'] >= resolution_pc
].copy()

n_cprops_removed = (
    n_cprops_before - len(cprops)
)


# Astrodendro

n_astro_before = len(astro)

astro = astro[
    astro['radius_pc'] >= resolution_pc
].copy()

n_astro_removed = (
    n_astro_before - len(astro)
)


print("\nEstructuras eliminadas por ser menores que la resolución:")

print(
    "CPROPS:",
    n_cprops_removed
)

print(
    "Astrodendro:",
    n_astro_removed
)


print("\nEstructuras restantes:")

print(
    "CPROPS:",
    len(cprops)
)

print(
    "Astrodendro:",
    len(astro)
)


# ============================================================
# FUNCIÓN PARA AJUSTE DE LEY DE POTENCIA
# ============================================================

def power_law_fit(x, y):

    # El ajuste se realiza en espacio logarítmico

    log_x = np.log10(x)
    log_y = np.log10(y)

    fit = linregress(
        log_x,
        log_y
    )

    alpha = fit.slope
    log_A = fit.intercept

    A = 10**log_A

    r_value = fit.rvalue

    r_squared = r_value**2

    return A, alpha, r_squared


# ============================================================
# DATOS PARA LOS FITEOS
# ============================================================

x_cprops = cprops['RAD_PC'].values
y_cprops = cprops['SIGV_KMS'].values

x_astro = astro['radius_pc'].values
y_astro = astro['v_rms'].values


# ============================================================
# FITEO CPROPS
# ============================================================

A_cprops, alpha_cprops, r2_cprops = power_law_fit(
    x_cprops,
    y_cprops
)


print("\n==============================")
print("FIT CPROPS")
print("==============================")

print(
    f"sigma_v = {A_cprops:.3f} * R^{alpha_cprops:.3f}"
)

print(
    f"alpha = {alpha_cprops:.3f}"
)

print(
    f"R² = {r2_cprops:.3f}"
)


# ============================================================
# FITEO ASTRODENDRO
# ============================================================

A_astro, alpha_astro, r2_astro = power_law_fit(
    x_astro,
    y_astro
)


print("\n==============================")
print("FIT ASTRODENDRO")
print("==============================")

print(
    f"sigma_v = {A_astro:.3f} * R^{alpha_astro:.3f}"
)

print(
    f"alpha = {alpha_astro:.3f}"
)

print(
    f"R² = {r2_astro:.3f}"
)


# ============================================================
# GENERAR CURVAS DE FIT
# ============================================================

# Rango de radios

R_min = min(
    x_cprops.min(),
    x_astro.min()
)

R_max = max(
    x_cprops.max(),
    x_astro.max()
)

R_fit = np.logspace(
    np.log10(R_min),
    np.log10(R_max),
    200
)


# CPROPS

sigma_fit_cprops = (
    A_cprops *
    R_fit**alpha_cprops
)


# Astrodendro

sigma_fit_astro = (
    A_astro *
    R_fit**alpha_astro
)


# ============================================================
# FIGURA
# ============================================================

fig, ax = plt.subplots(
    figsize=(8, 6)
)


# ============================================================
# CPROPS
# ============================================================

ax.scatter(
    x_cprops,
    y_cprops,
    marker='o',
    s=45,
    alpha=0.6,
    color='tab:blue',
    label='CPROPS'
)


# ============================================================
# ASTRODENDRO
# ============================================================

ax.scatter(
    x_astro,
    y_astro,
    marker='o',
    s=45,
    alpha=0.6,
    color='tab:orange',
    label='Astrodendro'
)


# ============================================================
# FIT CPROPS
# ============================================================

ax.plot(
    R_fit,
    sigma_fit_cprops,
    color='tab:blue',
    linewidth=2,
    label=(
        rf'CPROPS fit: '
        rf'$\sigma_v \propto R^{{{alpha_cprops:.2f}}}$'
    )
)


# ============================================================
# FIT ASTRODENDRO
# ============================================================

ax.plot(
    R_fit,
    sigma_fit_astro,
    color='tab:orange',
    linewidth=2,
    label=(
        rf'Astrodendro fit: '
        rf'$\sigma_v \propto R^{{{alpha_astro:.2f}}}$'
    )
)


# ============================================================
# ESCALA LOGARÍTMICA
# ============================================================

ax.set_xscale('log')
ax.set_yscale('log')


# ============================================================
# RESOLUCIÓN
# ============================================================

ax.axvline(
    resolution_pc,
    linestyle='--',
    color='gray',
    alpha=0.7,
    label=(
        rf'Resolution = {resolution_pc:.2f} pc'
    )
)


# ============================================================
# ETIQUETAS
# ============================================================

ax.set_xlabel(
    r'Radius [pc]'
)

ax.set_ylabel(
    r'$\sigma_v$ [km s$^{-1}$]'
)


ax.legend()

ax.grid(
    alpha=0.2,
    which='both'
)


plt.tight_layout()

plt.show()


