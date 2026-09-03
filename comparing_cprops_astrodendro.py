import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astrodendro import Dendrogram


# ============================================================
# PARÁMETROS
# ============================================================

# Distancia a la LMC
D_pc = 50000.0


# ============================================================
# LEER CATÁLOGOS
# ============================================================

# Catálogo CPROPS
cprops = pd.read_csv(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\NMoIR_12CO21_mask13CO_corrected.props.csv"
)

# Catálogo Astrodendro
astro = pd.read_csv(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\catalog_12CO_resolved.csv"
)


print("Columnas CPROPS:")
print(cprops.columns.tolist())

print("\nColumnas Astrodendro:")
print(astro.columns.tolist())


print("\nNúmero de estructuras CPROPS:", len(cprops))
print("Número de estructuras Astrodendro:", len(astro))


# ============================================================
# CARGAR DENDROGRAMA
# ============================================================

d = Dendrogram.load_from(
    r"C:\Users\adegr\OneDrive\Escritorio\Bibliografía Magister\Investigacion_LMC\dendrogram_12CO.fits"
)


# ============================================================
# IDENTIFICAR TIPO DE ESTRUCTURA ASTRODENDRO
# ============================================================

structure_types = []

for idx in astro['_idx']:

    s = d[int(idx)]

    if s.is_leaf:
        structure_types.append('leaf')

    elif s.parent is None:
        structure_types.append('trunk')

    else:
        structure_types.append('branch')


astro['structure_type'] = structure_types


print("\nEstructuras Astrodendro por tipo:")
print(astro['structure_type'].value_counts())


# ============================================================
# CONVERTIR RADIO ASTRODENDRO
# arcsec -> pc
# ============================================================

# radius está en arcsec.
#
# R_pc = theta_rad * D_pc
#
# theta_rad = radius_arcsec / 206265

astro['radius_pc'] = (
    astro['radius'] / 206265.0
) * D_pc


# ============================================================
# LIMPIAR VALORES NO FINITOS
# ============================================================

# CPROPS
cprops = cprops[
    np.isfinite(cprops['RAD_PC']) &
    np.isfinite(cprops['SIGV_KMS'])
].copy()


# Astrodendro
astro = astro[
    np.isfinite(astro['radius_pc']) &
    np.isfinite(astro['v_rms'])
].copy()


print("\nDespués de eliminar valores NaN:")
print("CPROPS:", len(cprops))
print("Astrodendro:", len(astro))


# ============================================================
# SEPARAR ASTRODENDRO POR TIPO
# ============================================================

leaves = astro[
    astro['structure_type'] == 'leaf'
]

branches = astro[
    astro['structure_type'] == 'branch'
]

trunks = astro[
    astro['structure_type'] == 'trunk'
]


print("\nNúmero de leaves:", len(leaves))
print("Número de branches:", len(branches))
print("Número de trunks:", len(trunks))


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
    cprops['RAD_PC'],
    cprops['SIGV_KMS'],
    marker='o',
    s=45,
    alpha=0.6,
    color='tab:blue',
    label='CPROPS'
)


# ============================================================
# ASTRODENDRO - LEAVES
# ============================================================

ax.scatter(
    leaves['radius_pc'],
    leaves['v_rms'],
    marker='o',
    s=45,
    alpha=0.6,
    color='tab:orange',
    label='Astrodendro - Leaves'
)


# ============================================================
# ASTRODENDRO - BRANCHES
# ============================================================

ax.scatter(
    branches['radius_pc'],
    branches['v_rms'],
    marker='^',
    s=55,
    alpha=0.6,
    color='tab:orange',
    label='Astrodendro - Branches'
)


# ============================================================
# ASTRODENDRO - TRUNKS
# ============================================================

ax.scatter(
    trunks['radius_pc'],
    trunks['v_rms'],
    marker='s',
    s=55,
    alpha=0.6,
    color='tab:orange',
    label='Astrodendro - Trunks'
)


# ============================================================
# FORMATO DEL GRÁFICO
# ============================================================

ax.set_xlabel(
    r'Radius [pc]'
)

ax.set_ylabel(
    r'$\sigma_v$ [km s$^{-1}$]'
)

ax.legend()

ax.grid(
    alpha=0.2
)

plt.tight_layout()

plt.show()
