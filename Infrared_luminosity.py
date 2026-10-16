import numpy as np

# ============================================================
# CONSTANTES
# ============================================================

# Masa del átomo de hidrógeno
m_H = 1.6735575e-24       # g

# Masa solar
M_sun = 1.98847e33        # g

# Conversión erg/s -> luminosidad solar
L_sun = 3.828e33          # erg/s


# ============================================================
# PARÁMETROS DEL MODELO DL07
# ============================================================

# Placeholder: mezcla de polvo SMC
# Dust-to-hydrogen mass ratio
phi_d_H = 0.00206

# Placeholder:
# potencia emitida por H para U = 1
#
# IMPORTANTE:
# Este valor debe ser reemplazado por el valor exacto
# que quieras adoptar de DL07 / Chastenet et al.
power_per_H = 4.0e-24     # erg s^-1 H^-1  <-- PLACEHOLDER

# Máxima intensidad de radiación utilizada para la
# componente tipo PDR
U_max = 1e3


# ============================================================
# CALCULAR P0
# ============================================================

# P0 = (power/H) / (phi_d/H * m_H)
#
# unidades:
# power/H       -> erg s^-1 H^-1
# phi_d_H*m_H   -> g/H
#
# resultado:
# P0 -> erg s^-1 g^-1

P0 = power_per_H / (phi_d_H * m_H)

print("P0 =", P0, "erg s^-1 g^-1")


# ============================================================
# FUNCION PARA CALCULAR L_TIR
# ============================================================

def calculate_L_TIR(U, M_dust_Msun):
    """
    Calcula la luminosidad infrarroja total (8-1000 um)
    utilizando:

        L_TIR = P0 * U * M_dust

    Parameters
    ----------
    U : float or array
        Intensidad media del campo de radiación.

    M_dust_Msun : float or array
        Masa de polvo en masas solares.

    Returns
    -------
    L_TIR_erg_s : float or array
        Luminosidad TIR en erg/s.

    L_TIR_Lsun : float or array
        Luminosidad TIR en L_sun.
    """

    # Convertir masa de polvo de Msun a gramos
    M_dust_g = M_dust_Msun * M_sun

    # Luminosidad TIR
    L_TIR_erg_s = P0 * U * M_dust_g

    # Convertir a luminosidades solares
    L_TIR_Lsun = L_TIR_erg_s / L_sun

    return L_TIR_erg_s, L_TIR_Lsun


# ============================================================
# EJEMPLO
# ============================================================

U = 5.0
M_dust = 1e4       # Msun

L_TIR_erg_s, L_TIR_Lsun = calculate_L_TIR(U, M_dust)

print()
print("Ejemplo:")
print("U =", U)
print("M_dust =", M_dust, "Msun")
print("L_TIR =", L_TIR_erg_s, "erg/s")
print("L_TIR =", L_TIR_Lsun, "Lsun")