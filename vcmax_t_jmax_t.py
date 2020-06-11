import numpy as np
import matplotlib.pyplot as plt

def xt3(k25, Tk, Ea, deltaS, Hd):
    RGAS = 8.31429958
    TREFK = 25.0 + 273.15 #298.200012

    coef = 1.0 + np.exp((deltaS * TREFK - Hd) / (RGAS * TREFK))
    num = coef * np.exp((Ea / (RGAS * TREFK)) * (1. - TREFK / Tk))
    den = 1.0 + np.exp((deltaS * Tk - Hd) / (RGAS * Tk))

    return k25 * np.maximum(0.0, num / den)

def xvcmxt3_acclim(k25, Tk, Ea, deltaS, Hd):
    # Acclimated temperature response
    # Kumarathunge et al., New. Phyt., 2019,
    # Eq 7

    RGAS = 8.31429958
    TREFK = 25.0 + 273.15 #298.200012

    coef = 1.0 + np.exp((deltaS * TREFK - Hd) / (RGAS * TREFK))
    num = coef * np.exp((Ea / (RGAS * TREFK)) * (1. - TREFK / Tk))
    den = 1.0 + np.exp((deltaS * Tk - Hd) / (RGAS * Tk))

    return k25 * np.maximum(0.0, num / den)

def arrh(k25, Ea, Tk):
    RGAS = 8.314
    return k25 * np.exp((Ea * (Tk - 298.15)) / (298.15 * RGAS * Tk))

def peaked_arrh(k25, Ea, Tk, deltaS, Hd):
    RGAS = 8.314
    arg1 = arrh(k25, Ea, Tk)
    arg2 = 1.0 + np.exp((298.15 * deltaS - Hd) / (298.15 * RGAS))
    arg3 = 1.0 + np.exp((Tk * deltaS - Hd) / (Tk * RGAS))

    return arg1 * arg2 / arg3

def peaked_arrh_acclim(k25, Ea, Tk, deltaS, Hd):
    RGAS = 8.314
    arg1 = arrh(k25, Ea, Tk)
    arg2 = 1.0 + np.exp((298.15 * deltaS - Hd) / (298.15 * RGAS))
    arg3 = 1.0 + np.exp((Tk * deltaS - Hd) / (Tk * RGAS))

    return arg1 * arg2 / arg3


kJ_to_J = 1000.0
Tleaf = np.arange(1.0, 50.0)
TleafK = Tleaf + 273.15

Vcmax25 = 80.52204642
Jmax25 = 138.

# CABLE EBF
Eaj = 50300.0
Eav = 73637.0
deltaSj = 495.0
deltaSv = 486.0
Hdv = 149252.0
Hdj = 152044.0

Vcmax_cable = xt3(Vcmax25, TleafK, Eav, deltaSv, Hdv)
Jmax_cable = xt3(Jmax25, TleafK, Eaj, deltaSj, Hdj)
#Vcmax_acclim = xvcmxt3_acclim(Vcmax25, TleafK, Eav_acclim, deltaSv_acclim,
#                              Hdv_acclim)


# Eucalyptus tereticornis
Eaj = 32292.6319602695
Eav = 66386.2007802811
deltaSj = 638.055941746436
deltaSv = 639.602873767181
Hdv = 200000.0
Hdj = 200000.0

Vcmax_euc = peaked_arrh(Vcmax25, Eav, TleafK, deltaSv, Hdv)
Jmax_euc = peaked_arrh(Jmax25, Eaj, TleafK, deltaSj, Hdj)

# Acclimation based on Kumarathunge
Tgrowth = 15.0 # deg C
Eav_acclim = (42.6 * kJ_to_J) + (1.14 * kJ_to_J) * Tgrowth  # J mol-1
Eaj_acclim = 40.71 * kJ_to_J # J mol-1
Hdv_acclim = 200000.0
Hdj_acclim = 200000.0
deltaSv_acclim = 645.13 - 0.38 * Tgrowth # J mol-1 K-1
deltaSj_acclim = 658.77 - 0.84 * Tgrowth - 0.52 # J mol-1 K-1

Vcmax_euc_acclim = peaked_arrh(Vcmax25, Eav_acclim, TleafK, deltaSv_acclim,
                              Hdv_acclim)
Jmax_euc_acclim = peaked_arrh(Jmax25, Eaj_acclim, TleafK, deltaSj_acclim,
                             Hdj_acclim)

fig = plt.figure(figsize=(14,6))
fig.subplots_adjust(hspace=0.1)
fig.subplots_adjust(wspace=0.2)
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Helvetica"
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.size'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14

ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.plot(Tleaf, Vcmax_cable, label="CABLE EBF")
ax1.plot(Tleaf, Vcmax_euc, label="Eucalyptus tereticornis")
ax1.plot(Tleaf, Vcmax_euc_acclim, label="Eucalyptus tereticornis - acclim")
ax1.set_xlabel("Leaf temperature ($^{\circ}\mathrm{C}$)", position=(1.1,0.5))
ax1.set_ylabel("$V_{\mathrm{cmax}}$ ($\mu$mol m$^{-2}$ s$^{-1}$)")
ax1.legend(numpoints=1, loc="best")

ax2.plot(Tleaf, Jmax_cable, label="CABLE EBF")
ax2.plot(Tleaf, Jmax_euc, label="Eucalyptus tereticornis")
ax2.plot(Tleaf, Jmax_euc_acclim, label="Eucalyptus tereticornis - acclim")
ax2.set_ylabel("$J_{\mathrm{max}}$ ($\mu$mol m$^{-2}$ s$^{-1}$)")
fig.savefig("CABLE_vs_HIE.png", dpi=150)
