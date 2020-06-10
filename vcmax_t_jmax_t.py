import numpy as np
import matplotlib.pyplot as plt

def xt3(k25, Tleaf, Ea, deltaS, Hd):
    RGAS = 8.31429958
    TREFK = 25.0 + 273.15 #298.200012

    coef = 1.0 + np.exp((deltaS * TREFK - Hd) / (RGAS * TREFK))
    num = coef * np.exp((Ea / (RGAS * TREFK)) * (1. - TREFK / Tleaf))
    den = 1.0 + np.exp((deltaS * Tleaf - Hd) / (RGAS * Tleaf))

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

Vcmax25 = 86
Jmax25 = 138
Tleaf = 25.0
TleafK = 25.0 + 273.15

Eaj = 50300.0
Eav = 73637.0
deltaSj = 495.0
deltaSv = 486.0
Hdv = 149252.0
Hdj = 152044.0

Tleaf = np.arange(1.0, 40.0)
TleafK = Tleaf + 273.15
Vcmax = xt3(Vcmax25, TleafK, Eav, deltaSv, Hdv)
Jcmax = xt3(Jmax25, TleafK, Eaj, deltaSj, Hdj)


Eaj = 32292.6319602695
Eav = 66386.2007802811
deltaSj = 638.055941746436
deltaSv = 639.602873767181
Hdv = 200000.0
Hdj = 200000.0

Vcmax_me = peaked_arrh(Vcmax25, Eav, TleafK, deltaSv, Hdv)
Jcmax_me = peaked_arrh(Jmax25, Eaj, TleafK, deltaSj, Hdj)

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

ax1.plot(Tleaf, Vcmax, label="Australia's LSM")
ax1.plot(Tleaf, Vcmax_me, label="Eucalyptus tereticornis")
ax1.set_xlabel("Leaf temperature ($^{\circ}\mathrm{C}$)", position=(1.1,0.5))
ax1.set_ylabel("$V_{\mathrm{cmax}}$ ($\mu$mol m$^{-2}$ s$^{-1}$)")
ax1.legend(numpoints=1, loc="best")

ax2.plot(Tleaf, Jcmax)
ax2.plot(Tleaf, Jcmax_me)
ax2.set_ylabel("$J_{\mathrm{max}}$ ($\mu$mol m$^{-2}$ s$^{-1}$)")
fig.savefig("CABLE_vs_HIE.png", dpi=150)
