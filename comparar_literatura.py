#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
 comparar_literatura.py   (version 2, corregida)
 Compara la simulacion contra los datos experimentales publicados.
.

 USO
 ---
     python3 comparar_literatura.py            # ultimo cuadro
     python3 comparar_literatura.py 7          # cuadro 7
     python3 comparar_literatura.py todos      # evolucion temporal
=============================================================================
"""

import glob
import os
import re
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# =========================================================================
#  CONFIGURACION
# =========================================================================

GUARDAR = True
MOSTRAR = True
DPI     = 130

NX_FORZADO, NY_FORZADO, NYS_FORZADO = None, None, None
NX = NY = NYS = 0

# Zona del grano [m]. Debe coincidir con mod_malla.f90.
X_INI_GRANO = 0.025
X_FIN_GRANO = 0.080

G_OX_SIM = 14.0          # kg/(m2 s), G_ox_fijo de mod_malla.f90

# Borde de capa limite: desviacion relativa de densidad respecto del gas
# libre. 0.01 = 1 %, el valor habitual.
UMBRAL_RHO = 0.01

JENS_FRAC_CLASICOS = 0.33
JENS_FRAC_LICUEFAC = 0.45

#      (test, G_ox [kg/m2 s], r [mm/s], p_camara [bar])
BIANCHI = [
    ( 2,  27.68, 0.47, 17.07), ( 4,  87.14, 0.64, 15.84),
    ( 9,  28.41, 0.47, 22.91), (10, 126.00, 0.82, 20.52),
    (13,  28.19, 0.38, 10.22), (18,  78.07, 0.55, 11.94),
    (19,  54.03, 0.50, 11.26),
]

RE_LAM, RE_TURB = 2300.0, 4000.0      # transicion en ducto
MU_REF = 3.0e-5                        # viscosidad de referencia [Pa s]

ESPECIES = ['C2H4', 'O2', 'CO2', 'H2O', 'CO', 'H2', 'H', 'O', 'OH']
COL_OH = 2 + ESPECIES.index('OH')      # = 10


# =========================================================================
#  LECTURA
# =========================================================================

def fatal(msg):
    print('\n' + '=' * 70)
    print(' ERROR')
    print('=' * 70)
    for l in str(msg).split('\n'):
        print(' ' + l)
    print('=' * 70 + '\n')
    raise SystemExit(1)


def numero_de(nombre, prefijo):
    m = re.search(rf'{prefijo}_(\d+)\.dat$', os.path.basename(nombre))
    return int(m.group(1)) if m else -1


def listar(prefijo):
    a = [f for f in glob.glob(f'{prefijo}_*.dat') if numero_de(f, prefijo) >= 0]
    return sorted(a, key=lambda f: numero_de(f, prefijo))


def detectar_malla():
    global NX, NY, NYS
    arch = listar('densidad') or listar('temperatura')
    if not arch:
        fatal("No hay archivos densidad_XXXX.dat ni temperatura_XXXX.dat.\n"
              f"Carpeta actual: {os.getcwd()}")
    d = np.loadtxt(arch[0])
    n, x = d.shape[0], d[:, 0]
    nx = next((i for i in range(1, n) if x[i] <= x[i - 1]), n)
    NX = NX_FORZADO or nx
    NY = NY_FORZADO or n // NX
    if NX * NY != n:
        fatal(f"La malla no cuadra: {n} nodos, NX*NY = {NX}*{NY}.")
    if NYS_FORZADO:
        NYS = NYS_FORZADO
    else:
        obs = listar('obstaculo')
        if obs:
            f = np.loadtxt(obs[0])[:, 2].reshape(NY, NX)
            fil = np.where(f.max(axis=1) > 0.5)[0]
            NYS = int(fil.max()) + 1 if fil.size else NY - 100
        else:
            NYS = NY - 100
    print(f"Malla: {NX} x {NY}, superficie en j = {NYS}")


def cuadros_disponibles():
    ind = [numero_de(a, 'densidad') for a in listar('densidad')]
    ind = [i for i in ind if i != 9999]
    if not ind:
        fatal("No hay archivos densidad_XXXX.dat en esta carpeta.\n"
              f"Carpeta actual: {os.getcwd()}")
    return sorted(ind)


def _col(prefijo, idx, col):
    f = f'{prefijo}_{idx:04d}.dat'
    if not os.path.exists(f):
        fatal(f"Falta '{f}'.")
    return np.loadtxt(f)[:, col].reshape(NY, NX)


def leer_tiempo(idx):
    f = f'tiempo_{idx:04d}.dat'
    if not os.path.exists(f):
        return None
    try:
        return float(np.loadtxt(f))
    except Exception:
        return None


def leer_rdot(idx):
    f = f'superficie_{idx:04d}.dat'
    if not os.path.exists(f):
        return None, None
    d = np.loadtxt(f)
    return d[:, 0], d[:, 3]           # x [m], rdot [mm/s]


# =========================================================================
#  EL CALCULO CORREGIDO
# =========================================================================

def analizar(idx):
    """delta por densidad, altura de llama por OH, Reynolds de ducto."""
    m = np.loadtxt(f'malla_nodos_{idx:04d}.dat')
    X = m[:, 0].reshape(NY, NX)
    Y = m[:, 1].reshape(NY, NX)
    rho = _col('densidad', idx, 2)
    T = _col('temperatura', idx, 2)
    u = _col('campo', idx, 0)
    OH = _col('especies', idx, COL_OH)

    js = NYS - 1
    Yg, rg, Tg = Y[js:, :], rho[js:, :], T[js:, :]
    ug, OHg = u[js:, :], OH[js:, :]
    x = X[0, :]

    delta = np.full(NX, np.nan)
    y_fl = np.zeros(NX)
    frac = np.full(NX, np.nan)
    ReDh = np.zeros(NX)
    Tmax = np.zeros(NX)
    ohmx = np.zeros(NX)
    Hp = np.zeros(NX)

    for i in range(NX):
        yy = Yg[:, i] - Yg[0, i]       # altura sobre la superficie
        rr = rg[:, i]
        Hp[i] = yy[-1]
        rho_inf = rr[-1]               # gas libre, en el techo

        # --- borde de capa limite por GRADIENTE DE DENSIDAD -------------
        if rho_inf > 1e-12:
            dev = np.abs(rr - rho_inf) / rho_inf
            k = np.where(dev > UMBRAL_RHO)[0]
            if k.size and yy[k[-1]] > 1e-9:
                delta[i] = yy[k[-1]]

        # --- altura de llama: maximo de OH ------------------------------
        kf = int(np.argmax(OHg[:, i]))
        ohmx[i] = OHg[kf, i]
        y_fl[i] = yy[kf]
        if np.isfinite(delta[i]) and delta[i] > 0:
            frac[i] = y_fl[i] / delta[i]

        Tmax[i] = Tg[:, i].max()

        # --- Reynolds por DIAMETRO HIDRAULICO: Dh = 2*Hp ----------------
        ReDh[i] = rho_inf * ug[:, i].max() * (2.0 * Hp[i]) / MU_REF

    return dict(x=x, delta=delta, y_llama=y_fl, frac=frac, Tmax=Tmax,
                Re=ReDh, oh=ohmx, Hp=Hp)


def mascara_grano(c):
    return ((c['x'] >= X_INI_GRANO) & (c['x'] <= X_FIN_GRANO)
            & np.isfinite(c['delta']) & (c['oh'] > 1e-6))


# =========================================================================
#  GRAFICOS
# =========================================================================

def graficar_cuadro(idx):
    c = analizar(idx)
    t = leer_tiempo(idx)
    g = mascara_grano(c)
    gg = (c['x'] >= X_INI_GRANO) & (c['x'] <= X_FIN_GRANO)
    xg_, rg_ = leer_rdot(idx)

    if not g.any():
        print(f"\n[AVISO] el cuadro {idx} no tiene llama sobre el grano.\n"
              "        Es anterior al encendido. Proba un cuadro posterior.\n")

    fig, axs = plt.subplots(2, 2, figsize=(13, 8), dpi=DPI)
    sub = f'cuadro {idx}' + (f'   t = {t:.4f} s' if t is not None else '')
    fig.suptitle(f'Comparacion contra la literatura  -  {sub}', fontsize=13)

    # --- (1) perfiles, SOLO sobre el grano ------------------------------
    ax = axs[0, 0]
    ax.plot(c['x'][gg] * 1e3, c['delta'][gg] * 1e3, color='tab:blue', lw=1.8,
            label=r'capa limite $\delta_\rho$ (gradiente de densidad)')
    if g.any():
        ax.plot(c['x'][g] * 1e3, c['y_llama'][g] * 1e3, color='tab:red',
                lw=1.8, label='altura de llama (max. de OH)')
    ax.plot(c['x'][gg] * 1e3, c['Hp'][gg] * 1e3, color='0.6', ls=':', lw=1.2,
            label='alto del puerto')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('altura sobre la superficie (mm)')
    ax.set_title('Capa limite y llama, solo sobre el grano')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # --- (2) la comparacion fuerte --------------------------------------
    ax = axs[0, 1]
    if g.any():
        ax.plot(c['x'][g] * 1e3, c['frac'][g], color='tab:red', lw=1.8,
                label='simulacion')
        fm = float(np.nanmean(c['frac'][g]))
        ax.axhline(fm, color='tab:red', ls=':', lw=1.2,
                   label=f'media simulacion = {fm:.2f}')
    ax.axhline(JENS_FRAC_CLASICOS, color='k', ls='--', lw=1.8,
               label=f'Jens 2019, HDPE/PMMA = {JENS_FRAC_CLASICOS}')
    ax.axhline(JENS_FRAC_LICUEFAC, color='0.5', ls='-.', lw=1.2,
               label=f'Jens 2019, parafina = {JENS_FRAC_LICUEFAC}')
    ax.set_ylim(0, 1.2)
    ax.set_xlim(X_INI_GRANO * 1e3, X_FIN_GRANO * 1e3)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel(r'$y_{llama}\,/\,\delta_\rho$')
    ax.set_title('CONFIANZA ALTA: posicion de la llama en la capa limite')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(alpha=0.3)

    # --- (3) regresion ---------------------------------------------------
    ax = axs[1, 0]
    G = np.array([b[1] for b in BIANCHI])
    R = np.array([b[2] for b in BIANCHI])
    n, lna = np.polyfit(np.log(G), np.log(R), 1)
    a_fit = np.exp(lna)
    Gs = np.logspace(np.log10(10), np.log10(140), 100)
    ax.plot(Gs, a_fit * Gs ** n, 'k--', lw=1.3,
            label=f'ajuste a Bianchi:  $r={a_fit:.3f}\\,G^{{{n:.2f}}}$')
    ax.plot(G, R, 'ks', ms=7, label='Bianchi 2019 (10-23 bar, cilindrico)')
    r_med = r_max = None
    if xg_ is not None:
        sel = (xg_ >= X_INI_GRANO) & (xg_ <= X_FIN_GRANO)
        if sel.any():
            r_med, r_max = float(rg_[sel].mean()), float(rg_[sel].max())
            ax.plot([G_OX_SIM], [r_med], 'r*', ms=15,
                    label=f'simulacion, MEDIA: {r_med:.3f} mm/s')
            ax.plot([G_OX_SIM], [r_max], 'rv', ms=9, mfc='none',
                    label=f'simulacion, MAXIMA: {r_max:.3f} mm/s')
    ax.axvspan(10, G.min(), color='orange', alpha=0.12)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'$G_{ox}$  (kg/m$^2$s)')
    ax.set_ylabel(r'$\dot r$  (mm/s)')
    ax.set_title('CONFIANZA MEDIA: regresion (naranja = extrapolado)')
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, which='both')

    # --- (4) Reynolds de DUCTO ------------------------------------------
    ax = axs[1, 1]
    ax.plot(c['x'][gg] * 1e3, c['Re'][gg], color='tab:green', lw=1.8,
            label=r'$Re_{D_h}$ simulado  ($D_h = 2H_p$)')
    ax.axhspan(RE_LAM, RE_TURB, color='orange', alpha=0.25,
               label='transicion en ducto (2300-4000)')
    ax.axhline(RE_TURB, color='k', ls='--', lw=1.2)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel(r'$Re_{D_h}$')
    ax.set_yscale('log')
    ax.set_title('CONTEXTO: por encima de 4000 el flujo real es turbulento')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which='both')

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    if GUARDAR:
        nom = f'comparacion_{idx:04d}.png'
        fig.savefig(nom, dpi=DPI, bbox_inches='tight')
        print(f"Guardado: {os.path.abspath(nom)}")

    print()
    print('=' * 70)
    print(f' RESUMEN  -  cuadro {idx}' +
          (f'   t = {t:.4f} s' if t is not None else ''))
    print('=' * 70)
    if g.any():
        print(f" delta (densidad), medio   : {np.nanmean(c['delta'][g])*1e3:7.2f} mm")
        print(f" altura de llama, media    : {np.nanmean(c['y_llama'][g])*1e3:7.2f} mm")
        print(f" y_llama / delta           : {np.nanmean(c['frac'][g]):7.3f}"
              f"   <-- Jens mide {JENS_FRAC_CLASICOS}")
        print(f" T maxima                  : {c['Tmax'][g].max():7.0f} K")
        print(f" columnas con llama        : {int(g.sum())} de {int(gg.sum())}"
              f" sobre el grano")
    else:
        print(" sin llama sobre el grano en este cuadro")
    print(f" Re_Dh medio sobre el grano: {np.nanmean(c['Re'][gg]):7.0f}"
          f"   (turbulento por encima de {RE_TURB:.0f})")
    if r_med is not None:
        print(f" r medio / maximo          : {r_med:7.3f} / {r_max:.3f} mm/s"
              f"   <-- Bianchi extrapolado {a_fit*G_OX_SIM**n:.3f}")
    print('=' * 70 + '\n')

    plt.show() if MOSTRAR else plt.close(fig)


def graficar_evolucion():
    ind = cuadros_disponibles()
    print(f"Leyendo {len(ind)} cuadros...")
    ts, fr, dd, rr, nn = [], [], [], [], []
    for i in ind:
        c = analizar(i)
        g = mascara_grano(c)
        t = leer_tiempo(i)
        if t is None or not g.any():
            continue
        ts.append(t)
        fr.append(float(np.nanmean(c['frac'][g])))
        dd.append(float(np.nanmean(c['delta'][g])) * 1e3)
        nn.append(int(g.sum()))
        xg_, rg_ = leer_rdot(i)
        if xg_ is not None:
            sel = (xg_ >= X_INI_GRANO) & (xg_ <= X_FIN_GRANO)
            rr.append(float(rg_[sel].mean()) if sel.any() else np.nan)
        else:
            rr.append(np.nan)

    if not ts:
        fatal("Ningun cuadro tiene llama sobre el grano.\n"
              "La simulacion no llego al encendido todavia.")

    fig, axs = plt.subplots(4, 1, figsize=(9, 11), dpi=DPI, sharex=True)
    fig.suptitle('Evolucion temporal', fontsize=13)

    axs[0].plot(ts, fr, 'o-', color='tab:red', ms=3)
    axs[0].axhline(JENS_FRAC_CLASICOS, color='k', ls='--',
                   label=f'Jens 2019 = {JENS_FRAC_CLASICOS}')
    axs[0].set_ylabel(r'$y_{llama}/\delta_\rho$')
    axs[0].legend(fontsize=8); axs[0].grid(alpha=0.3)

    axs[1].plot(ts, dd, 'o-', color='tab:blue', ms=3)
    axs[1].set_ylabel(r'$\delta_\rho$  (mm)'); axs[1].grid(alpha=0.3)

    axs[2].plot(ts, rr, 'o-', color='tab:green', ms=3)
    G = np.array([b[1] for b in BIANCHI]); R = np.array([b[2] for b in BIANCHI])
    n, lna = np.polyfit(np.log(G), np.log(R), 1)
    axs[2].axhline(np.exp(lna) * G_OX_SIM ** n, color='k', ls='--',
                   label='Bianchi extrapolado')
    axs[2].set_ylabel(r'$\dot r$  (mm/s)')
    axs[2].legend(fontsize=8); axs[2].grid(alpha=0.3)

    axs[3].plot(ts, nn, 'o-', color='tab:purple', ms=3)
    axs[3].set_ylabel('columnas\ncon llama')
    axs[3].set_xlabel('tiempo (s)'); axs[3].grid(alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    if GUARDAR:
        fig.savefig('comparacion_evolucion.png', dpi=DPI, bbox_inches='tight')
        print(f"Guardado: {os.path.abspath('comparacion_evolucion.png')}")
    plt.show() if MOSTRAR else plt.close(fig)


# =========================================================================
if __name__ == '__main__':
    detectar_malla()
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    if arg == 'todos':
        graficar_evolucion()
    elif arg is not None:
        graficar_cuadro(int(arg))
    else:
        graficar_cuadro(cuadros_disponibles()[-1])
