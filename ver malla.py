#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ver_malla_sola.py
Dibuja SOLO la malla que escribe el solver. No necesita nada mas que el
archivo malla_nodos; si encuentra el obstaculo del mismo instante, separa
solido de gas por color.

    python3 ver_malla_sola.py                 # malla_nodos.dat (la inicial)
    python3 ver_malla_sola.py 0012            # malla_nodos_0012.dat
    python3 ver_malla_sola.py 9998            # el volcado de aviso
    python3 ver_malla_sola.py --zoom 22 32    # acercar en x, en mm
    python3 ver_malla_sola.py --salto 1       # dibujar TODAS las lineas
    python3 ver_malla_sola.py --png           # guardar en vez de mostrar

Imprime ademas las cifras de calidad: espaciados, nodos por rampa,
pendiente maxima y angulo minimo entre lineas de malla.
"""
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

C_SOL = "#8c5a2b"      # solido
C_GAS = "#1f4e79"      # gas
C_SUP = "#000000"      # superficie


def leer(nombre):
    d = np.loadtxt(nombre)
    x = d[:, 0]
    # el solver escribe j por fuera e i por dentro: x crece y se reinicia
    NX = int(np.argmax(np.diff(x) < 0) + 1)
    if NX <= 1:                       # una sola fila, o archivo raro
        NX = len(x)
    NY, resto = divmod(len(x), NX)
    if resto:
        raise SystemExit(f"{nombre}: {len(x)} filas no son multiplo de NX={NX}")
    return d[:, 0].reshape(NY, NX), d[:, 1].reshape(NY, NX), NX, NY


def buscar_obstaculo(sufijo, NX, NY):
    nom = f"obstaculo{sufijo}.dat"
    if not os.path.exists(nom):
        return None
    try:
        ob = np.loadtxt(nom)[:, 2].reshape(NY, NX).astype(int)
    except Exception:
        return None
    filas = np.where((ob == 1).any(axis=1))[0]
    return int(filas.max()) + 1 if len(filas) else None     # NYS, base 1


def calidad(X, Y, NYS, NX, NY):
    x = X[0]
    dx = np.diff(x)
    e = Y[NYS - 1]                                    # superficie
    M = np.zeros(NX)
    M[1:-1] = (e[2:] - e[:-2]) / (x[2:] - x[:-2])
    print(f"  malla          : NX={NX}  NY={NY}  NYS={NYS}")
    print(f"  dx             : {dx.min()*1e6:7.1f} a {dx.max()*1e6:7.1f} um")
    if NYS < NY:
        dy1 = np.diff(Y[NYS - 1:, :], axis=0)[0]
        print(f"  1a celda de gas: {dy1.min()*1e6:7.1f} a {dy1.max()*1e6:7.1f} um")
    print(f"  espesor        : {e.min()*1e3:7.2f} a {e.max()*1e3:7.2f} mm")
    print(f"  |dY/dx| maximo : {np.abs(M).max():7.2f}"
          f"   -> angulo minimo entre lineas {np.degrees(np.arctan(1/max(np.abs(M).max(),1e-9))):.2f} grados")
    # rampas: columnas donde el espesor no esta ni en el minimo ni en el maximo
    tol = 0.02 * (e.max() - e.min())
    media = np.where((e > e.min() + tol) & (e < e.max() - tol))[0]
    if len(media):
        cortes = np.split(media, np.where(np.diff(media) > 1)[0] + 1)
        for k, c in enumerate(cortes, 1):
            print(f"  rampa {k}        : columnas {c[0]+1}-{c[-1]+1}"
                  f"  ({x[c[0]]*1e3:.2f} a {x[c[-1]]*1e3:.2f} mm)"
                  f"  {len(c)} nodos")


def main():
    arg = [a for a in sys.argv[1:] if not a.startswith("--")]
    png = "--png" in sys.argv
    salto = 4
    if "--salto" in sys.argv:
        salto = int(sys.argv[sys.argv.index("--salto") + 1])
    zoom = None
    if "--zoom" in sys.argv:
        k = sys.argv.index("--zoom")
        zoom = (float(sys.argv[k + 1]), float(sys.argv[k + 2]))
        salto = 1 if "--salto" not in sys.argv else salto

    if arg:
        sufijo = "_" + arg[0]
        nombre = f"malla_nodos{sufijo}.dat"
    else:
        sufijo, nombre = "", "malla_nodos.dat"
    if not os.path.exists(nombre):
        hay = sorted(glob.glob("malla_nodos*.dat"))
        raise SystemExit(f"No encuentro {nombre}.\nEn esta carpeta hay: "
                         + (", ".join(hay) if hay else "ninguno"))

    X, Y, NX, NY = leer(nombre)
    NYS = buscar_obstaculo(sufijo, NX, NY)
    if NYS is None:
        NYS = 41
        print(f"(sin obstaculo{sufijo}.dat: asumo NYS={NYS})")

    print(f"{nombre}:")
    calidad(X, Y, NYS, NX, NY)

    fig, ax = plt.subplots(figsize=(11, 5) if zoom is None else (7, 7))
    # lineas eta = cte (horizontales en el espacio computacional)
    for j in range(NY):
        if j != NYS - 1 and j % salto:
            continue
        c = C_SOL if j < NYS - 1 else C_GAS
        lw = 0.35
        if j == NYS - 1:
            c, lw = C_SUP, 1.2
        ax.plot(X[j] * 1e3, Y[j] * 1e3, lw=lw, color=c)
    # lineas xi = cte (verticales)
    for i in range(0, NX, salto):
        ax.plot(X[:NYS, i] * 1e3, Y[:NYS, i] * 1e3, lw=0.35, color=C_SOL)
        ax.plot(X[NYS - 1:, i] * 1e3, Y[NYS - 1:, i] * 1e3, lw=0.35, color=C_GAS)

    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_aspect("equal")
    if zoom:
        ax.set_xlim(*zoom)
        m = (X[0] * 1e3 >= zoom[0]) & (X[0] * 1e3 <= zoom[1])
        ax.set_ylim(0, Y[:, m].max() * 1e3 * 1.08)
        ax.set_title(f"{nombre}  --  zoom {zoom[0]:g} a {zoom[1]:g} mm"
                     f"  (1 de cada {salto} lineas)")
    else:
        ax.set_xlim(X.min() * 1e3, X.max() * 1e3)
        ax.set_ylim(0, Y.max() * 1e3)
        ax.set_title(f"{nombre}  --  {NX}x{NY} nodos"
                     f"  (1 de cada {salto} lineas)")
    fig.tight_layout()

    if png:
        sal = f"malla{sufijo or '_inicial'}.png"
        fig.savefig(sal, dpi=140)
        print(f"Guardado {sal}")
    plt.show()


if __name__ == "__main__":
    main()