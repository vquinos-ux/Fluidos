#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import os
import re
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation


NX_FORZADO, NY_FORZADO, NYS_FORZADO = None, None, None
NX, NY, NYS = 0, 0, 0

SEMILLAS = 'entrada'
N_SEMILLAS = 22
DENSIDAD   = 1.3
GROSOR_LINEA = 1.0

N_FLECHAS_X = 46
N_FLECHAS_Y = 16

NX_REG, NY_REG = 420, 190

MOSTRAR_SOLIDO = True

ASPECTO   = 'equal'
PERCENTIL = (1.0, 99.0)
FPS, DPI  = 8, 130
GUARDAR   = 'gif'
MOSTRAR   = True

CADA_N    = 1

ESPECIES = ['C2H4', 'O2', 'CO2', 'H2O', 'CO', 'H2', 'H', 'O', 'OH']


def fatal(msg):
    print('\n' + '=' * 68)
    print(' ERROR')
    print('=' * 68)
    for l in str(msg).split('\n'):
        print(' ' + l)
    print('=' * 68 + '\n')
    raise SystemExit(1)


def cargar(nombre):
    if not os.path.exists(nombre):
        fatal(f"No existe '{nombre}'.\nCarpeta actual: {os.getcwd()}")
    try:
        return np.loadtxt(nombre)
    except ValueError:
        patron = re.compile(r'(\d)([+-]\d\d\d)')
        filas = []
        with open(nombre) as fh:
            for linea in fh:
                if linea.strip() and not linea.lstrip().startswith('#'):
                    filas.append([float(v) for v in
                                  patron.sub(r'\1E\2', linea).split()])
        return np.array(filas)


def numero_de(n, p):
    m = re.search(rf'{p}_(\d+)\.dat$', os.path.basename(n))
    return int(m.group(1)) if m else -1


def listar(p):
    a = [f for f in glob.glob(f'{p}_*.dat') if numero_de(f, p) >= 0]
    return sorted(a, key=lambda f: numero_de(f, p))


def detectar_malla():
    global NX, NY, NYS
    fuente = None
    for p in ('temperatura', 'especies', 'densidad', 'obstaculo'):
        if listar(p):
            fuente = listar(p)[0]
            break
    if fuente is None and os.path.exists('malla_nodos.dat'):
        fuente = 'malla_nodos.dat'
    if fuente is None:
        fatal("No hay archivos con coordenadas (temperatura_, especies_,\n"
              "densidad_, obstaculo_ o malla_nodos.dat).\n"
              f"Carpeta actual: {os.getcwd()}")

    d = cargar(fuente)
    n, x = d.shape[0], d[:, 0]
    nx = next((i for i in range(1, n) if x[i] <= x[i - 1]), n)
    NX = NX_FORZADO or nx
    NY = NY_FORZADO or n // NX
    if NX * NY != n:
        fatal(f"No cuadra la malla de '{fuente}': {n} nodos, "
              f"NX*NY = {NX}*{NY} = {NX*NY}.")

    if NYS_FORZADO:
        NYS = NYS_FORZADO
    else:
        obs = listar('obstaculo')
        if obs:
            f = cargar(obs[0])[:, 2].reshape(NY, NX)
            fil = np.where(f.max(axis=1) > 0.5)[0]
            NYS = int(fil.max()) + 1 if fil.size else NY - 100
        else:
            NYS = NY - 100
    print(f"Malla: {NX} x {NY} = {n} nodos, superficie en j = {NYS}   [{fuente}]")


def leer_malla(idx):
    nom = f'malla_nodos_{idx:04d}.dat'
    if not os.path.exists(nom):
        nom = 'malla_nodos.dat'
    d = cargar(nom)
    return d[:, 0].reshape(NY, NX), d[:, 1].reshape(NY, NX)


def leer_tiempo(idx):
    nom = f'tiempo_{idx:04d}.dat'
    if not os.path.exists(nom):
        return None
    try:
        return float(np.loadtxt(nom))
    except Exception:
        return None


def leer_escalar(idx, campo):
    if campo in ('lineas', 'flechas', 'rapidez'):
        d = cargar(f'campo_{idx:04d}.dat')
        return d[:, 3].reshape(NY, NX), 'Rapidez |V|  (m/s)', 'turbo'
    if campo == 'temperatura':
        d = cargar(f'temperatura_{idx:04d}.dat')
        return d[:, 2].reshape(NY, NX), 'Temperatura (K)', 'inferno'
    if campo == 'oxidante':
        d = cargar(f'especies_{idx:04d}.dat')
        return (d[:, 2 + ESPECIES.index('O2')].reshape(NY, NX),
                'Fraccion masica de O$_2$', 'Blues')
    if campo == 'combustible':
        d = cargar(f'especies_{idx:04d}.dat')
        return (d[:, 2 + ESPECIES.index('C2H4')].reshape(NY, NX),
                'Fraccion masica de C$_2$H$_4$', 'YlOrBr')
    fatal(f"Modo desconocido: '{campo}'.\n"
          "Opciones: lineas, flechas, oxidante, combustible, temperatura")


def leer_velocidad(idx):
    d = cargar(f'campo_{idx:04d}.dat')
    return d[:, 0].reshape(NY, NX), d[:, 1].reshape(NY, NX)


def a_regular(X, Y, campos):
    ny_loc, nx_loc = X.shape
    x_col = X[0, :]

    y_reg = np.linspace(Y.min(), Y.max(), NY_REG)
    x_reg = np.linspace(x_col.min(), x_col.max(), NX_REG)

    tmp = [np.empty((NY_REG, nx_loc)) for _ in campos]
    for i in range(nx_loc):
        yc = Y[:, i]
        for t, c in zip(tmp, campos):
            t[:, i] = np.interp(y_reg, yc, c[:, i])

    salida = []
    for t in tmp:
        o = np.empty((NY_REG, NX_REG))
        for j in range(NY_REG):
            o[j, :] = np.interp(x_reg, x_col, t[j, :])
        salida.append(o)

    sup = np.interp(x_reg, x_col, Y[0, :])
    bajo = y_reg[:, None] < sup[None, :]

    return x_reg, y_reg, salida, bajo, sup


def animar(modo='lineas'):
    detectar_malla()

    archivos = listar('campo')
    if not archivos:
        fatal("No hay archivos campo_*.dat en esta carpeta.\n"
              f"Carpeta actual: {os.getcwd()}")
    indices = [numero_de(a, 'campo') for a in archivos]
    if 9999 in indices:
        print("   [AVISO] existe el cuadro 9999: la simulacion divergio.")
    print(f"Cuadros: {len(indices)} (del {min(indices)} al {max(indices)})")
    if CADA_N > 1:
        sel = indices[::CADA_N]
        if sel[-1] != indices[-1]:
            sel.append(indices[-1])
        print(f"CADA_N = {CADA_N}: se usan {len(sel)} de {len(indices)} cuadros")
        indices = sel
    print(f"Modo: {modo}   Interpolando a rejilla {NX_REG} x {NY_REG} ...")

    js = NYS - 1
    cuadros = []
    for k, idx in enumerate(indices, 1):
        X, Y = leer_malla(idx)
        u, v = leer_velocidad(idx)
        C, etiqueta, paleta = leer_escalar(idx, modo)

        X, Y = X[js:, :], Y[js:, :]
        u, v, C = u[js:, :], v[js:, :], C[js:, :]

        xr, yr, (ur, vr, Cr), bajo, sup = a_regular(X, Y, [u, v, C])
        ur[bajo] = 0.0
        vr[bajo] = 0.0
        Cr = np.ma.masked_array(Cr, mask=bajo)

        cuadros.append((xr, yr, ur, vr, Cr, sup, leer_tiempo(idx), idx))
        print(f"   {k}/{len(indices)}")

    todos = np.ma.compressed(np.ma.concatenate([np.ma.ravel(c[4])
                                                for c in cuadros]))
    todos = todos[np.isfinite(todos)]
    v_min = float(np.percentile(todos, PERCENTIL[0]))
    v_max = float(np.percentile(todos, PERCENTIL[1]))
    if v_max - v_min < 1e-12:
        v_max = v_min + 1.0
    print(f"-> escala de color: {v_min:.4g} a {v_max:.4g}")

    _, etiqueta, paleta = leer_escalar(indices[0], modo)
    norma = matplotlib.colors.Normalize(v_min, v_max)

    x_lo, x_hi = cuadros[0][0][0] * 1e3, cuadros[0][0][-1] * 1e3
    y_hi = max(c[1][-1] for c in cuadros) * 1e3
    if MOSTRAR_SOLIDO:
        y_lo = 0.0
    else:
        y_lo = min(c[1][0] for c in cuadros) * 1e3

    fig, ax = plt.subplots(figsize=(12, 5.2), dpi=DPI)
    sm = plt.cm.ScalarMappable(norm=norma, cmap=paleta)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label=etiqueta, pad=0.02, fraction=0.035)

    sm2 = None
    if modo == 'flechas':
        v_tot = np.ma.compressed(np.ma.concatenate(
            [np.ma.ravel(np.hypot(c[2], c[3])) for c in cuadros]))
        vq_max = float(np.percentile(v_tot, 99.0))
        norma_q = matplotlib.colors.Normalize(0.0, vq_max)
        sm2 = plt.cm.ScalarMappable(norm=norma_q, cmap='turbo')
        sm2.set_array([])

    def dibujar(n):
        xr, yr, ur, vr, Cr, sup, t, idx = cuadros[n]
        ax.clear()

        ax.pcolormesh(xr * 1e3, yr * 1e3, Cr, cmap=paleta,
                      norm=norma, shading='gouraud', zorder=1)

        if modo == 'flechas':
            ix = np.linspace(2, len(xr) - 3, N_FLECHAS_X).astype(int)
            iy = np.linspace(2, len(yr) - 3, N_FLECHAS_Y).astype(int)
            xq, yq = np.meshgrid(xr[ix] * 1e3, yr[iy] * 1e3)
            uq, vq = ur[np.ix_(iy, ix)], vr[np.ix_(iy, ix)]
            mag = np.hypot(uq, vq)
            ok = mag > 1e-6
            un = np.where(ok, uq / np.where(ok, mag, 1.0), 0.0)
            vn = np.where(ok, vq / np.where(ok, mag, 1.0), 0.0)
            ax.quiver(xq[ok], yq[ok], un[ok], vn[ok], mag[ok],
                      cmap='turbo', norm=sm2.norm, pivot='mid',
                      scale=N_FLECHAS_X * 1.15, width=0.0022,
                      headwidth=3.5, headlength=4.0, zorder=2)
        else:
            if SEMILLAS == 'entrada':
                y0 = np.linspace(sup[2] * 1.001, yr[-1] * 0.999, N_SEMILLAS)
                puntos = np.column_stack([np.full_like(y0, xr[2]) * 1e3,
                                          y0 * 1e3])
                ax.streamplot(xr * 1e3, yr * 1e3, ur, vr, color='white',
                              linewidth=GROSOR_LINEA, arrowsize=1.1,
                              start_points=puntos, zorder=2,
                              integration_direction='both')
            else:
                ax.streamplot(xr * 1e3, yr * 1e3, ur, vr, color='white',
                              linewidth=GROSOR_LINEA, arrowsize=1.1,
                              density=DENSIDAD, zorder=2)

        ax.fill_between(xr * 1e3, y_lo, sup * 1e3,
                        color='0.30', zorder=4, edgecolor='none')
        ax.plot(xr * 1e3, sup * 1e3, color='k', lw=1.6, zorder=5)
        if MOSTRAR_SOLIDO:
            ax.text(0.5 * (x_lo + x_hi), 0.5 * (y_lo + sup.mean() * 1e3),
                    'COMBUSTIBLE  (HDPE)', color='white', zorder=6,
                    ha='center', va='center', fontsize=11, alpha=0.75)

        titulo = f'{etiqueta}  -  cuadro {idx}'
        if t is not None:
            titulo += f'   t = {t:.4f} s'
        ax.set_title(titulo)
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_aspect(ASPECTO)
        return []

    ani = animation.FuncAnimation(fig, dibujar, frames=len(cuadros),
                                  interval=1000 / FPS, blit=False)

    if GUARDAR == 'gif':
        nom = f'flujo_{modo}.gif'
        print(f"Guardando '{nom}' ...")
        ani.save(nom, writer='pillow', fps=FPS, dpi=DPI)
        print(f"Listo: {os.path.abspath(nom)}")
    elif GUARDAR == 'mp4':
        nom = f'flujo_{modo}.mp4'
        print(f"Guardando '{nom}' ...")
        try:
            ani.save(nom, writer='ffmpeg', fps=FPS, dpi=DPI)
            print(f"Listo: {os.path.abspath(nom)}")
        except Exception as e:
            print(f"   [AVISO] no se pudo guardar el mp4: {e}")
            print("   Falta ffmpeg. Usa GUARDAR = 'gif'.")

    if MOSTRAR:
        plt.show()
    else:
        plt.close(fig)


if __name__ == '__main__':
    animar(sys.argv[1] if len(sys.argv) > 1 else 'lineas')
