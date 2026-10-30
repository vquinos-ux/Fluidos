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


NX_FORZADO  = None
NY_FORZADO  = None
NYS_FORZADO = None

NX, NY, NYS = 0, 0, 0
J_SUP   = 0
FORMATO = ''

REGION  = 'gas'
METODO  = 'pcolor'
ASPECTO = 'auto'

PALETA_T   = 'inferno'
PALETA_Y   = 'viridis'
PALETA_VEL = 'turbo'

NIVELES   = 100
FPS       = 12
DPI       = 130
PERCENTIL = (1.0, 99.5)

GUARDAR_GIF = True
MOSTRAR     = True

CADA_N = 1

MASCARA_SOLIDO = True

DESCARTAR_CUADROS_INCOMPATIBLES = True

ESPECIES = ['C2H4', 'O2', 'CO2', 'H2O', 'CO', 'H2', 'H', 'O', 'OH']

COLUMNAS_ESPERADAS = {
    'temperatura': 3,
    'especies':    11,
    'campo':       4,
    'densidad':    3,
    'obstaculo':   3,
    'malla_nodos': 2,
    'superficie':  6,
    'tiempo':      1,
}

FAMILIAS_DE_CAMPO = ('temperatura', 'especies', 'campo',
                     'densidad', 'obstaculo', 'malla_nodos')


def info(msg):
    print(f"   {msg}")


def aviso(msg):
    print(f"   [AVISO] {msg}")


def fatal(msg):
    print()
    print("=" * 70)
    print(" ERROR")
    print("=" * 70)
    for linea in str(msg).split('\n'):
        print(f" {linea}")
    print("=" * 70)
    print(" Sugerencia: ejecuta  python3 ver_resultados.py diagnostico")
    print(" para ver exactamente que hay en esta carpeta.")
    print()
    raise SystemExit(1)


def cargar(nombre):
    if not os.path.exists(nombre):
        fatal(f"No existe el archivo '{nombre}'.\n"
              f"Carpeta actual: {os.getcwd()}")
    try:
        return np.loadtxt(nombre)
    except ValueError:
        aviso(f"reparando exponentes malformados en {nombre}")
        patron = re.compile(r'(\d)([+-]\d\d\d)')
        filas = []
        with open(nombre) as fh:
            for linea in fh:
                if not linea.strip() or linea.lstrip().startswith('#'):
                    continue
                filas.append([float(v) for v in
                              patron.sub(r'\1E\2', linea).split()])
        if not filas:
            fatal(f"El archivo '{nombre}' esta vacio o no tiene numeros.")
        return np.array(filas)


def numero_de(nombre, prefijo):
    m = re.search(rf'{prefijo}_(\d+)\.dat$', os.path.basename(nombre))
    return int(m.group(1)) if m else -1


def listar(prefijo):
    archivos = glob.glob(f'{prefijo}_*.dat')
    archivos = [a for a in archivos if numero_de(a, prefijo) >= 0]
    return sorted(archivos, key=lambda a: numero_de(a, prefijo))


def medir(nombre):
    n, ncols = 0, None
    with open(nombre, 'rb') as fh:
        for linea in fh:
            s = linea.strip()
            if not s or s.startswith(b'#'):
                continue
            n += 1
            if ncols is None:
                ncols = len(s.split())
    return n, (ncols or 0)


def inventario():
    inv = {}
    for fam in COLUMNAS_ESPERADAS:
        archivos = listar(fam)
        if not archivos:
            continue
        inv[fam] = [(numero_de(a, fam), a) + medir(a) for a in archivos]
    return inv


def diagnostico():
    print()
    print("=" * 70)
    print(" DIAGNOSTICO DE LA CARPETA")
    print("=" * 70)
    print(f" Directorio: {os.getcwd()}")

    todos_dat = sorted(glob.glob('*.dat'))
    print(f" Archivos .dat presentes: {len(todos_dat)}")
    if not todos_dat:
        print()
        print(" No hay NINGUN archivo .dat aqui.")
        print(" Posibles causas:")
        print("   - el solver todavia no escribio nada (escribe cada")
        print("     n_salida = 2000 pasos, puede tardar bastante);")
        print("   - estas en otra carpeta: muevete a donde ejecutaste el")
        print("     binario, no a donde esta el codigo fuente.")
        print("=" * 70)
        return

    inv = inventario()
    if not inv:
        print()
        print(" Hay archivos .dat, pero ninguno con el patron esperado")
        print(" (temperatura_XXXX.dat, especies_XXXX.dat, etc.).")
        print(" Primeros archivos encontrados:")
        for a in todos_dat[:10]:
            print(f"   {a}")
        print("=" * 70)
        return

    print()
    print(f" {'familia':<14} {'cuadros':>8} {'filas':>12} {'cols':>6}  estado")
    print(" " + "-" * 66)

    nodos_vistos = {}
    for fam, entradas in sorted(inv.items()):
        filas = {e[2] for e in entradas}
        cols  = {e[3] for e in entradas}
        idxs  = [e[0] for e in entradas]

        esperado = COLUMNAS_ESPERADAS[fam]
        if len(cols) > 1:
            estado = 'COLUMNAS INCONSISTENTES'
        elif cols and next(iter(cols)) != esperado:
            estado = f'esperaba {esperado} columnas'
        elif len(filas) > 1:
            estado = 'FILAS INCONSISTENTES'
        else:
            estado = 'ok'

        s_filas = (str(next(iter(filas))) if len(filas) == 1
                   else f'{min(filas)}..{max(filas)}')
        s_cols = (str(next(iter(cols))) if len(cols) == 1 else 'varias')
        print(f" {fam:<14} {len(entradas):>8} {s_filas:>12} {s_cols:>6}  {estado}")

        if fam in FAMILIAS_DE_CAMPO:
            for e in entradas:
                nodos_vistos.setdefault(e[2], []).append(e[1])

        reales = [i for i in idxs if i != 9999]
        if reales:
            faltan = sorted(set(range(min(reales), max(reales) + 1))
                            - set(reales))
            if faltan:
                aviso(f"{fam}: faltan los cuadros {faltan[:12]}"
                      f"{' ...' if len(faltan) > 12 else ''}")
        if 9999 in idxs:
            aviso(f"{fam}: existe el cuadro 9999 -> la simulacion DIVERGIO. "
                  "Ese archivo es el ultimo estado antes del fallo.")

    print()
    if len(nodos_vistos) == 1:
        n = next(iter(nodos_vistos))
        print(f" Todos los archivos de campo tienen {n} nodos. Consistente.")
    elif len(nodos_vistos) > 1:
        print(" *** PROBLEMA: hay archivos con DISTINTO numero de nodos ***")
        print()
        for n, archivos in sorted(nodos_vistos.items()):
            print(f"   {n:>7} nodos  ->  {len(archivos):>4} archivos, "
                  f"p.ej. {archivos[0]}")
        print()
        print(" Casi siempre significa que en la misma carpeta quedaron")
        print(" resultados de DOS corridas con mallas distintas:")
        print("     20000 nodos = 200 x 100  -> program.f90 (malla fija)")
        print("     33981 nodos = 241 x 141  -> combustion_pe.f90")
        print()
        print(" Solucion: borra los .dat de la corrida vieja, o separa cada")
        print(" corrida en su propia carpeta, y vuelve a ejecutar.")

    print("=" * 70)
    print()


def detectar_malla(verboso=True):
    global NX, NY, NYS, J_SUP, FORMATO

    if verboso:
        print(f"Carpeta de trabajo: {os.getcwd()}")

    inv = inventario()
    if not inv:
        fatal("No se encontro ningun archivo con el patron esperado\n"
              "(temperatura_XXXX.dat, especies_XXXX.dat, campo_XXXX.dat...).\n"
              f"Carpeta actual: {os.getcwd()}\n"
              f"Archivos .dat aqui: {len(glob.glob('*.dat'))}")

    nodos = {}
    for fam in FAMILIAS_DE_CAMPO:
        for e in inv.get(fam, []):
            nodos.setdefault(e[2], []).append(e[1])

    if not nodos:
        fatal("Hay archivos, pero ninguno de campo (temperatura, especies,\n"
              "campo, densidad, obstaculo). No se puede deducir la malla.")

    if len(nodos) > 1:
        aviso("hay archivos con distinto numero de nodos:")
        for n, a in sorted(nodos.items()):
            info(f"  {n} nodos -> {len(a)} archivos, p.ej. {a[0]}")
        aviso("se usara la malla MAYORITARIA y se descartara el resto.")

    n_nodos = max(nodos, key=lambda k: len(nodos[k]))

    fuente = None
    for fam in ('temperatura', 'especies', 'densidad',
                'obstaculo', 'malla_nodos'):
        for e in inv.get(fam, []):
            if e[2] == n_nodos:
                fuente = e[1]
                break
        if fuente:
            break
    if fuente is None and os.path.exists('malla_nodos.dat'):
        fuente = 'malla_nodos.dat'
    if fuente is None:
        fatal("No hay ningun archivo con columnas de coordenadas (x, y).\n"
              "campo_XXXX.dat no las trae; hace falta temperatura_,\n"
              "especies_, densidad_, obstaculo_ o malla_nodos.dat.")

    d = cargar(fuente)
    if d.ndim != 2 or d.shape[1] < 2:
        fatal(f"'{fuente}' no tiene el formato esperado: se leyo una tabla\n"
              f"de forma {d.shape} en vez de al menos 2 columnas.")

    n = d.shape[0]
    x = d[:, 0]

    nx = n
    for i in range(1, n):
        if x[i] <= x[i - 1]:
            nx = i
            break

    NX = NX_FORZADO if NX_FORZADO else nx
    NY = NY_FORZADO if NY_FORZADO else n // NX

    if NX <= 1 or NY <= 1 or NX * NY != n:
        fatal(f"No se pudo deducir la malla de '{fuente}'.\n"
              f"  nodos en el archivo : {n}\n"
              f"  NX deducido         : {NX}\n"
              f"  NY deducido         : {NY}\n"
              f"  NX * NY             : {NX*NY}\n"
              "Si conoces los valores, ponlos en NX_FORZADO y NY_FORZADO\n"
              "arriba del script.")

    FORMATO = 'movil' if listar('malla_nodos') else 'fija'

    if FORMATO == 'movil':
        if NYS_FORZADO:
            NYS = NYS_FORZADO
        else:
            obs = [e for e in inv.get('obstaculo', []) if e[2] == n_nodos]
            if obs:
                flag = cargar(obs[0][1])[:, 2].reshape(NY, NX)
                filas = np.where(flag.max(axis=1) > 0.5)[0]
                NYS = int(filas.max()) + 1 if filas.size else NY - 100
            else:
                NYS = NY - 100
                aviso("no hay obstaculo_XXXX.dat; se asume NYG = 101 "
                      f"-> NYS = {NYS}")
    else:
        NYS = 1

    J_SUP = NYS - 1

    if verboso:
        print(f"Malla detectada: {NX} x {NY} = {n} nodos   [{fuente}]")
        if FORMATO == 'movil':
            print("Formato: body-fitted (combustion_pe.f90), "
                  f"superficie en la fila j = {NYS}")
        else:
            print("Formato: malla cartesiana fija (program.f90), "
                  "el solido va como mascara")


def asegurar_malla():
    if NX == 0 or NY == 0:
        detectar_malla()


def remodelar(d, col, nombre):
    if d.ndim != 2:
        fatal(f"'{nombre}' no se leyo como tabla (forma {d.shape}).\n"
              "Puede estar truncado o tener una sola fila.")
    if col >= d.shape[1]:
        fatal(f"'{nombre}' tiene {d.shape[1]} columnas y se pidio la "
              f"numero {col+1}.\n"
              "Puede ser de otro solver o de otra version del codigo.")
    n = d.shape[0]
    if n != NX * NY:
        fatal(f"'{nombre}' tiene {n} nodos, pero la malla detectada es\n"
              f"  NX x NY = {NX} x {NY} = {NX*NY} nodos.\n\n"
              "Esto pasa casi siempre cuando en la carpeta quedaron\n"
              "resultados de DOS corridas con mallas distintas:\n"
              "  20000 nodos = 200 x 100  -> program.f90 (malla fija)\n"
              "  33981 nodos = 241 x 141  -> combustion_pe.f90\n\n"
              "Borra los .dat de la corrida vieja, o separa cada corrida\n"
              "en su propia carpeta.")
    return d[:, col].reshape(NY, NX)


def leer_malla(idx):
    nombre = f'malla_nodos_{idx:04d}.dat'
    if not os.path.exists(nombre):
        nombre = 'malla_nodos.dat'
    if not os.path.exists(nombre):
        for pref in ('temperatura', 'especies', 'densidad', 'obstaculo'):
            alt = f'{pref}_{idx:04d}.dat'
            if os.path.exists(alt):
                nombre = alt
                break
    d = cargar(nombre)
    return remodelar(d, 0, nombre), remodelar(d, 1, nombre)


def leer_tiempo(idx):
    nombre = f'tiempo_{idx:04d}.dat'
    if not os.path.exists(nombre):
        return None
    try:
        return float(np.loadtxt(nombre))
    except Exception:
        return None


def leer_mascara(idx):
    nombre = f'obstaculo_{idx:04d}.dat'
    if not os.path.exists(nombre):
        return None
    return remodelar(cargar(nombre), 2, nombre) > 0.5


def leer_campo(idx, campo):
    if campo == 'temperatura':
        nombre = f'temperatura_{idx:04d}.dat'
        return remodelar(cargar(nombre), 2, nombre)

    if campo in ESPECIES:
        nombre = f'especies_{idx:04d}.dat'
        return remodelar(cargar(nombre), 2 + ESPECIES.index(campo), nombre)

    if campo in ('velocidad', 'u', 'v', 'presion'):
        nombre = f'campo_{idx:04d}.dat'
        col = {'u': 0, 'v': 1, 'presion': 2, 'velocidad': 3}[campo]
        return remodelar(cargar(nombre), col, nombre)

    if campo == 'densidad':
        nombre = f'densidad_{idx:04d}.dat'
        return remodelar(cargar(nombre), 2, nombre)

    fatal(f"Campo desconocido: '{campo}'.\n"
          "Opciones: temperatura, velocidad, u, v, presion, densidad,\n"
          f"          {', '.join(ESPECIES)}\n"
          "Modos: diagnostico, superficie, historia")


def archivo_de(campo):
    if campo == 'temperatura':
        return 'temperatura'
    if campo in ESPECIES:
        return 'especies'
    if campo in ('velocidad', 'u', 'v', 'presion'):
        return 'campo'
    if campo == 'densidad':
        return 'densidad'
    fatal(f"Campo desconocido: '{campo}'.\n"
          "Opciones: temperatura, velocidad, u, v, presion, densidad,\n"
          f"          {', '.join(ESPECIES)}\n"
          "Modos: diagnostico, superficie, historia")


def etiqueta_de(campo):
    if campo == 'temperatura':
        return 'Temperatura (K)', PALETA_T
    if campo in ESPECIES:
        return f'Fraccion masica de {campo}', PALETA_Y
    if campo == 'velocidad':
        return 'Modulo de la velocidad (m/s)', PALETA_VEL
    if campo == 'u':
        return 'u (m/s)', 'coolwarm'
    if campo == 'v':
        return 'v (m/s)', 'coolwarm'
    if campo == 'presion':
        return 'Presion dinamica (Pa)', 'coolwarm'
    if campo == 'densidad':
        return 'Densidad (kg/m3)', 'cividis'
    return campo, PALETA_T


def animar(campo='temperatura'):
    asegurar_malla()

    prefijo = archivo_de(campo)
    print(f"Buscando {prefijo}_*.dat ...")
    archivos = listar(prefijo)

    if not archivos:
        otras = [f for f in COLUMNAS_ESPERADAS if listar(f)]
        fatal(f"No hay archivos {prefijo}_*.dat en esta carpeta.\n"
              f"Carpeta actual: {os.getcwd()}\n"
              + (f"Si hay, en cambio: {', '.join(otras)}"
                 if otras else "Tampoco hay ninguna otra familia de datos."))

    indices = [numero_de(a, prefijo) for a in archivos]
    info(f"encontrados {len(indices)} cuadros "
         f"(del {min(indices)} al {max(indices)})")

    if 9999 in indices:
        aviso("existe el cuadro 9999: la simulacion DIVERGIO. "
              "Ese es el ultimo estado antes del fallo.")

    if DESCARTAR_CUADROS_INCOMPATIBLES:
        buenos, malos = [], []
        for a, idx in zip(archivos, indices):
            n, _ = medir(a)
            (buenos if n == NX * NY else malos).append((idx, a, n))
        if malos:
            aviso(f"se descartan {len(malos)} cuadros con otro numero de "
                  f"nodos (la malla detectada tiene {NX*NY}):")
            for idx, a, n in malos[:8]:
                info(f"  {a}  ->  {n} nodos")
            if len(malos) > 8:
                info(f"  ... y {len(malos)-8} mas")
        indices = [b[0] for b in buenos]
        if not indices:
            fatal("Despues de descartar, no queda ningun cuadro utilizable.\n"
                  "Ejecuta el modo 'diagnostico' para ver el detalle.")

    if CADA_N > 1:
        sel = indices[::CADA_N]
        if sel[-1] != indices[-1]:
            sel.append(indices[-1])
        info(f"CADA_N = {CADA_N}: se usan {len(sel)} de {len(indices)} cuadros")
        indices = sel
    print(f"Leyendo {len(indices)} cuadros en memoria...")
    cuadros = []
    hay_mascara = False
    for k, idx in enumerate(indices, 1):
        if k % 20 == 0 or k == len(indices):
            info(f"  {k}/{len(indices)}")
        X, Y = leer_malla(idx)
        C = leer_campo(idx, campo)
        t = leer_tiempo(idx)

        if FORMATO == 'movil':
            if REGION == 'gas':
                X, Y, C = X[J_SUP:, :], Y[J_SUP:, :], C[J_SUP:, :]
        elif MASCARA_SOLIDO:
            solido = leer_mascara(idx)
            if solido is not None:
                C = np.ma.masked_array(C, mask=solido)
                hay_mascara = True

        cuadros.append((np.asarray(X, np.float32),
                        np.asarray(Y, np.float32), C, t, idx))

    if FORMATO == 'fija' and MASCARA_SOLIDO and not hay_mascara:
        aviso("no se encontraron obstaculo_XXXX.dat: el bloque solido se "
              "dibujara junto al gas y puede dominar la escala de color.")

    todos = np.ma.concatenate([np.ma.ravel(c[2]) for c in cuadros])
    todos = np.ma.compressed(todos)
    todos = todos[np.isfinite(todos)]
    if todos.size == 0:
        fatal("Todos los valores del campo son invalidos (NaN o infinito).\n"
              "Probablemente la simulacion divergio antes de escribir datos\n"
              "utiles. Revisa la salida del solver y el cuadro 9999.")
    v_min = float(np.percentile(todos, PERCENTIL[0]))
    v_max = float(np.percentile(todos, PERCENTIL[1]))
    if v_max - v_min < 1e-12:
        aviso(f"el campo '{campo}' es practicamente constante "
              f"(valor {v_min:.4g}); la animacion no mostrara variacion.")
        v_max = v_min + 1.0
    del todos
    info(f"rango de color: {v_min:.4g} a {v_max:.4g}")

    x_lo = min(float(c[0].min()) for c in cuadros)
    x_hi = max(float(c[0].max()) for c in cuadros)
    y_lo = min(float(c[1].min()) for c in cuadros)
    if FORMATO == 'movil' and REGION == 'gas':
        y_lo = 0.0
    y_hi = max(float(c[1].max()) for c in cuadros)
    info(f"dominio: x = [{x_lo*1e3:.1f}, {x_hi*1e3:.1f}] mm, "
         f"y = [{y_lo*1e3:.2f}, {y_hi*1e3:.2f}] mm")

    etiqueta, paleta = etiqueta_de(campo)
    norma = matplotlib.colors.Normalize(vmin=v_min, vmax=v_max)
    niveles = np.linspace(v_min, v_max, NIVELES)

    fig, ax = plt.subplots(figsize=(11, 4.2), dpi=DPI)

    sm = plt.cm.ScalarMappable(norm=norma, cmap=paleta)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label=etiqueta, pad=0.02)

    def dibujar(n):
        X, Y, C, t, idx = cuadros[n]
        ax.clear()

        if METODO == 'contour':
            ax.contourf(X * 1e3, Y * 1e3, C, levels=niveles,
                        cmap=paleta, extend='both')
        else:
            sombra = 'nearest' if hay_mascara else 'gouraud'
            ax.pcolormesh(X * 1e3, Y * 1e3, C, cmap=paleta,
                          norm=norma, shading=sombra)

        if FORMATO == 'movil':
            js = 0 if REGION == 'gas' else J_SUP
            ax.fill_between(X[js, :] * 1e3, y_lo * 1e3, Y[js, :] * 1e3,
                            color='0.30', zorder=4, edgecolor='none')
            ax.plot(X[js, :] * 1e3, Y[js, :] * 1e3,
                    color='k', lw=1.4, zorder=5)
            ax.text(0.5 * (x_lo + x_hi) * 1e3,
                    0.5 * (y_lo * 1e3 + Y[js, :].mean() * 1e3),
                    'COMBUSTIBLE (HDPE)', color='white', zorder=6,
                    ha='center', va='center', fontsize=10, alpha=0.75)

        titulo = f'{etiqueta}  -  cuadro {idx}'
        if t is not None:
            titulo += f'   t = {t:.4f} s'
        ax.set_title(titulo)
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')
        ax.set_xlim(x_lo * 1e3, x_hi * 1e3)
        ax.set_ylim(y_lo * 1e3, y_hi * 1e3)
        ax.set_aspect(ASPECTO)
        return []

    ani = animation.FuncAnimation(fig, dibujar, frames=len(cuadros),
                                  interval=1000 / FPS, blit=False)

    if GUARDAR_GIF:
        nombre_gif = f'animacion_{campo}.gif'
        print(f"Guardando '{nombre_gif}' ...")
        try:
            ani.save(nombre_gif, writer='pillow', fps=FPS, dpi=DPI)
            print(f"Listo: {os.path.abspath(nombre_gif)}")
        except Exception as e:
            aviso(f"no se pudo guardar el GIF: {e}")
            aviso("falta el paquete 'pillow'?  ->  pip install pillow")

    if MOSTRAR:
        plt.show()
    else:
        plt.close(fig)


def graficar_superficie():
    asegurar_malla()
    print("Buscando superficie_*.dat ...")
    archivos = listar('superficie')
    if not archivos:
        fatal("No hay archivos superficie_*.dat en esta carpeta.\n"
              "Ese diagnostico lo escribe solo combustion_pe.f90\n"
              "(malla ajustada al cuerpo); program.f90 no lo genera.\n"
              f"Carpeta actual: {os.getcwd()}")

    idx = numero_de(archivos[-1], 'superficie')
    info(f"usando el ultimo cuadro: {archivos[-1]}")
    d = cargar(archivos[-1])
    if d.ndim != 2 or d.shape[1] < 6:
        fatal(f"'{archivos[-1]}' deberia tener 6 columnas\n"
              "(x, espesor, Tw, rdot, q_conv, q_rad).")

    x, esp, Tw, rdot, qc, qr = (d[:, k] for k in range(6))
    t = leer_tiempo(idx)

    fig, ejes = plt.subplots(2, 2, figsize=(11, 7), dpi=DPI)
    titulo = f'Diagnostico de la superficie  -  cuadro {idx}'
    if t is not None:
        titulo += f'   t = {t:.4f} s'
    fig.suptitle(titulo)

    ejes[0, 0].plot(x * 1e3, esp * 1e3, 'k-')
    ejes[0, 0].set_ylabel('espesor del bloque (mm)')
    ejes[0, 0].set_title('Regresion')

    ejes[0, 1].plot(x * 1e3, Tw, 'r-')
    ejes[0, 1].set_ylabel('$T_w$ (K)')
    ejes[0, 1].set_title('Temperatura de pared')

    ejes[1, 0].plot(x * 1e3, rdot, 'b-')
    ejes[1, 0].set_ylabel(r'$\dot r$ (mm/s)')
    ejes[1, 0].set_title('Velocidad de regresion')

    ejes[1, 1].plot(x * 1e3, qc * 1e-6, label='conveccion')
    ejes[1, 1].plot(x * 1e3, qr * 1e-6, label='radiacion')
    ejes[1, 1].set_ylabel('flujo de calor (MW/m$^2$)')
    ejes[1, 1].set_title('Flujos a la pared')
    ejes[1, 1].legend()

    for e in ejes.ravel():
        e.set_xlabel('x (mm)')
        e.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig('superficie.png', dpi=DPI)
    print(f"Guardado: {os.path.abspath('superficie.png')}")
    if MOSTRAR:
        plt.show()
    else:
        plt.close(fig)


def graficar_historia():
    asegurar_malla()
    print("Buscando superficie_*.dat ...")
    archivos = listar('superficie')
    if not archivos:
        fatal("No hay archivos superficie_*.dat en esta carpeta.\n"
              "La historia temporal necesita ese diagnostico, que solo\n"
              "escribe combustion_pe.f90.\n"
              f"Carpeta actual: {os.getcwd()}")

    info(f"encontrados {len(archivos)} cuadros")
    tiempos, Tw_max, r_max, esp_med, Tg_max = [], [], [], [], []

    for a in archivos:
        idx = numero_de(a, 'superficie')
        d = cargar(a)
        t = leer_tiempo(idx)
        tiempos.append(t if t is not None else idx)
        esp_med.append(d[:, 1].mean() * 1e3)
        Tw_max.append(d[:, 2].max())
        r_max.append(d[:, 3].max())

        nomT = f'temperatura_{idx:04d}.dat'
        if os.path.exists(nomT) and medir(nomT)[0] == NX * NY:
            T = remodelar(cargar(nomT), 2, nomT)
            Tg_max.append(T[J_SUP:, :].max())
        else:
            Tg_max.append(np.nan)

    if all(np.isnan(v) for v in Tg_max):
        aviso("no se pudo leer ninguna temperatura del gas; se grafica "
              "solo lo que viene de superficie_*.dat")

    fig, ejes = plt.subplots(3, 1, figsize=(9, 8), dpi=DPI, sharex=True)

    ejes[0].plot(tiempos, Tg_max, 'r-', label='$T$ max del gas')
    ejes[0].plot(tiempos, Tw_max, 'k--', label='$T_w$ max')
    ejes[0].set_ylabel('T (K)')
    ejes[0].legend()

    ejes[1].plot(tiempos, r_max, 'b-')
    ejes[1].set_ylabel(r'$\dot r$ max (mm/s)')

    ejes[2].plot(tiempos, esp_med, 'g-')
    ejes[2].set_ylabel('espesor medio (mm)')
    ejes[2].set_xlabel('tiempo (s)')

    for e in ejes:
        e.grid(alpha=0.3)

    fig.suptitle('Historia temporal')
    fig.tight_layout()
    fig.savefig('historia.png', dpi=DPI)
    print(f"Guardado: {os.path.abspath('historia.png')}")
    if MOSTRAR:
        plt.show()
    else:
        plt.close(fig)


def main(arg='temperatura'):
    if arg in ('diagnostico', 'diag', '-d'):
        diagnostico()
    elif arg in ('ayuda', '-h', '--help'):
        print(__doc__)
    elif arg == 'superficie':
        graficar_superficie()
    elif arg == 'historia':
        graficar_historia()
    else:
        animar(arg)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'temperatura')
