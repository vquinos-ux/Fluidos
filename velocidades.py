#!/usr/bin/env python3
"""
Grafica el perfil de velocidad u(y) en la mitad del dominio,
a partir de datos en malla no estructurada (triángulos).

Salida: perfil_velocidad.png
"""

import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

# --- 1. Cargar la Geometría (Nodos y Triángulos) ---
try:
    nodos = np.genfromtxt("malla_nodos.dat")
    x = nodos[:, 0]
    y = nodos[:, 1]
    
    # dtype=int es crucial porque los índices de los triángulos son números enteros
    triangulos = np.genfromtxt("malla_triangulos.dat", dtype=int)
except OSError:
    sys.exit("Error: No se encontraron 'malla_nodos.dat' o 'malla_triangulos.dat'.")

# --- 2. Cargar el último frame de simulación ---
archivos = sorted(glob.glob("campo_*.dat"))
if not archivos:
    sys.exit("No se encontraron archivos campo_*.dat. Ejecuta primero la simulación.")
nombre = archivos[-1]

datos = np.genfromtxt(nombre)
u = datos[:, 0]  # En nuestro nuevo Fortran, 'u' está en la primera columna (0)

# --- 3. Crear el Interpolador sobre la Malla ---
# Armamos el "esqueleto" para que Python entienda la geometría
malla = mtri.Triangulation(x, y, triangulos)
# Creamos una función que interpola la velocidad 'u' en cualquier punto (x,y)
interpolador = mtri.LinearTriInterpolator(malla, u)

# --- 4. Extraer el Perfil en la Mitad del Dominio ---
x_min, x_max = x.min(), x.max()
y_min, y_max = y.min(), y.max()
x_mid = (x_min + x_max) / 2.0

# Creamos 100 puntos artificiales formando una línea vertical en x_mid
y_linea = np.linspace(y_min, y_max, 100)
x_linea = np.full_like(y_linea, x_mid)

# Le pedimos al interpolador la velocidad exacta en esos 100 puntos
perfil = interpolador(x_linea, y_linea)

# --- 5. Gráfica ---
fig, ax = plt.subplots(figsize=(5, 6))
ax.plot(perfil, y_linea, "-", color="crimson", lw=2.5, label="Velocidad Interpolada")

# Líneas de las paredes inferior y superior
ax.axhline(y_min, color="black", lw=2)
ax.axhline(y_max, color="black", lw=2)
ax.axvline(0, color="0.7", lw=0.8, linestyle="--")

ax.set_xlabel("Velocidad  u  (m/s)")
ax.set_ylabel("y  (m)")
ax.set_title(f"Perfil de velocidad en x = {x_mid:.2f} m\n({nombre})")
ax.grid(True, alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig("perfil_velocidad.png", dpi=130)
print(f"Guardado: perfil_velocidad.png  (Interpolado en x = {x_mid:.2f} m)")