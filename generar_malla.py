#!/usr/bin/env python3
"""
(malla pasada no entra en el programa)
Genera la malla no estructurada (triángulos) y la guarda para Fortran y Python.
Salida: malla_nodos.dat, malla_triangulos.dat y malla.png
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

# --- 1. Parametros (iguales a tu simulation.f90 original) ---
NX, NY = 200, 100
DX, DY = 0.1, 0.1
Lx, Ly = (NX - 1) * DX, (NY - 1) * DY   # tamano fisico del dominio

# --- 2. Crear las coordenadas de los nodos ---
xn = np.linspace(0, Lx, NX)
yn = np.linspace(0, Ly, NY)
X, Y = np.meshgrid(xn, yn)

# Aplanamos las matrices a listas 1D (como le gusta a la malla no estructurada)
x_nodos = X.flatten()
y_nodos = Y.flatten()

# --- 3. Generar la Triangulación (Algoritmo Delaunay) ---
malla = mtri.Triangulation(x_nodos, y_nodos)
triangulos = malla.triangles # Esto es la matriz de conectividad

# --- 4. EXPORTAR LOS ARCHIVOS VITALES ---
# Guardamos las coordenadas para que Fortran y los animadores las lean
np.savetxt("malla_nodos.dat", np.column_stack((x_nodos, y_nodos)), fmt="%14.6e")
# Guardamos cómo se conectan los triángulos
np.savetxt("malla_triangulos.dat", triangulos, fmt="%d")

# --- 5. Graficar la malla para visualizarla ---
fig, ax = plt.subplots(figsize=(13, 5.5))

# ax.triplot dibuja el esqueleto de triángulos
ax.triplot(malla, color="steelblue", lw=0.4, alpha=0.7, zorder=2)

# --- Paredes solidas (no deslizamiento) ---
ax.plot([0, Lx], [0, 0], color="black", lw=3, zorder=4)
ax.plot([0, Lx], [Ly, Ly], color="black", lw=3, zorder=4)

# --- Etiquetas de las condiciones de contorno ---
bbox = dict(boxstyle="round", fc="white", ec="0.5")
ax.annotate("ENTRADA\n u=U_inf, v=0", (0, Ly / 2), (-1.9, Ly / 2),
            ha="center", va="center", fontsize=9, color="green",
            bbox=bbox, arrowprops=dict(arrowstyle="->", color="green"))
ax.annotate("SALIDA\n gradiente nulo", (Lx, Ly / 2), (Lx + 1.9, Ly / 2),
            ha="center", va="center", fontsize=9, color="purple",
            bbox=bbox, arrowprops=dict(arrowstyle="->", color="purple"))
ax.text(Lx / 2, Ly + 0.35, "PARED SUPERIOR  (no deslizamiento: u=v=0)",
        ha="center", fontsize=9, color="dimgray")
ax.text(Lx / 2, -0.55, "PARED INFERIOR  (no deslizamiento: u=v=0)",
        ha="center", fontsize=9, color="dimgray")

# --- Estetica ---
ax.set_xlim(-3.5, Lx + 3.5)
ax.set_ylim(-1.0, Ly + 0.9)
ax.set_aspect("equal")
ax.set_xlabel("x  (m)")
ax.set_ylabel("y  (m)")

n_nodos = len(x_nodos)
n_tri = len(triangulos)
ax.set_title(f"Malla Computacional Triangulada: {n_nodos} nodos, {n_tri} triángulos")

plt.tight_layout()
plt.savefig("malla.png", dpi=130)
print(f"ÉXITO: Se crearon malla_nodos.dat, malla_triangulos.dat y malla.png")
