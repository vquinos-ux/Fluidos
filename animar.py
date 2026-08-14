#!/usr/bin/env python3
"""
Visualiza el campo de velocidad en el tiempo a partir de los archivos
generados sobre una malla no estructurada (triángulos).

Uso:
    python animar_triangulos.py              # anima en pantalla
    python animar_triangulos.py --save mp4   # guarda campo.mp4 (requiere ffmpeg)
    python animar_triangulos.py --save gif   # guarda campo.gif (requiere pillow)
"""

import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.animation import FuncAnimation

# --- 1. Cargar la Geometría Estática ---
try:
    nodos = np.genfromtxt("malla_nodos.dat")
    x = nodos[:, 0]
    y = nodos[:, 1]
    triangulos = np.genfromtxt("malla_triangulos.dat", dtype=int)
    
    # Creamos el objeto de triangulación que Matplotlib usará de esqueleto
    malla = mtri.Triangulation(x, y, triangulos)
except OSError:
    sys.exit("Error: No se encontraron 'malla_nodos.dat' o 'malla_triangulos.dat'.")

# --- 2. Cargar todos los frames ordenados ---
archivos = sorted(glob.glob("campo_*.dat"))
if not archivos:
    sys.exit("No se encontraron archivos campo_*.dat.")

frames = []
for nombre in archivos:
    datos = np.genfromtxt(nombre)
    u = datos[:, 0]
    v = datos[:, 1]
    # presion = datos[:, 2]  # No la usamos en la animación por ahora
    speed = datos[:, 3]
    frames.append((u, v, speed))

print(f"Cargados {len(frames)} frames de {archivos[0]} a {archivos[-1]}")

u0, v0, s0 = frames[0]
vmax = max(f[2].max() for f in frames)

# --- 3. Configurar la figura ---
fig, ax = plt.subplots(figsize=(11, 5))

# tripcolor es la versión de pcolormesh pero para triángulos. 
# shading='gouraud' hace un degradado suave entre los nodos.
img = ax.tripcolor(malla, s0, shading='gouraud', cmap="turbo", vmin=0, vmax=vmax)
cb = fig.colorbar(img, ax=ax, label="Rapidez |V|  (m/s)")

# Vectores de velocidad
# Como x e y son listas 1D, un "paso" simplemente salta nodos en la lista.
paso = 5  # Cambia esto si las flechas se ven muy amontonadas o muy separadas
q = ax.quiver(x[::paso], y[::paso],
              u0[::paso], v0[::paso],
              color="white", scale=40, width=0.0015)

ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
titulo = ax.set_title("Campo de velocidad (Malla No Estructurada) - frame 1")

# --- 4. Función de animación ---
def actualizar(k):
    u, v, speed = frames[k]
    img.set_array(speed)  # Actualiza los colores en los nodos
    q.set_UVC(u[::paso], v[::paso])  # Actualiza las flechas
    titulo.set_text(f"Campo de velocidad - frame {k + 1}/{len(frames)}")
    return img, q, titulo

anim = FuncAnimation(fig, actualizar, frames=len(frames),
                     interval=150, blit=False)

# --- 5. Mostrar o guardar ---
if len(sys.argv) > 2 and sys.argv[1] == "--save":
    formato = sys.argv[2]
    if formato == "mp4":
        anim.save("campo_triangulos.mp4", writer="ffmpeg", fps=8, dpi=120)
        print("Guardado: campo_triangulos.mp4")
    elif formato == "gif":
        anim.save("campo_triangulos.gif", writer="pillow", fps=8)
        print("Guardado: campo_triangulos.gif")
    else:
        sys.exit("Formato no soportado. Usa 'mp4' o 'gif'.")
else:
    plt.tight_layout()
    plt.show()