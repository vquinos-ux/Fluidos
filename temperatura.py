import numpy as np
import matplotlib.pyplot as plt
import glob, re, os

# Parámetros de la malla (deben coincidir con Nx y Ny del código Fortran)
Nx = 200 #cambiados
Ny = 100 #cambiados

def leer_archivo(nombre_archivo):
    """Lee un archivo .dat con columnas x y T y devuelve array 2D de temperatura"""
    datos = np.loadtxt(nombre_archivo)   # columnas: x, y, T
    T = datos[:, 2].reshape(Ny, Nx)       # asumiendo orden lexicográfico (x rápido)
    return T

# Obtener lista de archivos ordenados por número de frame
archivos = glob.glob('temperatura_*.dat')
if not archivos:
    print("No se encontraron archivos temperatura_*.dat en la carpeta actual.")
    exit()

# Función para extraer el número de frame del nombre del archivo
def extraer_numero(nombre):
    m = re.search(r'temperatura_(\d+)\.dat', nombre)
    return int(m.group(1)) if m else 0

archivos_ordenados = sorted(archivos, key=extraer_numero)

# Cargar el primer archivo para obtener las coordenadas (son constantes para todos)
datos_primer = np.loadtxt(archivos_ordenados[0])
x_vec = datos_primer[:, 0].reshape(Ny, Nx)[0, :]   # vector x (1D)
y_vec = datos_primer[:, 1].reshape(Ny, Nx)[:, 0]   # vector y (1D)
T_primer = datos_primer[:, 2].reshape(Ny, Nx)

print(f"Se encontraron {len(archivos_ordenados)} archivos de temperatura.")

# Configurar la figura para la animación
plt.ion()   # modo interactivo
fig, ax = plt.subplots(figsize=(10, 4))
im = ax.contourf(x_vec, y_vec, T_primer, levels=50, cmap='hot')
cbar = plt.colorbar(im, ax=ax, label='Temperatura (K)')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title('Temperatura – Fotograma 0')
ax.axis('equal')

# Mostrar el primer fotograma
plt.show()

# Recorrer todos los fotogramas (incluyendo el primero para actualizar título)
for idx, archivo in enumerate(archivos_ordenados, start=1):
    T = leer_archivo(archivo)
    
    # EL ORDEN CORRECTO:
    # 1. Primero eliminamos la barra de color del frame anterior
    cbar.remove()
    
    # 2. Luego borramos el gráfico anterior
    ax.clear()   
    
    # 3. Dibujamos el nuevo contorno
    im = ax.contourf(x_vec, y_vec, T, levels=50, cmap='hot')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_title(f'Temperatura – Fotograma {idx} ({archivo})')
    ax.axis('equal')
    
    # 4. Creamos la nueva barra de color
    cbar = plt.colorbar(im, ax=ax, label='Temperatura (K)')
    
    plt.draw()
    plt.pause(0.05)   # pausa de 50 ms entre fotogramas

# Al finalizar la animación, mantener la ventana abierta hasta que el usuario la cierre
print("Animación completada. Cierra la ventana para salir.")
plt.ioff()
plt.show()