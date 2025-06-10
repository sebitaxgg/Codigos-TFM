# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 15:40:36 2025

@author: Sebas
"""
plantilla_archivo = "ETH_Rotacion_{}.txt"

ruta_directorio = r"C:\Users\Sebas\Desktop\Matlab_TFG\ETH_Scaling"

import numpy as np
import os
import glob
import re
from scipy.stats import linregress
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def modelo_1(x, A, B):
    return A / np.sqrt(x) + B

# Función 2: A/x^c + B
def modelo_2(x, A, c):
    return A / (x ** c)
# Define tus modelos si no lo has hecho antes:
# modelo_1 = lambda x, A, B: A / np.sqrt(x) + B
# modelo_2 = lambda x, A, c, B: A / (x**c) + B

energia_min = -3.1
energias_max = [-2.3]  # Puedes poner solo un valor o más
orden_deseado = [20, 28, 40, 56, 80, 112, 160, 224, 320]

# Configurar layout de subplots
num_valores = len(energias_max)
fig, axs = plt.subplots(2, num_valores, figsize=(6 * num_valores, 10))

if num_valores == 1:
    axs = np.array([[axs[0]], [axs[1]]])  # Asegura que axs[i][0] funcione si hay solo una columna

for idx, energia_max in enumerate(energias_max):
    desviacion_rot = []
    maximo_rot = []

    for numero in orden_deseado:
        nombre_archivo = plantilla_archivo.format(numero)
        ruta_archivo = os.path.join(ruta_directorio, nombre_archivo)
        if not os.path.exists(ruta_archivo):
            print(f"[{nombre_archivo}] → No encontrado.")
            continue

        try:
            datos = np.loadtxt(ruta_archivo)
            energias = datos[:, 0]
            auto_re_1 = datos[:, 8]
            auto_re_2 = datos[:, 9]
            auto_re_3 = datos[:, 10]

            filtro = (energias >= energia_min) & (energias <= energia_max)
            magnitudes_filtradas = np.concatenate([auto_re_1[filtro], auto_re_2[filtro], auto_re_3[filtro]])

            if len(magnitudes_filtradas) == 0:
                print(f"[{nombre_archivo}] → No hay datos en el rango especificado.")
            else:
                media = np.mean(magnitudes_filtradas)
                desviacion = np.std(magnitudes_filtradas)
                max_abs_rot = np.max(np.abs(magnitudes_filtradas))
                desviacion_rot.append(desviacion)
                maximo_rot.append(max_abs_rot)
                print(f"[energia_max = {energia_max}] [N = {numero}] → σ: {desviacion:.5f}, Máx: {max_abs_rot:.5f}")
        except Exception as e:
            print(f"[{nombre_archivo}] → Error: {e}")

    x_vals = np.linspace(min(orden_deseado), max(orden_deseado), 1000)

    # === σ vs N ===
    ax_sigma = axs[0, idx]
    ax_sigma.scatter(orden_deseado, desviacion_rot, color='blue', label='σ')

    popt1_sigma, _ = curve_fit(modelo_1, orden_deseado, desviacion_rot, maxfev=10000)
    A1, B1 = popt1_sigma
    popt2_sigma, _ = curve_fit(modelo_2, orden_deseado, desviacion_rot, maxfev=10000)
    A2, c2 = popt2_sigma

    ax_sigma.plot(x_vals, modelo_1(x_vals, *popt1_sigma), color='red',
                  label=f'A/sqrt(x)+B\nA={A1:.3f}, B={B1:.3f}')
    ax_sigma.plot(x_vals, modelo_2(x_vals, *popt2_sigma), color='green',
                  label=f'A/x^c + B\nA={A2:.3f}, c={c2:.3f}')
    ax_sigma.set_title(f"Autovalores Rotacion σ vs N\n(energia_max = {energia_max})")
    ax_sigma.set_xlabel('N')
    ax_sigma.set_ylabel('σ')
    ax_sigma.grid(True)
    ax_sigma.legend()

    # === Máximo vs N ===
    ax_max = axs[1, idx]
    ax_max.scatter(orden_deseado, maximo_rot, color='blue', label='Máximo')

    popt1_max, _ = curve_fit(modelo_1, orden_deseado, maximo_rot, maxfev=10000)
    A1, B1 = popt1_max
    popt2_max, _ = curve_fit(modelo_2, orden_deseado, maximo_rot, maxfev=10000)
    A2, c2 = popt2_max

    ax_max.plot(x_vals, modelo_1(x_vals, *popt1_max), color='red',
                label=f'A/sqrt(x)+B\nA={A1:.3f}, B={B1:.3f}')
    ax_max.plot(x_vals, modelo_2(x_vals, *popt2_max), color='green',
                label=f'A/x^c + B\nA={A2:.3f}, c={c2:.3f}')
    ax_max.set_title(f"Máximo Rotacion vs N\n(energia_max = {energia_max})")
    ax_max.set_xlabel('N')
    ax_max.set_ylabel('Máximo')
    ax_max.grid(True)
    ax_max.legend()

plt.tight_layout()
plt.show()
