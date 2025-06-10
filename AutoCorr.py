# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 15:37:16 2025

@author: Sebas
"""
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

#Corriente
plantilla_archivo = "Corriente_Positivos_{}.dat"

ruta_directorio = r"C:\Users\Sebas\Desktop\Matlab_TFG\ETH_Scaling"

energia_min = -3.1
energias_max = [-2.3,-2.4,-2.5]  # Puedes poner los valores que quieras, incluso solo uno
orden_deseado = [20, 28, 40, 56, 80, 112, 160, 224, 320]

# Calcular diseño de subplots: siempre 2 filas, columnas según número de energia_max
num_valores = len(energias_max)
num_filas = 2
num_columnas = num_valores

fig, axs = plt.subplots(num_filas, num_columnas, figsize=(6 * num_columnas, 10))

# Si solo hay un valor, axs no es una matriz 2D, así que lo convertimos en una
if num_valores == 1:
    axs = np.array([[axs[0]], [axs[1]]])  # Asegura que axs[0, 0] y axs[1, 0] funcionen

for idx, energia_max in enumerate(energias_max):
    desviacion_corr = []
    maximo_corr = []

    for numero in orden_deseado:
        nombre_archivo = plantilla_archivo.format(numero)
        ruta_archivo = os.path.join(ruta_directorio, nombre_archivo)
        if not os.path.exists(ruta_archivo):
            print(f"[{nombre_archivo}] → No encontrado.")
            continue

        try:
            datos = np.loadtxt(ruta_archivo)
            energias = datos[:, 0]
            auto_pos = datos[:, 1]
            auto_neg = datos[:, 2]
            
            filtro = (energias >= energia_min) & (energias <= energia_max)
            magnitudes_filtradas = np.concatenate([auto_pos[filtro], auto_neg[filtro]])

            if len(magnitudes_filtradas) == 0:
                print(f"[{nombre_archivo}] → No hay datos en el rango especificado.")
            else:
                media = np.mean(magnitudes_filtradas)
                desviacion = np.std(magnitudes_filtradas)
                desviacion_corr.append(desviacion)
                max_abs = np.max(np.abs(magnitudes_filtradas))
                maximo_corr.append(max_abs)
                print(f"[energia_max = {energia_max}] [N = {numero}] → Desviación: {desviacion:.5f}  Máximo: {max_abs:.5f}")
        except Exception as e:
            print(f"[{nombre_archivo}] → Error al procesar: {e}")

    x_vals = np.linspace(min(orden_deseado), max(orden_deseado), 1000)

    # --- Gráfico de Desviaciones ---
    ax1 = axs[0, idx]
    ax1.scatter(orden_deseado, desviacion_corr, color='blue', label='Sigma')

    popt1_sigma, _ = curve_fit(modelo_1, orden_deseado, desviacion_corr, maxfev=10000)
    popt2_sigma, _ = curve_fit(modelo_2, orden_deseado, desviacion_corr, maxfev=10000)

    A1, B1 = popt1_sigma
    A2, c2 = popt2_sigma

    ax1.plot(x_vals, modelo_1(x_vals, *popt1_sigma), color='red',
             label=f'A/sqrt(x) + B\nA={A1:.3f}, B={B1:.3f}')
    ax1.plot(x_vals, modelo_2(x_vals, *popt2_sigma), color='green',
             label=f'A/x^c + B\nA={A2:.3f}, c={c2:5f}')

    ax1.set_title(f"Autovalores Corriente σ vs N\n(energia_max = {energia_max})")
    ax1.set_xlabel('N')
    ax1.set_ylabel('σ')
    ax1.legend()
    ax1.grid(True)

    # --- Gráfico de Máximos ---
    ax2 = axs[1, idx]
    ax2.scatter(orden_deseado, maximo_corr, color='blue', label='Máximo')
    
    
    popt1_max, _ = curve_fit(modelo_1, orden_deseado, maximo_corr, maxfev=10000)
    popt2_max, _ = curve_fit(modelo_2, orden_deseado, maximo_corr, maxfev=10000)
    
    A1, B1 = popt1_max
    A2, c2 = popt2_max

    ax2.plot(x_vals, modelo_1(x_vals, *popt1_max), color='red',
             label=f'A/sqrt(x) + B\nA={A1:.3f}, B={B1:.3f}')
    ax2.plot(x_vals, modelo_2(x_vals, *popt2_max), color='green',
             label=f'A/x^c + B\nA={A2:.3f}, c={c2:.5f}')

    ax2.set_title(f"Máximo Corriente vs N\n(energia_max = {energia_max}) ")
    ax2.set_xlabel('N')
    ax2.set_ylabel('Máximo')
    ax2.legend()
    ax2.grid(True)

plt.tight_layout()
plt.show()