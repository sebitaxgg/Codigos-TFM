# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 15:40:39 2025

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
def modelo_2(x, A, c, B):
    return A / (x ** c) + B
plantilla_archivo = "ETH_Hopping_Autoval_Full_{}.dat"
energia_min = -3.1
energia_max = -2.3

# Rutas a los archivos
ruta_observado = r"C:\Users\Sebas\Desktop\Matlab_TFG\ETH_Scaling"
ruta_referencia = r"C:\Users\Sebas\Desktop\Matlab_TFG\ETH_Scaling\micro_hopping.dat"

energias_max = [-2.3, -2.4, -2.5]  # Puedes poner uno solo o más
orden_deseado = [20, 28, 40, 56, 80, 112, 160, 224, 320]

ref = np.loadtxt(ruta_referencia)
energia_ref = ref[:, 0]
magnitud_real = ref[:, 2]  # Tercera columna, magnitud real

# Layout de subgráficos
num_valores = len(energias_max)
fig, axs = plt.subplots(2, num_valores, figsize=(6 * num_valores, 10))
if num_valores == 1:
    axs = np.array([[axs[0]], [axs[1]]])  # Para el caso de una sola columna

for idx, energia_max in enumerate(energias_max):
    sigma_traza_Hopp = []
    maximo_traza_hopp = []

    for numero in orden_deseado:
        nombre_archivo = plantilla_archivo.format(numero)
        ruta_archivo = os.path.join(ruta_observado, nombre_archivo)

        if not os.path.exists(ruta_archivo):
            print(f"[{nombre_archivo}] → No encontrado.")
            continue

        datos_obs = np.loadtxt(ruta_archivo)
        energia_obs = datos_obs[:, 0]
        magnitud_obs = datos_obs[:, 1]  # Columna observada

        filtro = (energia_obs >= energia_min) & (energia_obs <= energia_max)
        energia_filtradas = energia_obs[filtro]
        magnitudes_filtradas = magnitud_obs[filtro]

        if len(energia_filtradas) == 0:
            print(f"[{nombre_archivo}] → Sin datos en el rango.")
            continue

        # Interpolación
        magnitud_interp = np.interp(energia_filtradas, energia_ref, magnitud_real)
        residuos = magnitud_interp - magnitudes_filtradas 

        sigma = np.std(residuos)
        maximo = np.max(np.abs(residuos))
        sigma_traza_Hopp.append(sigma)
        maximo_traza_hopp.append(maximo)

        print(f"[energia_max = {energia_max}] N = {numero} → σ: {sigma:.5f}  Max: {maximo:.5f}")

    x_vals = np.linspace(min(orden_deseado), max(orden_deseado), 1000)

    # ========== Sigma ==========
    ax_sigma = axs[0, idx]
    ax_sigma.scatter(orden_deseado, sigma_traza_Hopp, color='blue', label='σ')

    popt1, _ = curve_fit(modelo_1, orden_deseado, sigma_traza_Hopp, maxfev=10000)
    A1, B1 = popt1
    popt2, _ = curve_fit(modelo_2, orden_deseado, sigma_traza_Hopp, maxfev=10000)
    A2, c2, B2 = popt2

    ax_sigma.plot(x_vals, modelo_1(x_vals, *popt1), 'r-', label=f'A/sqrt(x)+B\nA={A1:.3f}, B={B1:.3f}')
    ax_sigma.plot(x_vals, modelo_2(x_vals, *popt2), 'g-', label=f'A/x^c + B\nA={A2:.3f}, c={c2:.3f}, B={B2:.3f}')
    ax_sigma.set_title(f"Autovalores Hopping σ \n(energia_max = {energia_max})")
    ax_sigma.set_xlabel('N')
    ax_sigma.set_ylabel('σ')
    ax_sigma.grid(True)
    ax_sigma.legend()

    # ========== Máximo ==========
    ax_max = axs[1, idx]
    ax_max.scatter(orden_deseado, maximo_traza_hopp, color='blue', label='Máximo')

    popt1, _ = curve_fit(modelo_1, orden_deseado, maximo_traza_hopp, maxfev=10000)
    A1, B1 = popt1
    popt2, _ = curve_fit(modelo_2, orden_deseado, maximo_traza_hopp, maxfev=10000)
    A2, c2, B2 = popt2

    ax_max.plot(x_vals, modelo_1(x_vals, *popt1), 'r-', label=f'A/sqrt(x)+B\nA={A1:.3f}, B={B1:.3f}')
    ax_max.plot(x_vals, modelo_2(x_vals, *popt2), 'g-', label=f'A/x^c + B\nA={A2:.3f}, c={c2:.3f}, B={B2:.3f}')
    ax_max.set_title(f"Máximo Autovalores Hopping \n(energia_max = {energia_max})")
    ax_max.set_xlabel('N')
    ax_max.set_ylabel('Máximo')
    ax_max.grid(True)
    ax_max.legend()

plt.tight_layout()
plt.show()