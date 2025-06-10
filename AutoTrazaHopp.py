#Autovalores Hopping 2 (Autoval- Traza)
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

plantilla_archivo = "ETH_Hopping_{}.txt"
energia_min = -3.1
energia_max = [-2.3]
orden_deseado = [20, 28, 40, 56, 80, 112, 160, 224, 320]
ruta_observado = r"C:\Users\Sebas\Desktop\Matlab_TFG\ETH_Scaling"

num_valores = len(energia_max)
fig, axs = plt.subplots(2, num_valores, figsize=(6 * num_valores, 10))
if num_valores == 1:
    axs = np.array([[axs[0]], [axs[1]]])  # Para el caso de una sola columna


# Leer los datos de referencia (estos no cambian)
for idx, energia_max in enumerate(energia_max):
    sigma_Auto_2 = []
    max_Auto_2 = []
    
    for numero in orden_deseado:
        nombre_archivo = plantilla_archivo.format(numero)
        ruta_archivo = os.path.join(ruta_observado, nombre_archivo)
        archivo = ruta_archivo
        if not os.path.exists(ruta_archivo):
            print(f"[{nombre_archivo}] → No encontrado.")
            continue
        # Leer datos del archivo observado
        datos_obs = np.loadtxt(archivo)
        energia_obs = datos_obs[:, 0]  # Columna de energías observadas
        traza = datos_obs[:, 5]  # Columna de magnitudes observadas
        autoval_1 = datos_obs[:, 6] 
        autoval_2 = datos_obs[:, 7]  
        autoval_3 = datos_obs[:, 8]  
        filtro = (energia_obs >= energia_min) & (energia_obs <= energia_max)
        traza_filtrada = traza[filtro]
        autoval_1_filtrada = autoval_1[filtro]
        autoval_2_filtrada = autoval_2[filtro]
        autoval_3_filtrada = autoval_3[filtro]
        
        residuo_1 = traza_filtrada/3-autoval_1_filtrada
        residuo_2 = traza_filtrada/3-autoval_2_filtrada
        residuo_3 = traza_filtrada/3-autoval_3_filtrada
        
        residuo = np.concatenate([residuo_1,residuo_2,residuo_3])
    
        # Calcular la desviación estándar de los residuos
        sigma_residuos = np.std(residuo)
        sigma_Auto_2.append(sigma_residuos)
        max_abs = np.max(np.abs(residuo))
        max_Auto_2.append(max_abs)
    
        # Imprimir los resultados para cada archivo
        print(f" N = {numero} Sigma : {sigma_residuos:.5f} Maximo : {max_abs:.5f}")
    
    ax_sigma = axs[0, idx]
    ax_sigma.scatter(orden_deseado, sigma_Auto_2, color='blue', label='σ')

    popt_1, _ = curve_fit(modelo_1, orden_deseado, sigma_Auto_2,maxfev = 10000)    
    A_1, B_1 = popt_1  # Parámetros ajustados de la función 1
    
    # Ajustar a la función 2
    popt_2, _ = curve_fit(modelo_2, orden_deseado, sigma_Auto_2, maxfev = 10000)
    A_2, c_2, B_2 = popt_2  # Parámetros ajustados de la función 2
    
    
    # Graficar los ajustes
    x_vals = np.linspace(min(orden_deseado), max(orden_deseado), 1000)
    
    # Ajuste de la función 1
    ax_sigma.plot(x_vals, modelo_1(x_vals, *popt_1), color='red', label=f'A/sqrt(x) + B: A={A_1:.3f}, B={B_1:.3f}')
    
    # Ajuste de la función 2
    ax_sigma.plot(x_vals, modelo_2(x_vals, *popt_2), color='green', label=f'A/x^c + B: A={A_2:.3f}, c={c_2:.3f}, B={B_2:.3f}')
    
    # Etiquetas y título
    ax_sigma.set_xlabel('N')
    ax_sigma.set_ylabel('Sigma')
    ax_sigma.set_title('Auto Hopp - Traza/3 ')
    ax_sigma.legend()
    
    # Mostrar gráfico
    ax_sigma.grid(True)
    ax_sigma.legend()
    
    # ---------------
    
    ax_max = axs[1, idx]
    ax_max.scatter(orden_deseado, max_Auto_2, color='blue', label='Máximo')


    popt_1, _ = curve_fit(modelo_1, orden_deseado, max_Auto_2,maxfev = 10000)    
    A_1, B_1 = popt_1  # Parámetros ajustados de la función 1
    
    # Ajustar a la función 2
    popt_2, _ = curve_fit(modelo_2, orden_deseado, max_Auto_2, maxfev = 10000)
    A_2, c_2, B_2 = popt_2  # Parámetros ajustados de la función 2

    # Graficar los ajustes
    x_vals = np.linspace(min(orden_deseado), max(orden_deseado), 1000)
    
    # Ajuste de la función 1
    ax_max.plot(x_vals, modelo_1(x_vals, *popt_1), color='red', label=f'A/sqrt(x) + B: A={A_1:.3f}, B={B_1:.3f}')
    
    # Ajuste de la función 2
    ax_max.plot(x_vals, modelo_2(x_vals, *popt_2), color='green', label=f'A/x^c + B: A={A_2:.3f}, c={c_2:.3f}, B={B_2:.3f}')
    
    # Etiquetas y título
    ax_max.set_xlabel('N')
    ax_max.set_ylabel('Maximo')
    ax_max.set_title('Maximo Auto Hopp - Traza/3 ')
    ax_max.legend()
    ax_max.grid(True)

# Mostrar gráfico

plt.tight_layout()
plt.show()   