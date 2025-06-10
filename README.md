# Codigos-TFM

A continuación presento lo que hace cada código y el orden en el que se deben correr para obtener los resultados mostrados en el TFM. La mayortía son códigos en Matlab, pero también hay códigos en ROOT y Python
Los siguientes códigos no hace falta correrlos ya que Scaling.m los llama directamente. El orden de llamado es:

- Base.m: Calcula la base |n1,n2,n3> del modelo.
- Basesim.m: Calcula la base de autoestados de R.
- Energiasim_mejoras.m Calcula la matriz de H en la base conjunta con R, separado el espectro en los 3 sectores.
- Corriente_mejorada.m: Calcula el operador Corriente C en la base conjunta (ecuación 11 TFM)
- CreacDest12_Mejorado.m: Calcula el operador h_{12} en la base conjunta (ecuación 13 TFM)
- Memoria.m: Calcula el operador I en la base conjunta (ecuación 12 TFM)
- Scaling.m: Calcula la traza y autovalores de los operadores C, I y h_{12} para diferentes valores de N.
Los datos obtenidos por Scaling se guardan en los ficheros:
- ETH_Rotacion_N.txt con N = 20,28,40,56,80,112,160,224,320 (Significado de las columnas: 1- autovalor de H con r=1, 2- autovalor de H con r=e^{2pi/3}, 3- traza subespacio 2x2 degenrado de rotación, 4- parte real del primer autovalor en subespacio 2x2, 5- parte real del segundo autovalor en el subespacio 2x2, 6- parte imaginaria del primer autovalor en subespacio 2x2, 7- parte imaginaria del segundo autovalor en subespacio 2x2, 8- traza subespacio 3x3, 9- parte real del primer autovalor en subespacio 3x3, 10- parte real del segundo autovalor en subespacio 3x3, 11- parte real del tercer autovalor en subespacio 3x3, 12- parte imaginaria del primer autovalor en subespacio 3x3, 13- parte imaginaria del segundo autovalor en subespacio 3x3, 14- parte imaginaria del tercer autovalor en subespacio 3x3, 15- un 3 (antes tenía significado pero conforme avanzó el código se quedó)).
- ETH_Hopping_N.txt con N = 20,28,40,56,80,112,160,224,320 (Significado de las columnas: 1- autovalor de H con r=1, 2- autovalor de H con r=e^{2pi/3}, 3- traza subespacio 2x2 degenrado de rotación, 4- primer autovalor en subespacio 2x2, 5- segundo autovalor en el subespacio 2x2, 6- traza subespacio 3x3, 7- primer autovalor en subespacio 3x3, 8- segundo autovalor en subespacio 3x3, 9- tercer autovalor en subespacio 3x3, 10- un 3 (lo mismo))
- ETH_Corriente_N.txt con N = 20,28,40,56,80,112,160,224,320 (Mismo significado en las columnas que el Hopping)

Utilizamos el código Ajuste_Traza.C escrito para ROOT, para simplificar la información de los ficheros anteriores, centrándonos en la necesaria para realizar la figura (4). En particular de todos las funciones que hay en Ajuste_Traza los códigos necesarios son
-Radio_Rot: Ordena los autovalores de Rotacion, creando los ficheros Orden_N.dat (Significado de las columnas: 1- autovalor de H con r=1, 2.- parte real del autovalor con parte imaginaria casi nula, 3.- Su correspondiente parte imaginaria, 4.- Parte real del autovalor con parte imaginaria claramente positiva, 5.-  Su correspondiente parte imaginaria, 6.- Parte real del autovalor con parte imaginaria claramente negativa, 7.-  Su correspondiente parte imaginaria)
- Positivos_Corriente: Ordena los autovalores de la corriente, creando los ficheros Corriente_Positivos_N.dat, (Significado de las columnas: 1- autovalor de H con r=1, 2.- Autovalor positivo, 3.- Autovalores negativo (El tercer autovalor es nulo))
- Ordenar_Autoval_Hopping: Lo mismo con el Hopping, creando los ficheros ETH_Hopping_Autoval_Full_N.dat, (Significado de las columnas: 1- autovalor de H con r=1, 2.- Autovalor del Hopping) Me coloca los 3 autovalores en una misma columna.

El fichero de valores microcanónicos aportado por el profe es
- micro_hopping.dat
   
Hasta aquí ya podemos representar la figura 2 y 3 completa. Y la parte cuántica de la figura 6

A través de estos ficheros se obtienen los resultados de la figura 4 (a) con los códigos:
- AutoCorr.py: Calcula la sigma y los maximos de la distribucion de autovalores de la corriente en la zona térmica y ajuste a una potencia de N.
- AutoRot.py: Lo mismo con la Rotación
- AutoTrazaHopp.py: Compara la ecuación (8) para el Hopping, haciendo lo mismo que los otros códigos.
- TrazaHopp.py: Compara la ecuación (7) para el Hopping, haciendo lo mismo que los otros códigos.
- AutoHopp.py: Compara los autovalores con el valor microcanónico del Hopping, haciendo lo mismo que los otros códigos.
- TFM.ipynb: Lo mismo que todos los anteriores pero en un notebook.

Para obtener las figuras 4 (b), (c) y (d) usamos
- ETH_Comparacion_Fuerte. Además también calcula una figura análoga para la condición (8) del Hopping.

Por último, los ratios se calculan de los datos
-espectro.t1_mas.u-5.n1000.dat
con el código
Ratios_TFM.m
