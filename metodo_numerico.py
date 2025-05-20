from ecuaciones_modelo import *
import pandas as pd
import os

# Constantes actualizadas con k = 0.51
diccionario_constantes  = {
    'L': 38,
    'd': 20,
    'rho_0': 5,
    'rho_b': 9.36,
    'g': 9.8184,
    'E': 0.4,
    'U': 1,
    'rho_p': 2.8,
    'Pv': 600,
    'k': 0.51,
    'P61': 0.02
}

# === CONFIGURACIÓN DEL ARCHIVO DE ENTRADA ===
archivo_entrada = 'data_2024.csv'  # Cambia esto si tu archivo tiene otro nombre o ruta
archivo_salida = 'resultados_2024.csv'

# Leer el archivo de entrada
data = pd.read_csv(archivo_entrada)

# Borrar el archivo de salida si ya existe
if os.path.exists(archivo_salida):
    os.remove(archivo_salida)

# Procesar cada fila
for index, row in data.iterrows():
    try:
        Pbruta = pd.to_numeric(row['Potencia(kW)'], errors='coerce')
        rpm = 12.2  # Fijo, ya definido por ti
        p_solidos = pd.to_numeric(row['Porcentaje_de_solidos'], errors='coerce') / 100
        fecha = row.get('Fecha', index)

        if pd.isna(p_solidos) or p_solidos == 0 or pd.isna(Pbruta):
            continue

        results = funcion_morrell(fecha, rpm, Pbruta, p_solidos, diccionario_constantes)

        pd.DataFrame([results]).to_csv(archivo_salida, mode='a', header=not os.path.exists(archivo_salida), index=False)

    except (ZeroDivisionError, ValueError, TypeError) as e:
        print(f"Error en fila {index}: {e}")
        continue

print(f"Proceso finalizado. Resultados en: {archivo_salida}")
