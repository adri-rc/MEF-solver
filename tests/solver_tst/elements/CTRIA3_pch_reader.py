"""
@file CTRIA3_stiffness_matrix_pch_tests.py
@brief Validación de la Matriz de Rigidez Elemental contra MSC Nastran (.pch).
"""
import numpy as np
import csv
import os

def leer_dmig_pch(ruta_archivo, nombre_matriz='KAAX', num_nudos=3):
    """
    Lee una matriz en formato DMIG Punch (Large Field) de Nastran 
    y la reconstruye en una matriz NumPy de (N*6) x (N*6).
    """
    total_dofs = num_nudos * 6
    K = np.zeros((total_dofs, total_dofs))
    
    with open(ruta_archivo, 'r') as f:
        lineas = f.readlines()
        
    en_matriz = False
    col_nodo = None
    col_dof = None
    
    for linea in lineas:
        # Detectar inicio o fin de la matriz buscada
        if linea.startswith('DMIG '):
            partes = linea.split()
            if len(partes) > 1 and partes[1] == nombre_matriz:
                en_matriz = True
            elif en_matriz and len(partes) > 1 and partes[1] != nombre_matriz:
                en_matriz = False
                
        # Leer Cabecera de Columna (formato Large Field 16 caracteres)
        elif en_matriz and linea.startswith('DMIG*'):
            col_nodo = int(linea[32:48].strip())
            col_dof = int(linea[48:64].strip())
            
        # Leer Valores de Fila
        elif en_matriz and linea.startswith('*'):
            row_nodo = int(linea[8:24].strip())
            row_dof = int(linea[24:40].strip())
            val_str = linea[40:56].strip()
            
            # Ajuste por si Nastran imprime 1.0-4 en lugar de 1.0E-4
            if 'E' not in val_str:
                if val_str.find('-', 1) != -1:
                    idx = val_str.rfind('-')
                    val_str = val_str[:idx] + 'E' + val_str[idx:]
                elif val_str.find('+', 1) != -1:
                    idx = val_str.rfind('+')
                    val_str = val_str[:idx] + 'E' + val_str[idx:]
            
            valor = float(val_str)
            
            # Calcular índices en Python (0 a 17)
            idx_col = (col_nodo - 1) * 6 + (col_dof - 1)
            idx_row = (row_nodo - 1) * 6 + (row_dof - 1)
            
            # Asignar a la matriz (es simétrica)
            K[idx_row, idx_col] = valor
            K[idx_col, idx_row] = valor 
            
    return K

def exportar_matriz_pch_csv(K, num_nudos=3, nombre_archivo="KAAX_Global.csv"):
    """Exporta la matriz extraída al formato CSV con etiquetas completas."""

    '''
    dof_names =  ['u1', 'v1', 'w1', 'theta_x1', 'theta_y1', 'theta_z1',
                'u2', 'v2', 'w2', 'theta_x2', 'theta_y2', 'theta_z2',
                'u3', 'v3', 'w3', 'theta_x3', 'theta_y3', 'theta_z3']
    '''

    dof_names = []
    for i in range(1, num_nudos + 1):
        dof_names.extend([f'u{i}', f'v{i}', f'w{i}', f'theta_x{i}', f'theta_y{i}', f'theta_z{i}'])

    # Lógica de rutas
    if os.path.isabs(nombre_archivo):
        ruta_completa = nombre_archivo
    else:
        directorio_script = os.path.dirname(os.path.abspath(__file__))
        ruta_completa = os.path.join(directorio_script, nombre_archivo)
        
    directorio_destino = os.path.dirname(ruta_completa)
    if directorio_destino and not os.path.exists(directorio_destino):
        os.makedirs(directorio_destino)
        
    with open(ruta_completa, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter=',')
        
        # 1. Cabecera (Index, DoF, u1, v1, ...)
        writer.writerow(['Index', 'DoF'] + dof_names)
        
        # 2. Filas de datos
        for i in range(len(K)):
            fila_valores = [f"{val:.3e}" if abs(val) > 1e-10 else "0.0" for val in K[i]]
            writer.writerow([i, dof_names[i]] + fila_valores)
            
    print(f"✅ Matriz global extraída del .pch y exportada a: \n   --> '{ruta_completa}'")


# ==========================================
# EJECUCIÓN PRINCIPAL
# ==========================================

input_pch_file_path = r"C:\Users\ADRIAN\Desktop\TFM-repo\MEF-solver\tests\solver_tst\validation_cases\case_manual_cquad\cquad_kglobal_flex_pure.pch"
output_K_csv_file_path = r"C:\Users\ADRIAN\Desktop\TFM-repo\MEF-solver\tests\solver_tst\validation_cases\case_manual_cquad\draft_cquad_kglobal_flex_pure_Nastran.csv"
CASE = 'CQUAD' # 'CTRIA' or 'CQUAD'


# ==========================================
if __name__ == "__main__":
    
    archivo_punch = input_pch_file_path

    if CASE == 'CTRIA':
        nodes_num = 3
    elif CASE == 'CQUAD':
        nodes_num = 4
    else:
        raise ValueError('¿Es un "CTRIA" o un "CQUAD" ?')
    
    # 1. Leemos la matriz Global de Rigidez (KAAX) asumiendo 3 nudos
    K_global_nastran = leer_dmig_pch(archivo_punch, nombre_matriz='KAAX', num_nudos=nodes_num)
    
    # 2. Exportamos a CSV
    exportar_matriz_pch_csv(K_global_nastran, num_nudos=nodes_num, nombre_archivo=output_K_csv_file_path)