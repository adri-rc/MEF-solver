"""
@file run_cquad4_validation.py
@brief Ejecuta el solver CQUAD4 Python y exporta la matriz K a CSV para comparar con Nastran.
"""

import numpy as np
import csv
import os
import sys
from pathlib import Path

# --- CONFIGURACIÓN DEL PATH PARA IMPORTAR 'src' ---
# 1. Obtenemos la ruta absoluta de donde está este script
current_dir = os.path.dirname(os.path.abspath(__file__))

# 2. Subimos niveles hasta llegar a la raíz del proyecto (MEF-solver).
# Como tu script está en: MEF-solver/tests/solver_tst/run_cquad4_validation.py (2 niveles abajo)
# o en MEF-solver/tests/solver_tst/validation_cases/case_manual/... (4 niveles abajo)
# vamos a usar una forma segura: buscar la carpeta que contiene 'src'
project_root = current_dir
while not os.path.exists(os.path.join(project_root, 'src')):
    parent = os.path.dirname(project_root)
    if parent == project_root: # Llegamos a la raíz del disco duro sin encontrarlo
        raise RuntimeError("No se pudo encontrar el directorio raíz del proyecto que contiene 'src'.")
    project_root = parent

# 3. Añadimos la raíz al inicio del path de Python
if project_root not in sys.path:
    sys.path.insert(0, project_root)

VALIDATION_ROOT = Path(project_root) / "tests" / "solver_tst" / "validation_cases"

# Importa tu solver (ajusta esta ruta según la estructura de tu proyecto)
from src.solver.elements.cquad4_element import CQUAD4ElementSolver

# --- MOCKS PARA REPLICAR TU DATA_MODEL ---
class MockNode:
    def __init__(self, nid, coords):
        self.nid = nid
        self._coords = np.array(coords, dtype=float)
    def get_coords(self):
        return self._coords

class MockMat:
    def __init__(self, E, nu):
        self.E = E
        self.nu = nu

class MockProp:
    def __init__(self, t):
        self.t = t

def exportar_matriz_python_csv(K, num_nudos=4, nombre_archivo="K_cquad4_python.csv"):
    """Exporta la matriz K del solver Python al mismo formato CSV que el conversor de Nastran."""
    
    dof_names = []
    for i in range(1, num_nudos + 1):
        dof_names.extend([f'u{i}', f'v{i}', f'w{i}', f'theta_x{i}', f'theta_y{i}', f'theta_z{i}'])

    # Lógica de rutas
    if not os.path.isabs(nombre_archivo):
        directorio_script = os.path.dirname(os.path.abspath(__file__))
        nombre_archivo = os.path.join(directorio_script, nombre_archivo)
        
    with open(nombre_archivo, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter=',')
        
        # Cabecera
        writer.writerow(['Index', 'DoF'] + dof_names)
        
        # Filas
        for i in range(len(K)):
            # Formateamos a notación científica y filtramos ruido numérico
            fila_valores = [f"{val:.3e}" if abs(val) > 1e-10 else "0.0" for val in K[i]]
            writer.writerow([i, dof_names[i]] + fila_valores)
            
    print(f"✅ Matriz de Python exportada con éxito a: \n   --> '{nombre_archivo}'")

def run_validation():
    print("--- INICIANDO VALIDACIÓN CQUAD4 (T=I, Membrana) ---")
    
    # 1. Definir los Nodos (Cuadrado 2x2 en plano XY)
    n1 = MockNode(1, [0.0, 0.0, 0.0])
    n2 = MockNode(2, [2.0, 0.0, 0.0])
    n3 = MockNode(3, [2.0, 2.0, 0.0])
    n4 = MockNode(4, [0.0, 2.0, 0.0])
    
    # 2. Definir Propiedades y Materiales
    t = 1.0
    prop = MockProp(t=t)
    mat_mem = MockMat(E=0.0, nu=0.3)  # E=210000, nu=0.3
    mat_flex = MockMat(E=2.1E5, nu=0.3)  # E=210000, nu=0.3
    mat_shear = MockMat(E=0.0, nu=0.3)  # E=210000, nu=0.3
    
    # 3. Inicializar tu Solver
    # Pasamos solo mat_membrane para emular que MID2 y MID3 están en blanco
    solver = CQUAD4ElementSolver(
        node1=n1, node2=n2, node3=n3, node4=n4,
        prop=prop, mat_membrane=mat_mem, mat_bending=mat_flex, mat_shear=mat_shear,
        use_membrane_SRI=False
    )
    
    # 4. Obtener la Matriz de Rigidez
    K_global = solver.get_stiffness_matrix()
    
    # 5. AJUSTE K6ROT=0.0
    # En tu solver se añade un k_diag_small en los DOFs 5, 11, 17, 23.
    # Como en Nastran definimos PARAM, K6ROT, 0.0, anulamos estos términos 
    # en Python para que el CSV coincida exactamente a nivel numérico.
    drill_idx = [5, 11, 17, 23]
    for idx in drill_idx:
        K_global[idx, idx] = 0.0
        
    # 6. Exportar a CSV
    output_path = VALIDATION_ROOT / "case_manual_cquad" / "K_cquad4_pure_flex_integracion_completa_python_v4.csv"
    exportar_matriz_python_csv(K_global, num_nudos=4, nombre_archivo=output_path)



def run_case_51():
    print("--- CASE 51 ---")
    output_path = VALIDATION_ROOT / "case_51" / "draft_ej51_pure_membrane_integracion_completa_python.csv"
    
    # 1. Definir los Nodos (Cuadrado 2x2 en plano XY)
    n1 = MockNode(1, [0.1, 0.2, 1.0])
    n2 = MockNode(2, [2.3, 0.1, 1.0])
    n3 = MockNode(3, [2.1, 2.2, 1.0])
    n4 = MockNode(4, [0.2, 1.9, 1.0])
    
    # 2. Definir Propiedades y Materiales
    t = 1.0
    
    mat_mem = MockMat(E=210000, nu=0.3)  # E=210000, nu=0.3
    mat_flex = MockMat(E=0.0, nu=0.3)  # E=0.0, nu=0.3
    mat_shear = MockMat(E=0.0, nu=0.3)  # E=0.0, nu=0.3
    use_SRI = False
    
    prop = MockProp(t=t)
    # 3. Inicializar tu Solver
    # Pasamos solo mat_membrane para emular que MID2 y MID3 están en blanco
    solver = CQUAD4ElementSolver(
        node1=n1, node2=n2, node3=n3, node4=n4,
        prop=prop, mat_membrane=mat_mem, mat_bending=mat_flex, mat_shear=mat_shear,
        use_membrane_SRI=use_SRI
    )
    
    # 4. Obtener la Matriz de Rigidez
    K_global = solver.get_stiffness_matrix()
    
    # 5. AJUSTE K6ROT=0.0
    # En tu solver se añade un k_diag_small en los DOFs 5, 11, 17, 23.
    # Como en Nastran definimos PARAM, K6ROT, 0.0, anulamos estos términos 
    # en Python para que el CSV coincida exactamente a nivel numérico.
    drill_idx = [5, 11, 17, 23]
    for idx in drill_idx:
        K_global[idx, idx] = 0.0
        
    # 6. Exportar a CSV
    exportar_matriz_python_csv(K_global, num_nudos=4, nombre_archivo=output_path)


def run_case_52():
    print("--- CASE 52 ---")
    output_path = VALIDATION_ROOT / "case_52" / "draft_ej52_k_pure_flex_mid3_0_k6rot_0_int_completa_python_CAMBIO2.csv"
    
    # 1. Definir los Nodos (Cuadrado 2x2 en plano XY)
    n1 = MockNode(1, [0.1, 0.2, 1.0])
    n2 = MockNode(2, [2.3, 0.1, 1.0])
    n3 = MockNode(3, [2.1, 2.2, 1.0])
    n4 = MockNode(4, [0.2, 1.9, 1.0])
    
    # 2. Definir Propiedades y Materiales
    t = 1.0
    
    mat_mem = MockMat(E=0.0, nu=0.3)  # E=0.0, nu=0.3
    mat_flex = MockMat(E=210000, nu=0.3)  # E=210000, nu=0.3
    mat_shear = MockMat(E=0.0, nu=0.3)  # E=0.0, nu=0.3
    use_SRI = False
    
    prop = MockProp(t=t)
    # 3. Inicializar tu Solver
    # Pasamos solo mat_membrane para emular que MID2 y MID3 están en blanco
    solver = CQUAD4ElementSolver(
        node1=n1, node2=n2, node3=n3, node4=n4,
        prop=prop, mat_membrane=mat_mem, mat_bending=mat_flex, mat_shear=mat_shear,
        use_membrane_SRI=use_SRI
    )
    
    # 4. Obtener la Matriz de Rigidez
    K_global = solver.get_stiffness_matrix()
    
    # 5. AJUSTE K6ROT=0.0
    # En tu solver se añade un k_diag_small en los DOFs 5, 11, 17, 23.
    # Como en Nastran definimos PARAM, K6ROT, 0.0, anulamos estos términos 
    # en Python para que el CSV coincida exactamente a nivel numérico.
    drill_idx = [5, 11, 17, 23]
    for idx in drill_idx:
        K_global[idx, idx] = 0.0
        
    # 6. Exportar a CSV
    exportar_matriz_python_csv(K_global, num_nudos=4, nombre_archivo=output_path)



if __name__ == "__main__":
    run_case_52()
