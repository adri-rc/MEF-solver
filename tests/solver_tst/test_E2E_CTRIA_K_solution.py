import pytest
import numpy as np
import pandas as pd

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver

# ==========================================
# DATOS DE VALIDACIÓN DE MATRICES (K_global)
# ==========================================
STIFFNESS_VALIDATION_DATA = [
    #region --- CASO manual: CON K6ROT = 0.0 (MEMBRANA + FLEXION PURA (MID3-> E=0)) ---
    {
        "case_name": "Ej_MANUAL_KGlobal_CTRIA3_K6ROT_0.0",
        "folder": "case_k6rot_ctria", 
        "filename": "Ej_manual_Kglobal_MID3_cero_K6rot_0.0.dat",
        "k_csv_file": "ej_manual_kglobal_mid3_cero_k6rot_0.0_Nastran.csv", # CSV generado con tu script desde el .pch
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, # 
        "debug_export": True, # Activa/desactiva la exportación para este caso
        "output_k_file": "draft_Ej_manual_Kglobal_MID3_cero_K6rot_0.0_python_test.csv"   # Nombre del archivo donde se guardará la matriz
    
    },
    #endregion
     #region --- CASO manual: CON K6ROT = 100.0 (MEMBRANA + FLEXION PURA (MID3-> E=0)) ---
    {
        "case_name": "Ej_MANUAL_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_k6rot_ctria", 
        "filename": "Ej_manual_Kglobal_MID3_cero_K6rot_100.0.dat",
        "k_csv_file": "ej_manual_kglobal_mid3_cero_k6rot_100.0_Nastran.csv", # CSV generado con tu script desde el .pch
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, # 
        "debug_export": True, # Activa/desactiva la exportación para este caso
        "output_k_file": "draft_Ej_manual_Kglobal_MID3_cero_K6rot_100.0_python_test.csv"   # Nombre del archivo donde se guardará la matriz
    
    },
    #endregion
    #region --- CASO 32: CON K6ROT = 0.0 ---
    {
        "case_name": "Ej32_KGlobal_CTRIA3_K6ROT_0.0",
        "folder": "case_32", 
        "filename": "Ej32_Ksolution_K6ROT_0.0.dat",
        "k_csv_file": "ej32_ksolution_k6rot_0.0_Nastran.csv", # CSV generado con tu script desde el .pch
        "description": "Validación de K_global para 2 elementos en plano YZ.",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
        "debug_export": True,                  # Activa/desactiva la exportación para este caso
        "output_k_file": "draft_Ej32_KGlobal_CTRIA3_K6ROT_0.0_python_test.csv"   # Nombre del archivo donde se guardará la matriz
    
    },
    #endregion
    #region --- CASO 32: CON K6ROT = 100.0 ---
    {
        "case_name": "Ej32_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_32", 
        "filename": "Ej32_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej32_ksolution_k6rot_100_Nastran.csv", # CSV generado con tu script desde el .pch
        "description": "Validación de K_global para 2 elementos en plano YZ.",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
        "debug_export": True,                  # Activa/desactiva la exportación para este caso
        "output_k_file": "draft_Ej32_KGlobal_CTRIA3_K6ROT_100.0_python_test.csv"   # Nombre del archivo donde se guardará la matriz
    
    },
    #endregion
    #region --- CASO 33: CON K6ROT = 100.0 ---
    {
        "case_name": "Ej33_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_33", 
        "filename": "Ej33_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej33_ksolution_k6rot_100_Nastran.csv", # CSV generado con tu script desde el .pch
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
        "debug_export": False,                  # Activa/desactiva la exportación para este caso
        "output_k_file": "draft_Ej33_KGlobal_CTRIA3_K6ROT_100.0_python_test.csv"   # Nombre del archivo donde se guardará la matriz
    
    },
    #endregion
    #region --- CASO 34: CON K6ROT = 100.0 ---
    {
        "case_name": "Ej34_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_34", 
        "filename": "Ej34_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej34_ksolution_k6rot_100_Nastran.csv", # CSV generado con tu script desde el .pch
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
        "debug_export": True,                  # Activa/desactiva la exportación para este caso
        "output_k_file": "draft_Ej34_KGlobal_CTRIA3_K6ROT_100.0_python_test.csv"   # Nombre del archivo donde se guardará la matriz
    
    },
    #endregion
    #region --- CASO 35: CON K6ROT = 100.0 ---
    {
        "case_name": "Ej35_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_35", 
        "filename": "Ej35_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej35_ksolution_k6rot_100_Nastran.csv", # CSV generado con tu script desde el .pch
        "description": "Test del K6ROT",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
        "debug_export": True,                  # Activa/desactiva la exportación para este caso
        "output_k_file": "draft_Ej35_KGlobal_CTRIA3_K6ROT_100_python_test.csv"   # Nombre del archivo donde se guardará la matriz
    
    },
    #endregion
    #region --- CASO 36: CON K6ROT = 100.0 ---
    {
        "case_name": "Ej36_KGlobal_CTRIA3_OFFSET_K6ROT_100.0",
        "folder": "case_36", 
        "filename": "Ej36_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej36_ksolution_k6rot_100_Nastran.csv", # CSV generado con tu script desde el .pch
        "description": "Ejemplo con Offset de cordenadas, y enormes dimensiones frente al espesor.",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
        "debug_export": True,                  # Activa/desactiva la exportación para este caso
        "output_k_file": "draft_Ej36_KGlobal_CTRIA3_K6ROT_100_python_test.csv"   # Nombre del archivo donde se guardará la matriz
    },
    #endregion
    
    # --- CASO 37 (tiene más de un ctria): CON K6ROT = 100.0 ---
    # --- CASO 38 (tiene más de un ctria): CON K6ROT = 100.0 ---
    # --- CASO 39 (tiene más de un ctria): CON K6ROT = 100.0 ---
    # --- CASO 40 (tiene más de un ctria): CON K6ROT = 100.0 ---
    # --- CASO 41 (tiene más de un ctria): CON K6ROT = 100.0 ---
    # --- CASO 42 (tiene más de un ctria): CON K6ROT = 100.0 ---
    # --- CASO 43 (tiene más de un ctria): CON K6ROT = 100.0 ---
    #region --- CASO 44: CON K6ROT = 100.0 ---
    {
        "case_name": "Ej44_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_44", 
        "filename": "Ej44_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej44_ksolution_k6rot_100_Nastran.csv", # CSV generado con tu script desde el .pch
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
        "debug_export": True,                  # Activa/desactiva la exportación para este caso
        "output_k_file": "draft_Ej44_KGlobal_CTRIA3_K6ROT_100.0_python_test.csv"   # Nombre del archivo donde se guardará la matriz
    
    },
    #endregion
    # --- CASO 45 (tiene más de un ctria): CON K6ROT = 100.0 ---
    #region --- CASO 46: CON K6ROT = 100.0 ---
    {
        "case_name": "Ej46_KGlobal_CTRIA3_K6ROT_100",
        "folder": "case_46", 
        "filename": "Ej46_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej46_ksolution_k6rot_100_Nastran.csv", # CSV generado con tu script desde el .pch
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
        "debug_export": True,                  # Activa/desactiva la exportación para este caso
        "output_k_file": "draft_Ej46_KGlobal_CTRIA3_K6ROT_100_python_test.csv"   # Nombre del archivo donde se guardará la matriz
    
    },
    #endregion    
]

class TestStiffnessMatrix:
    """
    @class TestStiffnessMatrix
    @brief Compara la Matriz de Rigidez Global ensamblada por el solver 
           de Python frente a la matriz DMIGPCH exportada por MSC Nastran.
    """

    def _exportar_matriz_a_csv(self, K, ruta_salida, num_nudos=3):
        """
        @brief Exporta la matriz K en formato CSV puro (Index,DoF,u1,v1,...)
               sin banners ni adornos para comparación directa.
        """
        dof_names = []
        for i in range(1, num_nudos + 1):
            dof_names.extend([f'u{i}', f'v{i}', f'w{i}', f'theta_x{i}', f'theta_y{i}', f'theta_z{i}'])

        with open(ruta_salida, 'w', encoding='utf-8') as f:
            # 1. Cabecera limpia
            header = "Index,DoF," + ",".join(dof_names)
            f.write(header + "\n")
            
            # 2. Filas separadas estrictamente por comas
            for i in range(len(K)):
                fila_valores = [f"{val:.9e}" if abs(val) > 1e-10 else "0.0" for val in K[i]]
                f.write(f"{i},{dof_names[i]}," + ",".join(fila_valores) + "\n")
                
        print(f"\n[DEBUG] 💾 Matriz CSV de Python exportada a: {ruta_salida}")
    
    
    @pytest.mark.parametrize("case_data", STIFFNESS_VALIDATION_DATA, ids=lambda x: x["case_name"])
    def test_global_stiffness_matrix(self, case_data, validation_data_dir):
        """
        @brief Ejecuta el ensamblaje y compara K_global contra Nastran.
        """
        folder = case_data.get("folder", "")
        filename = case_data["filename"]
        csv_filename = case_data["k_csv_file"]
        
        dat_file_path = validation_data_dir / folder / filename
        csv_file_path = validation_data_dir / folder / csv_filename
        
        if not dat_file_path.exists():
            pytest.fail(f"Falta archivo .dat: {dat_file_path}")
        if not csv_file_path.exists():
            pytest.fail(f"Falta archivo .csv de Nastran: {csv_file_path}")

        # 1. Leer el CSV exportado de Nastran
        # Asumimos que tu script exporta con las columnas ['Index', 'DoF', 'u1', 'v1', ...]
        df_nastran = pd.read_csv(csv_file_path)
        # Nos quedamos solo con los valores numéricos de la matriz (quitamos 'Index' y 'DoF')
        K_nastran = df_nastran.iloc[:, 2:].to_numpy(dtype=float)

        # 2. Parsear el modelo en Python
        try:
            model = parse_nastran_input(str(dat_file_path))
        except Exception as e:
            pytest.fail(f"El parser falló en {filename}: {e}")

        # 3. Ensamblar la matriz con tu Solver
        try:
            solver = LinearStaticSolver(model)
            solver.assemble()  # Esto debe construir solver.K_global
        except Exception as e:
            pytest.fail(f"Fallo al ensamblar en {filename}: {e}")

        # Extraer K_global del solver (si es sparse lil_matrix o csr_matrix, la pasamos a densa)
        if hasattr(solver, 'K'):
            K_python = solver.K.toarray() if hasattr(solver.K, 'toarray') else solver.K
        elif hasattr(solver, 'K_global'):
            K_python = solver.K_global.toarray() if hasattr(solver.K_global, 'toarray') else solver.K_global
        else:
            pytest.fail("No se encontró el atributo de la matriz global en el solver (ej. solver.K o solver.K_global).")

        # [NUEVO] COMPROBACIÓN INDIVIDUAL DE DEPURACIÓN Y EXPORTACIÓN
        # ---------------------------------------------------------------------
        if case_data.get("debug_export", False) and case_data.get("output_k_file"):
            ruta_salida = validation_data_dir / folder / case_data["output_k_file"]
            nodes_num = K_python.shape[0] // 6
            # Llamamos al nuevo método de volcado CSV
            self._exportar_matriz_a_csv(K_python, ruta_salida, num_nudos=nodes_num)
        # ---------------------------------------------------------------------
        
        # 4. Verificaciones de dimensiones
        assert K_python.shape == K_nastran.shape, \
            f"Dimensión de matriz K no coincide. Python: {K_python.shape}, Nastran: {K_nastran.shape}"

        # 5. Comparación elemento a elemento
        rel_tol = case_data.get("rel_tol", 1e-4)
        abs_tol = case_data.get("abs_tol", 1e-6)

        # Identificamos los índices donde hay diferencias significativas
        diferencias = []
        n_rows, n_cols = K_python.shape
        
        for i in range(n_rows):
            for j in range(n_cols):
                val_py = K_python[i, j]
                val_nas = K_nastran[i, j]
                
                # Usamos isclose de numpy para tolerancia robusta
                if not np.isclose(val_py, val_nas, rtol=rel_tol, atol=abs_tol):
                    diferencias.append({
                        "Fila_Global": i,
                        "Columna_Global": j,
                        "Val_Python": val_py,
                        "Val_Nastran": val_nas,
                        "Error_Abs": abs(val_py - val_nas)
                    })

        # Si hay diferencias, fallamos el test mostrando los peores términos
        if diferencias:
            # Ordenar por el error absoluto de mayor a menor para depurar mejor
            diferencias.sort(key=lambda x: x["Error_Abs"], reverse=True)
            top_diffs = diferencias[:10]  # Mostrar los 10 peores fallos
            
            error_msg = f"Las matrices K_global no coinciden en {case_data['case_name']}.\n"
            error_msg += f"Total de términos con error mayor a tolerancia: {len(diferencias)}\n"
            error_msg += "Top 10 mayores diferencias:\n"
            for diff in top_diffs:
                error_msg += f"Idx[{diff['Fila_Global']},{diff['Columna_Global']}] -> "
                error_msg += f"Py: {diff['Val_Python']:.4e} | Nastran: {diff['Val_Nastran']:.4e} | Diff: {diff['Error_Abs']:.4e}\n"
            
            pytest.fail(error_msg)