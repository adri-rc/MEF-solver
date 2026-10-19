import pytest
import numpy as np
import pandas as pd

from src.preprocessing.nastran_reader import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver

STIFFNESS_VALIDATION_DATA = [
    # Manual case: K6ROT=0, membrane plus pure bending (MID3 E=0).
    {
        "case_name": "Ej_MANUAL_KGlobal_CTRIA3_K6ROT_0.0",
        "folder": "case_k6rot_ctria", 
        "filename": "Ej_manual_Kglobal_MID3_cero_K6rot_0.0.dat",
        "k_csv_file": "ej_manual_kglobal_mid3_cero_k6rot_0.0_Nastran.csv", # CSV generated from the PCH reference
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1,
    
    },
    # Manual case: K6ROT=100, membrane plus pure bending (MID3 E=0).
    {
        "case_name": "Ej_MANUAL_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_k6rot_ctria", 
        "filename": "Ej_manual_Kglobal_MID3_cero_K6rot_100.0.dat",
        "k_csv_file": "ej_manual_kglobal_mid3_cero_k6rot_100.0_Nastran.csv", # CSV generated from the PCH reference
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1,
    
    },
    # Case 32 with K6ROT=0.
    {
        "case_name": "Ej32_KGlobal_CTRIA3_K6ROT_0.0",
        "folder": "case_32", 
        "filename": "Ej32_Ksolution_K6ROT_0.0.dat",
        "k_csv_file": "ej32_ksolution_k6rot_0.0_Nastran.csv", # CSV generated from the PCH reference
        "description": "Global stiffness of two elements in the YZ plane.",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
    
    },
    # Case 32 with K6ROT=100.
    {
        "case_name": "Ej32_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_32", 
        "filename": "Ej32_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej32_ksolution_k6rot_100_Nastran.csv", # CSV generated from the PCH reference
        "description": "Global stiffness of two elements in the YZ plane.",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
    
    },
    # Case 33 with K6ROT=100.
    {
        "case_name": "Ej33_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_33", 
        "filename": "Ej33_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej33_ksolution_k6rot_100_Nastran.csv", # CSV generated from the PCH reference
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
    
    },
    # Case 34 with K6ROT=100.
    {
        "case_name": "Ej34_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_34", 
        "filename": "Ej34_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej34_ksolution_k6rot_100_Nastran.csv", # CSV generated from the PCH reference
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
    
    },
    # Case 35 with K6ROT=100.
    {
        "case_name": "Ej35_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_35", 
        "filename": "Ej35_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej35_ksolution_k6rot_100_Nastran.csv", # CSV generated from the PCH reference
        "description": "K6ROT validation.",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
    
    },
    # Case 36 with K6ROT=100.
    {
        "case_name": "Ej36_KGlobal_CTRIA3_OFFSET_K6ROT_100.0",
        "folder": "case_36", 
        "filename": "Ej36_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej36_ksolution_k6rot_100_Nastran.csv", # CSV generated from the PCH reference
        "description": "Large coordinate offset and dimensions relative to thickness.",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
    },
    # Case 44 with K6ROT=100.
    {
        "case_name": "Ej44_KGlobal_CTRIA3_K6ROT_100.0",
        "folder": "case_44", 
        "filename": "Ej44_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej44_ksolution_k6rot_100_Nastran.csv", # CSV generated from the PCH reference
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
    
    },
    # Case 46 with K6ROT=100.
    {
        "case_name": "Ej46_KGlobal_CTRIA3_K6ROT_100",
        "folder": "case_46", 
        "filename": "Ej46_Ksolution_K6ROT_100.dat",
        "k_csv_file": "ej46_ksolution_k6rot_100_Nastran.csv", # CSV generated from the PCH reference
        "description": "",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1, 
    
    },
    {
        "case_name": "CTRIA3_PSHELL_TS_T_0.25",
        "folder": "case_ctria3_pshell_factors_k",
        "filename": "ctria3_k_ts_t_025.dat",
        "k_csv_file": "ctria3_k_ts_t_025_Nastran.csv",
        "description": "DMIG control for PSHELL TS/T=0.25.",
        "rel_tol": 2e-3,
        "abs_tol": 1e-1,
        "max_relative_frobenius_error": 1.1e-4,
        "matrix_terms": [
            (3, 3), (3, 4), (3, 9), (3, 10),
            (4, 4), (4, 15), (4, 16), (10, 10), (15, 15),
        ],
    },
    {
        "case_name": "CTRIA3_PSHELL_TS_T_1.25",
        "folder": "case_ctria3_pshell_factors_k",
        "filename": "ctria3_k_ts_t_125.dat",
        "k_csv_file": "ctria3_k_ts_t_125_Nastran.csv",
        "description": "DMIG control for PSHELL TS/T=1.25.",
        "rel_tol": 2e-3,
        "abs_tol": 1e-1,
        "max_relative_frobenius_error": 5.1e-4,
        "matrix_terms": [
            (3, 3), (3, 4), (3, 9), (3, 10),
            (4, 4), (4, 15), (4, 16), (10, 10), (15, 15),
        ],
    },
    {
        "case_name": "CTRIA3_PSHELL_12I_T3_0.40",
        "folder": "case_ctria3_pshell_factors_k",
        "filename": "ctria3_k_12i_t3_040.dat",
        "k_csv_file": "ctria3_k_12i_t3_040_Nastran.csv",
        "description": "DMIG control for PSHELL 12I/T^3=0.40.",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1,
    },
    {
        "case_name": "CTRIA3_PSHELL_12I_T3_1.60",
        "folder": "case_ctria3_pshell_factors_k",
        "filename": "ctria3_k_12i_t3_160.dat",
        "k_csv_file": "ctria3_k_12i_t3_160_Nastran.csv",
        "description": "DMIG control for PSHELL 12I/T^3=1.60.",
        "rel_tol": 1e-4,
        "abs_tol": 1e-1,
    },
]

class TestStiffnessMatrix:
    """Compare assembled stiffness matrices with MSC Nastran DMIG data."""
    
    
    @pytest.mark.parametrize("case_data", STIFFNESS_VALIDATION_DATA, ids=lambda x: x["case_name"])
    def test_global_stiffness_matrix(self, case_data, validation_data_dir):
        """Assemble ``K_global`` and compare it with the Nastran reference."""
        folder = case_data.get("folder", "")
        filename = case_data["filename"]
        csv_filename = case_data["k_csv_file"]
        
        case_folder_path = validation_data_dir.parent.parent /'validation_cases'/ folder
        dat_file_path = case_folder_path / filename
        csv_file_path = case_folder_path / csv_filename
        
        if not dat_file_path.exists():
            pytest.fail(f"Missing DAT file: {dat_file_path}")
        if not csv_file_path.exists():
            pytest.fail(f"Missing Nastran CSV file: {csv_file_path}")

        # The reference CSV begins with the Index and DoF label columns.
        df_nastran = pd.read_csv(csv_file_path)
        K_nastran = df_nastran.iloc[:, 2:].to_numpy(dtype=float)

        # Parse and assemble the model.
        try:
            model = parse_nastran_input(str(dat_file_path))
        except Exception as e:
            pytest.fail(f"Parser failed for {filename}: {e}")

        try:
            solver = LinearStaticSolver(model)
            solver.assemble()
        except Exception as e:
            pytest.fail(f"Assembly failed for {filename}: {e}")

        K_python = solver.K_global.toarray()
        
        # Matrix dimensions.
        assert K_python.shape == K_nastran.shape, \
            f"Stiffness shape mismatch. Python: {K_python.shape}, Nastran: {K_nastran.shape}"

        max_frobenius_error = case_data.get("max_relative_frobenius_error")
        if max_frobenius_error is not None:
            reference_norm = np.linalg.norm(K_nastran)
            relative_frobenius_error = (
                np.linalg.norm(K_python - K_nastran) / reference_norm
            )
            assert relative_frobenius_error <= max_frobenius_error, (
                f"Excessive relative Frobenius error in {case_data['case_name']}: "
                f"{relative_frobenius_error:.6e} > {max_frobenius_error:.6e}"
            )

        # Term-by-term comparison.
        rel_tol = case_data.get("rel_tol", 1e-4)
        abs_tol = case_data.get("abs_tol", 1e-6)

        differences = []
        n_rows, n_cols = K_python.shape
        matrix_terms = case_data.get("matrix_terms")
        if matrix_terms is None:
            matrix_terms = (
                (i, j) for i in range(n_rows) for j in range(n_cols)
            )

        for i, j in matrix_terms:
            val_py = K_python[i, j]
            val_nas = K_nastran[i, j]

            if not np.isclose(val_py, val_nas, rtol=rel_tol, atol=abs_tol):
                differences.append({
                    "global_row": i,
                    "global_column": j,
                    "python_value": val_py,
                    "nastran_value": val_nas,
                    "absolute_error": abs(val_py - val_nas),
                })

        if differences:
            differences.sort(key=lambda item: item["absolute_error"], reverse=True)
            top_diffs = differences[:10]
            
            error_msg = f"K_global mismatch in {case_data['case_name']}.\n"
            error_msg += (
                "Terms outside tolerance: "
                f"{len(differences)}\nLargest 10 differences:\n"
            )
            for diff in top_diffs:
                error_msg += (
                    f"Idx[{diff['global_row']},{diff['global_column']}] -> "
                    f"Py: {diff['python_value']:.4e} | "
                    f"Nastran: {diff['nastran_value']:.4e} | "
                    f"Diff: {diff['absolute_error']:.4e}\n"
                )
            
            pytest.fail(error_msg)
