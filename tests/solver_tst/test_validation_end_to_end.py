"""
@file test_validation_end_to_end.py
@brief System tests validation against Nastran reference cases.
@details
Runs full analysis cycles (Parse -> Assemble -> Solve) on external .dat files
and compares the resulting displacements against reference values extracted 
from MSC Nastran .f06 output files.
"""
import pytest
import os
from pathlib import Path

# Imports de tu proyecto
from src.preprocessing.nastran_reader import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver

# --- DATABASE OF VALIDATION CASES ---
# Format:
# {
#   "case_name": Unique ID,
#   "folder": Subfolder name inside 'tests/validation_cases' (e.g., "case_1"),
#   "filename": .dat file name (must exist in tests/validation_cases),
#   "expected_displacements": { NodeID: { DOF_Index: Value } }
# }
# DOF Index map: 0=Tx, 1=Ty, 2=Tz, 3=Rx, 4=Ry, 5=Rz

VALIDATION_DATA = [
    {
        "case_name": "Ej1_v1_Simple_Tension",
        "folder": "case_1",
        "filename": "Ej1_v1.dat",
        "expected_displacements": {
            1: {0: 0.0},          # Node 1, T1 (Tx)
            2: {0: 1.190476e-3,   # Node 2, T1 (Tx)
                1: 0.0}           # Node 2, T2 (Ty)
        },
        "description": "Single Rod element. Fixed at Node 1. Load at Node 2."
    },
    {
        "case_name": "Ej2_Series_Rods",
        "folder": "case_2",
        "filename": "Ej2.dat",
        "expected_displacements": {
            2: {0: 1.250000e-3}, # Node 2, T1 (Tx)
            3: {0: 4.424603e-3}  # Node 3, T1 (Tx)
        },
        "description": "Two Rod elements in series with different materials."
    },
    {
        "case_name": "Ej3_Planar_Truss",
        "folder": "case_3",
        "filename": "Ej3.dat",
        "expected_displacements": {
            5: {
                0: 9.998600e-4, # Node 5, T1 (Tx)
                1: 1.444375e-4  # Node 5, T2 (Ty)
            }
        },
        "description": "2D Truss (Triangle). Validates rotated coordinate assembly."
    }
]
class TestNastranEndToEnd:
    """
    @class TestNastranEndToEnd
    @brief Automated system testing using parametrization.
    """

    @pytest.mark.parametrize("case_data", VALIDATION_DATA, ids=lambda x: x["case_name"])
    def test_run_nastran_case(self, case_data, validation_data_dir):
        """
        @brief Executes a complete analysis flow for a given Nastran case.
        @param case_data Dictionary containing file info and expected results.
        @param validation_data_dir Pytest fixture pointing to the data folder.
        """
        # 1. Resolve File Path (Supporting Subfolders)
        # Construct path: tests/validation_cases / {folder} / {filename}
        folder = case_data.get("folder", "")
        filename = case_data["filename"]
        
        file_path = validation_data_dir / folder / filename
        
        if not file_path.exists():
            pytest.fail(f"Test file not found: {file_path}")

        # 2. PARSE (Implicit integration test for nastran_reader)
        try:
            model = parse_nastran_input(str(file_path))
        except Exception as e:
            pytest.fail(f"Parser failed on {case_data['filename']}: {e}")

        # 3. SOLVE
        # Note: We assume the .dat files are now "Python-Solver-Ready" 
        # (containing necessary SPCs to avoid singularity), so we don't need manual stabilization here.
        try:
            solver = LinearStaticSolver(model)
            solver.assemble()
            displacements = solver.solve()
        except Exception as e:
            pytest.fail(f"Solver failed on {case_data['filename']}: {e}")

        # 4. VERIFY RESULTS
        expected_disp = case_data.get("expected_displacements", {})
        
        for nid, dofs_dict in expected_disp.items():
            # Get Global Indices
            try:
                global_indices = solver.mapper.get_dof_indices(nid)
            except ValueError:
                pytest.fail(f"Node {nid} expected in results but not found in Solver Mapper.")
            
            for local_dof, val_ref in dofs_dict.items():
                idx = global_indices[local_dof]
                val_calc = displacements[idx]
                
                # Check with tolerance (0.1% relative error is generous but safe for FEM)
                # Absolute tolerance 1e-9 handles near-zero values.
                assert val_calc == pytest.approx(val_ref, rel=1e-3, abs=1e-9), \
                    f"Mismatch in {case_data['case_name']}, Node {nid}, DOF {local_dof}. " \
                    f"Expected {val_ref}, Got {val_calc}"
                

                