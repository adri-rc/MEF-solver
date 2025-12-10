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
#   "expected_reactions":     { NodeID: { DOF_Index: Value } }
#   "expected_element_forces":{ ElemID: { 'axial_force': Val, 'axial_stress': Val } }
# }
# DOF Index map: 0=Tx, 1=Ty, 2=Tz, 3=Rx, 4=Ry, 5=Rz

VALIDATION_DATA = [
   # --- CASO 1: Barra Simple (Reference) ---
    {
        "case_name": "Ej1_v1_Simple_Tension",
        "folder": "case_1",
        "filename": "Ej1_v1.dat",
        "expected_displacements": {
            1: {0: 0.0},          # Node 1, T1 (Tx)
            2: {0: 1.190476e-3,   # Node 2, T1 (Tx)
                1: 0.0}           # Node 2, T2 (Ty)
        },
        "expected_reactions": {
            1: {0: -5.000000E+02} # Node 1, Tx. F06: -500.0
        },
        "expected_element_forces": {
            101: {"axial_force": 5.000000E+02, "axial_stress": 5.000000E+00}
        },
        "description": "Single Rod element. Fixed at Node 1. Load at Node 2."
    },

    # --- CASO 2: Barras en Serie (Multi-material) ---
    {
        "case_name": "Ej2_Series_Rods",
        "folder": "case_2",
        "filename": "Ej2.dat",
        "expected_displacements": {
            2: {0: 1.250000e-3}, # Node 2, T1 (Tx)
            3: {0: 4.424603e-3}  # Node 3, T1 (Tx)
        },
        "expected_reactions": {
            1: {0: -5.000000E+02} # Node 1 takes full load
        },
        "expected_element_forces": {
            101: {"axial_force": 5.000000E+02, "axial_stress": 5.000000E+00},
            102: {"axial_force": 5.000000E+02, "axial_stress": 2.222222E+00}
        },
        "description": "Two Rod elements in series with different materials."
    },

    # --- CASO 3: Cercha Plana (Rotated Coords) ---
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
        "expected_reactions": {
            # Reference: ej3.f06 PAGE 12
            # Node 1: T1=-1.250108E+02, T2=-2.165187E+02
            1: {0: -1.250108E+02, 1: -2.165187E+02},
            # Node 4: T1=-3.749892E+02, T2= 2.165187E+02
            4: {0: -3.749892E+02, 1:  2.165187E+02}
        },
        "expected_element_forces": {
            # Reference: ej3.f06 PAGE 13 & 14
            104: {"axial_force": 2.500162E+02, "axial_stress": 2.500162E+00},
            105: {"axial_force": 4.330095E+02, "axial_stress": 1.924487E+00}
        },
        "description": "2D Truss (Triangle). Validates rotated coordinate assembly."
    },

    # --- CASO 4: El Riel 3D (3D Rail Effect) ---
    # Validates projection of stiffness on a constrained sliding node.
    {
        "case_name": "Ej4_Rail_3D",
        "folder": "case_4",
        "filename": "Ej4.dat",
        "expected_displacements": { 
            2: {0: 5.633333E+00, 1: 0.0, 2: 0.0} # X moves, Y/Z fixed
        },
        "expected_reactions": { 
            1: {0: -300.0, 1: -400.0, 2: -1200.0}, # Wall takes full load
            2: {1: 0.0, 2: 0.0} # Rail takes NO load (perfect alignment)
        },
        "expected_element_forces": { 
            104: {"axial_force": 1.300000E+03, "axial_stress": 1.300000E+03}
        }
    },

    # --- CASO 5: CSV Parser Test ---
    # Physically identical to Ej2 but input is Comma-Separated.
    {
        "case_name": "Ej5_CSV_Parser",
        "folder": "case_5",
        "filename": "Ej5.dat",
        "expected_displacements": { 
            2: {0: 2.380952E-03}, 
            3: {0: 4.761905E-03} 
        },
        # Force is 1000 in both elements
        "expected_element_forces": {
            105: {"axial_force": 1.000000E+03, "axial_stress": 1.000000E+01},
            115: {"axial_force": 1.000000E+03, "axial_stress": 1.000000E+01}
        }
    },

    # --- CASO 6: Torsion Oblicua ---
    # Validates Moment -> Torque projection on rotated element.
    {
        "case_name": "Ej6_Oblique_Torsion",
        "folder": "case_6",
        "filename": "Ej6.dat",
        "expected_displacements": { 
            2: {3: 4.062500E-01} # Rotation X (Radiants)
        },
        "expected_element_forces": { 
            # Element 106. Torque is amplified by 1/cos(alpha) = 1/0.8 = 1.25
            106: {
                "torque": 1.250000E+03, 
                "torsional_stress": 0.0 # Zero because C=0 in PROD entry
            } 
        }
    },

    # --- CASO 7: Tripode (Full 3D Structure) ---
    {
        "case_name": "Ej7_Tripod_3D",
        "folder": "case_7",
        "filename": "Ej7.dat",
        "expected_displacements": { 
            2: {0: 0.0, 1: 0.0, 2: 0.0},
            # Node 4: T1, T3 (T2 is 0)
            4: {0: 1.860011E-02, 1: 0.0, 2: -8.266668E-03}
        },
        "expected_reactions": { 
            1: {0: -5.555556E+02, 2: 8.333333E+02},
            2: {0:  2.777778E+01, 1: -4.811111E+01},
            3: {0:  2.777778E+01, 1:  4.811111E+01}
        },
        "expected_element_forces": { 
            # Element 171 (Connects 1-4)
            171: {"axial_force": -1.001542E+03, "axial_stress": -2.003084E+02},
            # Element 172 (Connects 2-4)
            172: {"axial_force": -1.001535E+02, "axial_stress": -2.003070E+01},
            # Element 173 (Connects 3-4)
            173: {"axial_force": -1.001535E+02, "axial_stress": -2.003070E+01}
        }
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
            results = solver.solve()
        except Exception as e:
            pytest.fail(f"Solver failed on {case_data['filename']}: {e}")

        # Unpack results
        displacements = results["displacements"]
        reactions = results["reactions"]
        elem_stresses = results["stresses"]

        # 4. VERIFY DISPLACEMENTS
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
        
        # 5. VERIFY REACTIONS
        expected_react = case_data.get("expected_reactions", {})
        
        for nid, dofs_dict in expected_react.items():
            if nid not in reactions:
                # Nastran reports reactions only for constrained nodes.
                pytest.fail(f"Reaction for Node {nid} expected but not found in results.")
            
            # Reactions are stored as a numpy array of 6 floats per node [Fx, Fy, Fz, Mx, My, Mz]
            reac_vec = reactions[nid] 
            
            for local_dof, val_ref in dofs_dict.items():
                val_calc = reac_vec[local_dof]
                
                # Tolerancia: 0.1% o 1e-4 absoluto para valores cercanos a cero
                assert val_calc == pytest.approx(val_ref, rel=1e-3, abs=1e-4), \
                    f"Reaction Mismatch {case_data['case_name']}, Node {nid}, DOF {local_dof}"
        
        # 6. VERIFY ELEMENT STRESSES / FORCES
        expected_forces = case_data.get("expected_element_forces", {})
        
        for eid, quantities in expected_forces.items():
            if eid not in elem_stresses:
                pytest.fail(f"Stresses for Element {eid} expected but not found in results.")
            
            calc_data = elem_stresses[eid] # Dictionary {'axial_force': ..., 'axial_stress': ...}
            
            for qty_name, val_ref in quantities.items():
                if qty_name not in calc_data:
                    pytest.fail(f"Quantity '{qty_name}' not found for Element {eid}")
                
                val_calc = calc_data[qty_name]
                
                assert val_calc == pytest.approx(val_ref, rel=1e-3, abs=1e-4), \
                    f"Stress Mismatch {case_data['case_name']}, Elem {eid}, {qty_name}"

