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
    },

# --- CASO 8: Riel 3D en Eje Y (Ej8.dat) ---
    {
        "case_name": "Ej8_Rail_Y_Axis",
        "folder": "case_8",
        "filename": "Ej8.dat",
        "expected_displacements": {
            8: {1: 1.919646E-04}
        },
        "expected_reactions": {
            7: {0: 3.516000E+01, 1: -3.516000E+01, 2: -2.344000E+01}
        },
        "expected_element_forces": {
            104: {"axial_force": 5.497167E+01, "axial_stress": 5.497167E-01}
        },
        "description": "3D Rod with Node 8 constrained to slide only along Global Y axis."
    },

    # --- CASO 9: Riel 3D en Eje X (Ej9.dat) ---
    {
        "case_name": "Ej9_Rail_X_Axis",
        "folder": "case_9",
        "filename": "Ej9.dat",
        "expected_displacements": {
            8: {0: 7.916667E-03}
        },
        "expected_reactions": {
            7: {0: -1.734220E+02, 1: -3.901995E+02, 2: -2.601330E+02}
        },
        "expected_element_forces": {
            101: {"axial_force": 5.000000E+02,"axial_stress": 5.000000E+00}
        },
        "description": "3D Rod with Node 8 constrained to slide only along Global X axis."
    },

    # --- CASO 10: Cercha Clase DAA (Planar Truss Multi-bay) ---
    {
        "case_name": "Ej10_DAA_Class_Truss",
        "folder": "case_10",
        "filename": "Ej_DAA_clase.dat",
        
        "expected_displacements": {
            1: {0: 0.0, 1: 0.0, 2:0.0},
            2: {0: -4.523520E-03, 1: -8.945180E-03},
            3: {0: -8.994993E-03, 1: -1.793394E-02},
            4: {0: -1.347874E-02, 1: -3.143701E-02},
            5: {0: -2.735212E-03, 1: -4.040562E-02},
            6: {0: 0.0, 1: 0.0, 2:0.0}
            
        },
        "expected_reactions": {
            1: {0: 7.533646E+04, 1: 6.000000E+04, 2:0.0, 3:0.0, 4:0.0, 5:0.0},
            6: {0: -7.533646E+04, 1: 8.000000E+04, 2:0.0}
        },
        
        "expected_element_forces": {
            1: {"axial_force": -8.485281E+04, "axial_stress": -3.367175E+08},
            3: {"axial_force": 2.000000E+04, "axial_stress": 4.494382E+08},
            4: {"axial_force": -6.000000E+04, "axial_stress": -4.477612E+08},
            6: {"axial_force": 4.663537E+03, "axial_stress": 3.129890E+08},
            7: {"axial_force": 1.000000E+05, "axial_stress": 4.484305E+08},
            9: {"axial_force": -1.131371E+05, "axial_stress": -4.489567E+08}
            
        },
        "description": ("Multi-bay planar truss. Validates asymmetric loading and reactions.\n"
        "Ejercicio de clase DAA. Solo ha sido necesario adaptar las SPC de los nodos libres restringiendo todos los girros y movimiento en eje Z.")
    },
    # --- CASO 11: Spider Heterogénea y Barra Torsión (Ej11.dat) ---
    # Reference: ej11.f06
    {
        "case_name": "Ej11_Heterogeneous_Spider",
        "folder": "case_11",
        "filename": "Ej11.dat",
        "expected_displacements": {
            100: {0: 1.366455E-03, 1: -9.998412E-03, 2: 1.794456E-02},
            12: {5: 2.285560E+00} 
        },
        "expected_reactions": {
            1: {0: -4.457096E+02, 1: -8.914191E+02, 2: -6.685643E+02},
            2: {0: 1.843799E+03, 1: -5.900157E+02, 2: -1.696295E+03},
            4: {0: -1.525898E+03, 1: -7.629488E+02, 2: -7.629488E+02, 3: 0.0, 4: 0.0, 5: 0.0},
            11: {0: 0.0, 1: 0.0, 3: -2.000000E+04, 4: -1.000000E+04, 5: -5.000000E+03}
        },
        "expected_element_forces": {
            101: {"axial_force": 1.200110E+03, "axial_stress": 2.400219E+01},
            103: {"axial_force": 4.585915E+03, "axial_stress": 9.171829E+01},
            104: {"axial_force": -1.868835E+03, "axial_stress": -3.737671E+01},
            200: {"axial_force": 0.0, "torque": 2.291288E+04, "axial_stress": 0.0, "torsional_stress": 1.432055E+02}
        },
        "description": "Complex 3D assembly with heterogeneous materials and properties. Validates 3D torque projection."
    },

    # --- CASO 12: Cadena 3D con Desplazamientos Impuestos (Ej12.dat) ---
    # PENDIENTE INPLEMENTAR LECTURA DE DESPLAZAMIENTOS IMPUESTOS (SPC)

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

