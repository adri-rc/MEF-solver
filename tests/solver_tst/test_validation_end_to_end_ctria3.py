"""
@file test_validation_end_to_end_ctria3.py
@brief System tests validation against Nastran reference cases with CTRIA3 elements.
@details
Runs full analysis cycles (Parse -> Assemble -> Solve) on external .dat files
and compares the resulting displacements against reference values extracted 
from MSC Nastran .f06 output files.
"""
import pytest
import os
from pathlib import Path

# Imports de tu proyecto
from src.preprocessing.nastran_reader_v2 import parse_nastran_input
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
   

    #region --- CASO 18: Tripode Hibrido con todos los Stress Points definidos posibles (C-F) ---
    # Descripcion: Variante del caso 17 pero validando tensiones en 4 puntos (C,D,E,F).
    #              Valida la recuperacion de tensiones en esquinas de secciones CBAR.
    {
        "case_name": "Ej18_Hybrid_Tripod_Stress_Recovery_Points",
        "folder": "case_18",
        "filename": "Ej18.dat",
        "description": "Hybrid tripod structure with defined PBAR stress points. Validates detailed stress recovery at cross-section corners (C, D, E, F) and total stress envelopes.",
        "expected_displacements": {
            # Node 3
            3: {
                0:  0.0, # T1
                1:  0.0, # T2
                2:  0.0, # T3
                3:  0.0, # R1
                4:  0.0, # R2
                5:  0.0  # R3
            },
            # Node 4 (Union Central) - Valores extraidos de ej18.f06 (PAGE 11)
            4: {
                0:  1.786003E-04, # T1
                1: -9.267777E-05, # T2
                2:  5.951593E-04, # T3
                3:  3.652847E-05, # R1
                4:  1.530955E-04, # R2
                5:  1.264375E-04  # R3
            },
            # Node 5 (Top CROD)
            5: {
                1:  0.0,           # T1
                2:  2.023731E-03,  # T3 (Elongacion axial dominante)
                5:  1.264375E-04   # R3 
            }
        },
        "expected_reactions": {
            1: {
                0: -2.768665E+02,
                1:  3.374797E+02,
                2: -1.393156E+03,
                3: -2.676242E+03,
                4: -2.106779E+03,
                5: 4.265762E+01
            }, 
            3: {
                0: 6.290351E+01,
                1:  3.530363E+01,
                2: -7.565560E+01,
                3: 1.832815E+02,
                4: -7.837672E-01, 
                5: 1.058018E+02
            }
        },
        "expected_element_forces": {
            # --- ELEMENT 101 (CBAR Pata 1) ---
            101: {
                "axial_force":  1.008341E+03,
                "axial_stress": 2.016681E+01, # AXIAL STRESS column
                "torque":       2.982355E+02,
                
                # Validation Points (SA columns in f06 = Bending Only in this solution)
                # P0 (Point C), P1 (Point D), P2 (Point E), P3 (Point F)
                "stress_points.P0.sigma_bending":  1.195720E+02, # SA1
                "stress_points.P1.sigma_bending": -9.727516E+00, # SA2
                "stress_points.P2.sigma_bending": -1.195720E+02, # SA3
                "stress_points.P3.sigma_bending":  9.727516E+00, # SA4
                
                # Check Total (Axial + Bending) for Point C
                # Total = 20.17 + 119.57 = 139.74
                "stress_points.P0.sigma_total":    1.397388E+02, # Total Point C (=SA-MAX)

                # Check Total (Axial + Bending) for Point D
                # Total = 20.17 + (-9.727) = 10.439
                "stress_points.P1.sigma_total":    1.0439294E+01, # Total Point D

                # Check Total (Axial + Bending) for Point E
                # Total = 20.17 + (-119.572) = -99.402
                "stress_points.P2.sigma_total":    -9.940515E+01, # Total Point E (=SA-MIN)

                # MAX and MIN total stress across all recovery points (SA-MAX, SA-MIN)
                "max_total_stress":  1.397388E+02, # Matches SA-MAX column
                "min_total_stress": -9.940515E+01  # Matches SA-MIN column
            },

            102: {
                "axial_force": 4.694745E+02,
                "axial_stress": 1.564915E+01, # AXIAL STRESS column
                "torque":       3.775480E+01,
                
                # Validation Points (SA columns in f06 = Bending Only in this solution)
                # P0 (Point C), P1 (Point D), P2 (Point E), P3 (Point F)
                "stress_points.P0.sigma_bending":  8.467422E+01, # SA1
                "stress_points.P1.sigma_bending": -8.467422E+01, # SA2
                "stress_points.P2.sigma_bending": -2.327506E+01, # SA3
                "stress_points.P3.sigma_bending":  2.327506E+01, # SA4
                
                # Check Total (Axial + Bending) for Point C
                # Total = 15.64 + 84.67 = 139.74
                "stress_points.P0.sigma_total":    1.003234E+02, # Total Point C (=SA-MAX)

                # Check Total (Axial + Bending) for Point E
                # Total = 15.64 + (-23.27506) = 139.74
                "stress_points.P2.sigma_total":    -7.62591E+00, # Total Point E

                # MAX and MIN total stress across all recovery points (SA-MAX, SA-MIN)
                "max_total_stress":  1.003234E+02, # Matches SA-MAX column
                "min_total_stress": -6.902507E+01  # Matches SA-MIN column
            },
            
            # --- ELEMENT 103 (CBAR Pata 3) ---
            103: {
                "axial_force":  8.132725E+01,
                "axial_stress": 8.132725E+00,
                "torque":       -7.479646E+00,

                # SA1 and SA2
                "stress_points.P0.sigma_bending": -2.214997E+01, # SA1
                "stress_points.P1.sigma_bending":  4.019926E+01, # SA2
                "stress_points.P2.sigma_bending":  0.0, # SA3
                "stress_points.P3.sigma_bending":  0.0, # SA4
                # Check Total (Axial + Bending) for Point C
                # Total = 8.13 + (-22.15) = -14.02
                "stress_points.P0.sigma_total":    -1.401725E+01, # Total Point C (=SA-MIN)

                # MAX and MIN total stress across all recovery points (SA-MAX, SA-MIN)
                "max_total_stress":  4.833198E+01, # Matches SA-MAX
                "min_total_stress": -1.401725E+01  # Matches SA-MIN
            },
            
            # --- ELEMENT 104 (CROD) ---
            104: {
                "axial_force":  2.000000E+03,
                "torque":       7.105427E-15,
                "axial_stress": 1.000000E+02
            }
        }
    },
    #endregion

    #region --- CASO 20: Single CTRIA3 Element (Membrane Tension) ---
    {
        "case_name": "Ej20_Single_CTRIA3",
        "folder": "case_20",
        "filename": "Ej20.dat",
        "description": "Single CTRIA3 element under uniaxial tension. Validates Constant Strain Triangle (CST) formulation.",
        "expected_displacements": {
            2: {0: 2.600000E-05, 1: 0.0, 2: 0.0},
            3: {0: 0.0, 1: 0.0, 2: 0.0}
        },
        "expected_reactions": {
            1: {0: -1.000000E+03, 1: -3.000000E+02},
            3: {1: 3.000000E+02}
        },
        "expected_element_forces": {
            1: {
                "sigma_x": 2.000000E+06,
                "sigma_y": 6.000000E+05,
                "tau_xy": 0.0, 
                "von_mises": 1.777639E+06,
                "major_principal": 2.000000E+06,
                "minor_principal": 6.000000E+05,
                
                "membrane_force_x": 2.000000E+03,
                "membrane_force_y": 6.000000E+02,
                "membrane_force_xy": 0.0,
                "bending_moment_x": 0.0,
                "bending_moment_y": 0.0,
                "bending_moment_xy": 0.0,
                "shear_force_x": 0.0,
                "shear_force_y": 0.0
            }
        }
    },
    #endregion
    #region --- CASO 32: Single CTRIA3 Element (Membrane Tension) ---
    {
        "case_name": "Ej32_3D_Skewed_Membrane",
        "folder": "case_32",
        "filename": "Ej32.dat",
        "description": "3D Skewed CTRIA3 with in-plane edge load. Validates membrane stiffness in rotated coords.",
        "expected_displacements": {
            2: {
                0: 1.323787E-05, # T1
                1: 2.329880E-05, # T2
                2: 0.0           # T3 (Debe ser 0 en membrana pura local si la carga está bien alineada)
            },
            # Node 3 (Restringido en Z en este caso según F06)
            3: {
                0: -5.701283E-07, # T1
                1:  1.547821E-05, # T2
                2:  0.0           # T3
            }
        },
        "expected_reactions": {
            1: {
                0: -2.000000E+03, # T1 (FX react)
                1: -5.000000E+03, # T2 (FY react)
                2: -1.000000E+03  # T3 (FZ react)
            }
        },
        "expected_element_forces": {
            1: {
                "membrane_force_x":  2.621112E+03,
                "membrane_force_y":  1.352662E-06, # Casi cero
                "membrane_force_xy": -4.920298E-08, # Casi cero
                
                # Momentos deben ser despreciables (ruido numérico < 1e-5)
                "bending_moment_x":  2.679009E-06 
            }
        }
    },
    #endregion

    #region --- CASO 33: 3D Skewed - 'Pure Bending' Load Configuration ---
    # Descripcion: Carga normal al plano en Nodo 3.
    # NOTA: Según el output F06 (Page 10), el nodo 3 tiene T3=0.0, lo que indica
    # que está restringido en Z. Esto genera fuerzas de membrana significativas
    # (FX, FY) para reaccionar a la carga inclinada, en lugar de flexión pura de voladizo.
    # Se valida contra los resultados exactos del F06 proporcionado.
    {
        "case_name": "Ej33_3D_Skewed_Bending_Mixed",
        "folder": "case_33",
        "filename": "Ej33.dat",
        "description": "3D Skewed CTRIA3 with normal load. Validates mixed coupling/reactions.",
        "expected_displacements": {
            2: {
                0: 0.0, 
                1: 0.0, 
                2: 0.0
            },
            # Node 3 (Cargado, T3=0 en output)
            3: {
                0: 6.921123E-06, # T1
                1: 1.528502E-06, # T2
                2: 0.0           # T3
            }
        },
        "expected_reactions": {
            # Node 1 - Page 11
            1: {
                0: -1.095923E+03, # T1
                1: -1.398952E+02, # T2
                2: -8.368408E+02,  # T3
                3:  1.358457E-01, # R1
                4: -1.941154E-02, # R2
                5: -1.746568E-01 # R3
            },
            2: {
                0: -3.040769E+02, # T1
                1: 3.398952E+02, # T2
                2: -2.742704E+02,  # T3
                3: 0.0, # R1
                4: 0.0, # R2
                5: 0.0  # R3
            },
            3: {
                2: 2.911111E+03  # T3
            }
        },
        "expected_element_forces": {
            # Element 1 - Page 12 (Resultados dominados por Membrana debido a constraints)
            1: {
                "membrane_force_x":  2.070368E+02,
                "membrane_force_y":  6.273842E+02,
                "membrane_force_xy": 1.940644E+02,
                
                # Momentos muy bajos
                #"bending_moment_x":  7.299514E-07
            }
        }
    },
    #endregion

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
        
        if case_data['case_name'] in ['Ej21_Multi_CTRIA3_Patch', 'Ej22_Bending_Simple', 'Ej23_Mixed_Materials']:
            print(f"\nCAPTURING RESULTS FOR {case_data['case_name']}:")
            print(f"DISPLACEMENTS TYPE: {type(displacements)}")
            print(f"DISPLACEMENTS: {displacements}")
            
            # If it is a dict, iterate. If array, just printed.
            if isinstance(displacements, dict):
                for nid, disps in displacements.items():
                    print(f"Node {nid}: {disps}")
            
            print("ELEMENT STRESSES:")
            # stresses might be a dict of dicts
            print(elem_stresses)

        # 4. COMPARE RESULTS
        rel_tol = case_data.get("rel_tol", 1e-3)
        abs_tol_disp = case_data.get("abs_tol_disp", 1e-9)
        abs_tol_force = case_data.get("abs_tol_force", 1e-4)

        # 4. VERIFY DISPLACEMENTS

        expected_disp = case_data.get("expected_displacements", {})
        
        for nid, dofs_dict in expected_disp.items():
            # Get Global Indices
            try:
                global_indices = solver.dof_mapper.get_dof_indices(nid)
            except ValueError:
                pytest.fail(f"Node {nid} expected in results but not found in Solver Mapper.")
            
            for local_dof, val_ref in dofs_dict.items():
                idx = global_indices[local_dof]
                val_calc = displacements[idx]
                
                # Check with tolerance (0.1% relative error is generous but safe for FEM)
                # Absolute tolerance 1e-9 handles near-zero values.
                assert val_calc == pytest.approx(val_ref, rel=rel_tol, abs=abs_tol_disp), \
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
                assert val_calc == pytest.approx(val_ref, rel=rel_tol, abs=abs_tol_force), \
                    f"Reaction Mismatch {case_data['case_name']}, Node {nid}, DOF {local_dof}"

        
        # 6. VERIFY ELEMENT STRESSES / FORCES
        expected_forces = case_data.get("expected_element_forces", {})
        
        for eid, quantities in expected_forces.items():
            if eid not in elem_stresses:
                pytest.fail(f"Stresses for Element {eid} expected but not found in results.")
            
            calc_data = elem_stresses[eid] # Dictionary {'axial_force': ..., 'axial_stress': ...}
            
            for qty_name, val_ref in quantities.items():
             
                #--- UPDATE FOR NESTED DICTIONARIES (CBAR) ---
                val_calc = calc_data
                try:
                    if "." in qty_name:
                        keys = qty_name.split(".")
                        for k in keys:
                            # Verificamos si es indexable (dict o list) antes de acceder
                            if not hasattr(val_calc, '__getitem__'):
                                raise ValueError(f"Cannot access key '{k}' in scalar value: {val_calc}")
                            val_calc = val_calc[k]
                    else:
                        val_calc = calc_data[qty_name]
                except (KeyError, IndexError, ValueError, TypeError) as e:
                    pytest.fail(f"Error accessing '{qty_name}' in Element {eid}. Details: {e}")

                # Verificamos que sea un número (float o numpy scalar)
                assert val_calc == pytest.approx(val_ref, rel=rel_tol, abs=abs_tol_force), \
                    f"Stress Mismatch {case_data['case_name']}, Elem {eid}, {qty_name}"

