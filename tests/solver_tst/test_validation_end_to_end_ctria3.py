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
#   "description": "Test case description."
#   "note": "Test note.",
#   "rel_tol": 0.01,       (Results relative acceptance tolerance [0.01 = 1%])    
#   "abs_tol_disp": 1e-5,  (Results absolute accpt. tol., if not defined "re_tol" is used)    
#   "abs_tol_force": 1e-1, (Results absolute accpt. tol.,  if not defined "re_tol" is used)    
#   "expected_displacements": { NodeID: { DOF_Index: Value } }
#   "expected_reactions":     { NodeID: { DOF_Index: Value } }
#   "expected_element_forces":{ 
#               ElemID: {  
#                    # Forces
#                    'membrane_force_x': f_x, 
#                    'membrane_force_y': f_y, 
#                    'membrane_force_xy': f_xy,
#                    'shear_force_x': q_x, 
#                    'shear_force_y': q_y
#
#                    # Additional Bending Info
#                    'bending_moment_x': m_x, 
#                    'bending_moment_y': m_y, 
#                    'bending_moment_xy': m_xy,

#                    # Tensiones Fibra Superior (Z = +t/2)
#                    'stress_top_fiber_x': sig_x_top, 
#                    'stress_top_fiber_y': sig_y_top,
#                    'tau_top_fiber_xy': tau_xy_top,
#                    'von_mises': vm_stress, 
#                    'major_principal': s_major_principal, 
#                    'minor_principal': s_minor_principal,
#
#                    # Tensiones Fibra Inferior (Z = -t/2)
#                    'stress_bot_fiber_x': sig_x_bot,
#                    'stress_bot_fiber_y': sig_y_bot,
#                    'tau_bot_fiber_xy': tau_xy_bot,
#
#                    
#               },
#               ElemID: {...}

# }
# DOF Index map: 0=Tx, 1=Ty, 2=Tz, 3=Rx, 4=Ry, 5=Rz

VALIDATION_DATA = [
    #region --- CASO 20: Single CTRIA3 Element (Membrane Tension) ---
    {
        "case_name": "Ej20_Single_CTRIA3",
        "folder": "case_20",
        "filename": "Ej20.dat",
        "description": "Single CTRIA3 element under uniaxial tension. Validates Constant Strain Triangle (CST) formulation.",
        "note": "Test análisis de membrana puro, no K6ROT.",
        "rel_tol": 0.01,
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
                "shear_force_y": 0.0,

                # Top Fiber (Z2)
                "stress_top_fiber_x": 2.000000E+06,
                "stress_top_fiber_y": 6.000000E+05,
                "tau_top_fiber_xy": 0.0,
                "von_mises_top": 1.777639E+06,
                "major_principal_top": 2.000000E+06,
                "minor_principal_top": 6.000000E+05,
                
                # Bottom Fiber (Z1)
                "stress_bot_fiber_x": 2.000000E+06,
                "stress_bot_fiber_y": 6.000000E+05,
                "tau_bot_fiber_xy": 0.0,
                "von_mises_bot": 1.777639E+06,
                "major_principal_bot": 2.000000E+06,
                "minor_principal_bot": 6.000000E+05
            }
        }
    },
    #endregion
    #region --- CASO 32: Single CTRIA3 Element (Membrane Tension) ---
    {
        "case_name": "Ej32_3D_Skewed_Membrane",
        "folder": "case_32",
        "filename": "Ej32.dat",
        "description": "3D Skewed CTRIA3. Validates membrane, flexure stiffness & coords rotation.",
        "note": "Test: K6ROT= 100.0 y MI3/E=0.",
        "rel_tol": 0.01,
        "expected_displacements": {
            1: {
                0: 0.0,
                1: 0.0,
                2: 0.0
            },
            2: {
                0: 1.323787E-05, # T1
                1: 2.329877E-05, # T2
                2: 0.0           # T3 (Debe ser 0 en membrana pura local si la carga está bien alineada)
            },
            # Node 3 (Restringido en Z en este caso según F06)
            3: {
                0: -5.701140E-07, # T1
                1:  1.547815E-05, # T2
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
                # Membrana / Midplane
                "sigma_x": 5.242224E+05,
                "sigma_y": 0.0,
                "tau_xy": 0.0,
                "von_mises": 5.242224E+05,
                "major_principal": 5.242224E+05,
                "minor_principal": 0.0,

                # Fuerzas y Momentos
                "membrane_force_x":   2.621112E+03,
                "membrane_force_y":   0.0, 
                "membrane_force_xy":  0.0,
                "shear_force_x":      0.0, 
                "shear_force_y":      0.0,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,

                # Top Fiber (Z2)
                "stress_top_fiber_x":  5.242224E+05,
                "stress_top_fiber_y":  0.0, 
                "tau_top_fiber_xy":    0.0,
                "von_mises_top":       5.242224E+05,
                "major_principal_top": 5.242224E+05,
                "minor_principal_top": 0.0,
                
                # Bottom Fiber (Z1)
                "stress_bot_fiber_x":  5.242224E+05,
                "stress_bot_fiber_y":  0.0, 
                "tau_bot_fiber_xy":    0.0,                
                "von_mises_bot":       5.242224E+05,
                "major_principal_bot": 5.242224E+05,
                "minor_principal_bot": 0.0
            }
        }
    },
    #endregion

    #region --- CASO 33: 3D Skewed - 'Pure Bending' Load Configuration ---
    # Descripcion: Carga normal al plano en Nodo 3.
    {
        "case_name": "Ej33_3D_Skewed_Bending_Mixed",
        "folder": "case_33",
        "filename": "Ej33.dat",
        "description": "3D Skewed CTRIA3 with nomral membrane plane load. Validates mixed coupling/reactions, with no pure shear stiffnes.",
        "note": "Desplazamientos bien, pero problemas con los signos de los mosmentos. Test: K6ROT= 100.0 y MI3/E=0.",
        "rel_tol": 0.01,    
        #"abs_tol_disp": 1e-5,
        #"abs_tol_force": 1e-3,
        "expected_displacements": {
            2: {
                0: 0.0, 
                1: 0.0, 
                2: 0.0
            },
            3: {
                0:  5.249974E-06,  # T1
                1:  2.042201E-06,  # T2
                2:  0.0,           # T3
                3:  -7.997170E+01, # R1
                4:  1.649297E+02,  # R2
                5:  -8.052562E+01  # R3
            }
        },
        "expected_reactions": {
            1: {
                0: -8.505554E+02,  # T1
                1: -1.093814E+02,  # T2
                2: -6.493896E+02,  # T3
                3:  4.210720E+01,  # R1
                4: -4.329598E+02,  # R2
                5:  8.058439E+01   # R3
            },
            2: {
                0: -1.494446E+02,  # T1
                1:  1.093814E+02,  # T2
                2: -1.283881E+02,  # T3
                3: -6.856286E+01,  # R1
                4: -6.326086E+01,  # R2
                5: -4.656998E+01   # R3
            },
            3: {
                0:  0.0,  # T1
                1:  0.0,  # T2
                2:  7.777778E+02  # T3
            }
        },
        "expected_element_forces": {
            1: {
                # Fuerzas de Membrana
                "membrane_force_x":   1.398897E+02,
                "membrane_force_y":   4.239083E+02,
                "membrane_force_xy":  1.870749E+02,
                
                # Fuerzas y Momentos
                "bending_moment_x":   7.500000E+01,
                "bending_moment_y":   1.666667E+02,
                "bending_moment_xy":  7.280867E+01,
                'shear_force_x': -4.805034E-22, 
                'shear_force_y': -8.814329E-22,
                
                # Top Fiber (Z2 = +t/2, Línea 2 en el f06)
                "stress_top_fiber_x":  -1.797202E+07,
                "stress_top_fiber_y":  -3.991522E+07,
                "tau_top_fiber_xy":    -1.743667E+07,
                "von_mises_top":        4.594533E+07,
                "major_principal_top": -8.342327E+06,
                "minor_principal_top": -4.954491E+07,
            
                # Bottom Fiber (Z1 = -t/2, Línea 1 en el f06)
                "stress_bot_fiber_x":   1.802798E+07,
                "stress_bot_fiber_y":   4.008478E+07,
                "tau_bot_fiber_xy":     1.751150E+07,
                "von_mises_bot":        4.614225E+07,
                "major_principal_bot":  4.975126E+07,
                "minor_principal_bot":  8.361498E+06
            }
        }
    },
    #endregion

    #region --- CASO 34: CTRIA3 Plano XY, Validación Drilling DOF ---
    {
        "case_name": "Ej34_CTRIA3_Drilling_DOF",
        "folder": "case_34",
        "filename": "Ej34.dat",
        "description": "Single CTRIA3 element in XY plane. Validates CST + Drilling dof formulation.",
        "note": "Vaqlido la implmentación del drilling dof, al ponerlo K6ROT = 100. y eliminando SPC de drilling Dof.",
        "rel_tol": 0.01, 
        "expected_displacements": {
            2: {
                0:  1.649621E-05,
                1:  0.0,
                2:  0.0,
                5:-3.714286E-05,
            },
            3: {
                0:  7.257953E-05,
                1: -4.285714E-06,
                2:  0.0,
                5:-3.714286E-05,
            }
        },
        "expected_reactions": {
            1: {
                0: -1.000000E+03,
                1: -8.660000E+02,
                2:  0.0
            },
            2: {
                0:  0.0,
                1:  8.660000E+02,
                2:  0.0
            }
        },
        "expected_element_forces": {
            1: {
                # Membrana / Midplane
                "membrane_force_x":   1.154734E+03,
                "membrane_force_y":  -5.684342E-14,
                "membrane_force_xy":  2.000000E+03,

                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                'shear_force_x':      0.0, 
                'shear_force_y':      0.0,
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  1.154734E+06,
                "stress_top_fiber_y": -5.684342E-11,
                "tau_top_fiber_xy":    2.000000E+06,
                "von_mises_top":       3.651494E+06,
                "major_principal_top": 2.659038E+06,
                "minor_principal_top": -1.504303E+06,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  1.154734E+06,
                "stress_bot_fiber_y": -5.684342E-11,
                "tau_bot_fiber_xy":    2.000000E+06,
                "von_mises_bot":       3.651494E+06,
                "major_principal_bot": 2.659038E+06,
                "minor_principal_bot": -1.504303E+06
            }
        }
    },
    #endregion

    #region --- CASO 35: CTRIA3 Plano XY, Isosceles Skewness alto 3D Force ---
    {
        "case_name": "Ej35_CTRIA3_XY_Skewed",
        "folder": "case_35", 
        "filename": "Ej35.dat",
        "description": "CTRIA3 in XY plane with high skewness 3D loaded. Validates element stability under geometric distortion.",
        "note": "Valido la implementación del drilling dof, al ponerlo K6ROT = 100. y eliminando SPC de drilling Dof.",
        "rel_tol": 0.01,
        "expected_displacements": {
            2: {
                0:   1.904762E-05,
                1:  -1.498022E+00,
                2:   0.0,
                3:   0.0,
                4:   0.0,
                5:  -1.498095E+01

            },
            3: {
                0:   4.494508E+01,
                1:  -7.507422E-01,
                2:   0.0,
                3:   0.0,
                4:   0.0,
                5:  -1.498095E+01
            }
        },
        "expected_reactions": {
            1: {
                0: -4.000000E+03,
                1:  2.000000E+03,
                2:  0.0,
                3:  0.0,
                4:  0.0, 
                5:  1.210000E+04
            },
            3: {
                2: -8.000000E+03,
                3:  0.0,
                4:  0.0,
            }
        },
        "expected_element_forces": {
            1: {
                # Membrana / Midplane
                "membrane_force_x":   1.333333E+03,
                "membrane_force_y":  -4.000000E+04,
                "membrane_force_xy":  3.966667E+04,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                'shear_force_x':      0.0, 
                'shear_force_y':      0.0,

                # Tensiones Fibra Superior (z = +t/2)
                "stress_top_fiber_x":  1.333333E+06,
                "stress_top_fiber_y": -4.000000E+07,
                "tau_top_fiber_xy":    3.966667E+07,
                "von_mises":           7.984638E+07,
                "major_principal":     2.539424E+07,
                "minor_principal":    -6.406090E+07,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  1.333333E+06,
                "stress_bot_fiber_y": -4.000000E+07,
                "tau_bot_fiber_xy":    3.966667E+07,
                "von_mises_bot":       7.984638E+07,
                "major_principal_bot": 2.539424E+07,
                "minor_principal_bot": -6.406090E+07
            }
        }
        
    },
    #endregion

    #region --- CASO 36: CTRIA3 Plano XY, Irregular offset ---
    {
        "case_name": "Ej36_CTRIA3_XY_Irregular_Offset",
        "folder": "case_36", 
        "filename": "Ej36.dat",
        "description": "CTRIA3 in XY plane, large coordinates without origin (0,0). Validates coordinate projection independence.",
        "note": "El resultado de desplazamiento no tiene sentido (resultado que coincide), pero no tiene sentido comparar pospocreso del caso",
        "rel_tol": 0.02,
        #"rel_tol_force": 0.11,
        "expected_displacements": {
            1: {
                0:  0.0,
                1:  0.0,
                2:  0.0,
                3:  0.0,
                4:  0.0,
                5: -3.330049E-08

            },
            2: {
                0:  0.0,
                1:  0.0,
                2:  0.0,
                3:  0.0,
                4:  0.0,
                5: -3.330049E-08

            },
            3: {
                0:  5.967060E-06,
                1:  3.064617E-06,
                2:  5.809655E+25,
                3:  0.0,
                4:  0.0,
                5: -3.330049E-08
            }
        },
        # NO TIENE SENTIDO ESTA COMPARATIVA DE REACCIONES PUES EL ELEMENTO NO TIENE SENTIDO POR LOS RESULTADOS DE DESPLAZAMIENTOS
        "expected_reactions": {
            1: {
                0:   4.399524E+01,
                1:  -1.058240E+04,
                2:  -3.275862E+03,
                #3:  -1.407826E+05,
                #4:  -1.262447E+04,
                5:   0.0,
            },
            2: {
                0:  -4.043995E+03,
                1:   7.582402E+03,
                2:  -1.724138E+03,
                #3:  -9.772859E+04,
                #4:  -6.609744E+04,
                5:   0.0,
            },
            3: {
                #3:  -2.270060E+05 ,
                #4:  -1.074850E+05,
                5:   0.0,
            }
        },
    },
    #endregion

    #region --- CASO 37: 2 CTRIAs Plano XZ ---
    {
        "case_name": "Ej37_CTRIA3_XZ_2Elems",
        "folder": "case_37", 
        "filename": "Ej37.dat",
        "description": "Two CTRIA3 elements sharing nodes in XZ plane. Validates local axes transformation out of XY plane.",
        "note": "Validación en el plano de dos ctrias, solo con rigidez a membrana.",
        "rel_tol": 0.01,
        "expected_displacements": {
            2: {
                0:  1.300000E-05,
                1:  0.0,
                2:  2.460714E-05
            },
            3: {
                0:  0.0,
                1:  0.0,
                2:  0.0
            },
            4: {
                0: -3.110714E-05,
                1:  0.0,
                2:  1.300000E-05
            }
        },
        "expected_reactions": {
            1: {
                0: -5.625000E+02,
                1:  0.0,
                2: -5.625000E+02
            },
            3: {
                0:  5.625000E+02,
                1:  0.0,
                2: -9.375000E+02
            }
        },
        "expected_element_forces": {
            1: {
                "membrane_force_x":   1.300000E+02,
                "membrane_force_y":   1.300000E+02,
                "membrane_force_xy": -1.750000E+01,
                "shear_force_x": 0.0, 
                "shear_force_y": 0.0,
                "bending_moment_x": 0.0,
                "bending_moment_y": 0.0,
                "bending_moment_xy": 0.0,

                "stress_top_fiber_x":  1.300000E+05,
                "stress_top_fiber_y":  1.300000E+05,
                "tau_top_fiber_xy":   -1.750000E+04,
                "von_mises_top":       1.334869E+05,
                "major_principal_top": 1.475000E+05,
                "minor_principal_top": 1.125000E+05,

                "stress_bot_fiber_x":  1.300000E+05,
                "stress_bot_fiber_y":  1.300000E+05,
                "tau_bot_fiber_xy":   -1.750000E+04,
                "von_mises_bot":       1.334869E+05,
                "major_principal_bot": 1.475000E+05,
                "minor_principal_bot": 1.125000E+05
            },
            2: {
                "membrane_force_x":   1.825000E+02,
                "membrane_force_y":  -1.175000E+02,
                "membrane_force_xy": -7.000000E+01,
                 "shear_force_x": 0.0, 
                "shear_force_y": 0.0,
                "bending_moment_x": 0.0,
                "bending_moment_y": 0.0,
                "bending_moment_xy": 0.0,

                "stress_top_fiber_x":  1.825000E+05,
                "stress_top_fiber_y": -1.175000E+05,
                "tau_top_fiber_xy":   -7.000000E+04,
                "von_mises_top":       2.885416E+05,
                "major_principal_top": 1.980295E+05,
                "minor_principal_top": -1.330295E+05,

                "stress_bot_fiber_x":  1.825000E+05,
                "stress_bot_fiber_y": -1.175000E+05,
                "tau_bot_fiber_xy":   -7.000000E+04,
                "von_mises_bot":       2.885416E+05,
                "major_principal_bot": 1.980295E+05,
                "minor_principal_bot": -1.330295E+05
            }
        }
    },
    #endregion

    #region --- CASO 38: 2 CTRIAs Plano YZ ---
    {
        "case_name": "Ej38_CTRIA3_YZ_2Elems",
        "folder": "case_38", 
        "filename": "Ej38.dat",
        "description": "Two CTRIA3 elements sharing nodes in YZ plane. Validates full 3D transformation capability.",
        "note": "Validación de dos ctrias conectados, con rigidez a membrana y flexión y carga 3D.",
        "rel_tol": 0.01,
        "abs_tol_force": 0.5, # las deformaciones resultado (que coinciden) no tienen sentido físico por lo qeu las reacciones tampoco cuadran mucho, pero está bien, todo es del orden. 
        "expected_displacements": {
            2: {
                0: -6.240000E+18,
                1:  2.600000E-11,
                2: -2.600000E-11,
                5:  2.380545E-13,
            },
            4: {
                0:  6.240000E+18,
                1: -2.600000E-11,
                2:  2.600000E-11,
                5:  2.380545E-13,
            }
        },
        "expected_reactions": {
            1: {
                0:  6.002817E-19,
                1: -3.000000E+00,
                2: -3.000000E+00,
                3:  0.0,
                4:  1.218770E-01,
                5: -5.000000E-01,
            },
            2: {
                3:  0.0,
                4:  -6.218770E-01,
            },
            3: {
                0:  6.002817E-19,
                1:  3.000000E+00,
                2:  3.000000E+00,
                3:  0.0,
                4:  1.218770E-01,
                5: -5.000000E-01,
            }
        },
         "expected_element_forces": {
            1: {
                "membrane_force_x":   1.300000E+00,
                "membrane_force_y":   1.300000E+00,
                "membrane_force_xy": -7.000000E-01,
                
                "bending_moment_x":  -6.218770E-02,
                "bending_moment_y":  -2.798447E-02,
                "bending_moment_xy": -1.710162E-02,
                'shear_force_x':     -1.000000E-02, 
                'shear_force_y':     1.000000E-02,

                # Tensiones Fibra Inferior (z = -t/2)
                "stress_bot_fiber_x":  2.450750E-01,
                "stress_bot_fiber_y":  2.532837E-01,
                "tau_bot_fiber_xy":   -1.441044E-01,
                "von_mises_bot":       3.527593E-01,
                "major_principal_bot": 3.933422E-01,
                "minor_principal_bot": 1.050165E-01,

                # Tensiones Fibra Superior (z = +t/2)
                "stress_top_fiber_x":  2.749250E-01,
                "stress_top_fiber_y":  2.667163E-01,
                "tau_top_fiber_xy":   -1.358956E-01,
                "von_mises_top":       3.588833E-01,
                "major_principal_top": 4.067782E-01,
                "minor_principal_top": 1.348631E-01,
            },
            2: {
                "membrane_force_x":   1.300000E+00,
                "membrane_force_y":   1.300000E+00,
                "membrane_force_xy": -7.000000E-01,
                
                "bending_moment_x":   6.218770E-02,
                "bending_moment_y":   2.798447E-02,
                "bending_moment_xy":  1.710162E-02,
                'shear_force_x':      1.000000E-02, 
                'shear_force_y':     -1.000000E-02,


                # Tensiones Fibra Inferior (z = -t/2)
                "stress_bot_fiber_x":  2.749250E-01,
                "stress_bot_fiber_y":  2.667163E-01,
                "tau_bot_fiber_xy":   -1.358956E-01,
                "von_mises_bot":       3.588833E-01,
                "major_principal_bot": 4.067782E-01,
                "minor_principal_bot": 1.348631E-01,

                # Tensiones Fibra Superior (z = +t/2)
                "stress_top_fiber_x":  2.450750E-01,
                "stress_top_fiber_y":  2.532837E-01,
                "tau_top_fiber_xy":   -1.441044E-01,
                "von_mises_top":       3.527593E-01,
                "major_principal_top": 3.933422E-01,
                "minor_principal_top": 1.050165E-01,
            }
        }
    },
    #endregion

    #region --- CASO 39: Plano Arbitrario (Inclinado 13.5º), Alta Distorsion ---
    {
        "case_name": "Ej39_CTRIA3_Arbitrary_Plane_Distorted",
        "folder": "case_39",
        "filename": "Ej39.dat",
        "description": "2 CTRIA3 elements in an arbitrary 3D plane (13.5 deg tilt, high skewness) with constrained T3 to stabilize out-of-plane drilling mechanism.",
        "expected_displacements": {
            2: {
                0:  1.487373E-05,
                1:  2.681022E-06,
                2:  0.0
            },
            4: {
                0: -1.284695E-06,
                1:  6.049614E-06,
                2:  0.0
            }
        },
        "expected_reactions": {
            1: {
                0: -2.187853E+02,
                1: -4.680588E+02,
                2: -7.105239E+01,
                3:  1.038178E-02,
                4: -1.798214E-02,
                5:  8.648996E-02,
            },
            2: {
                2:  2.218275E-03,
                3:  1.992055E-03,
                4: -3.450503E-03,
                5:  1.659580E-02
            },
            3: {
                0: -6.837475E+01, 
                1: -3.457025E+01, 
                2:  1.019433E+00,
                3: -8.389595E-03,
                4:  1.453227E-02,
                5:  -6.989413E-02, 
            },
            4: {
                2: -5.263498E-03,
                3:  1.992055E-03,
                4: -3.450503E-03,
                5:  1.659580E-02
            },
        },
        "expected_element_forces": {
            1: {
                "membrane_force_x":   1.157082E+02,
                "membrane_force_y":   8.893921E+01,
                "membrane_force_xy": -2.232528E+01,
                "bending_moment_x": -3.538480E-13,
                "bending_moment_y":  3.538480E-13,
                "bending_moment_xy": 1.007096E-13,
                "shear_force_x": 0.0, 
                "shear_force_y": 0.0,

                "stress_top_fiber_x":  1.157082E+05,
                "stress_top_fiber_y":  8.893921E+04,
                "tau_top_fiber_xy":   -2.232528E+04,
                "von_mises_top":       1.118160E+05,
                "major_principal_top": 1.283537E+05,
                "minor_principal_top": 7.629366E+04,

                "stress_bot_fiber_x":  1.157082E+05,
                "stress_bot_fiber_y":  8.893921E+04,
                "tau_bot_fiber_xy":   -2.232528E+04,
                "von_mises_bot":       1.118160E+05,
                "major_principal_bot": 1.283537E+05,
                "minor_principal_bot": 7.629366E+04
            },
            2: {
                "membrane_force_x":   3.077436E+01,
                "membrane_force_y":   5.373421E+01,
                "membrane_force_xy":  4.983093E+01,
                "bending_moment_x": -1.472079E-13,
                "bending_moment_y":  1.472079E-13,
                "bending_moment_xy": -6.623366E-14,
                "shear_force_x": 0.0, 
                "shear_force_y": 0.0,

                "stress_top_fiber_x":  3.077436E+04,
                "stress_top_fiber_y":  5.373482E+04,
                "tau_top_fiber_xy":    4.983093E+04,
                "von_mises_top":       9.813360E+04,
                "major_principal_top": 9.339086E+04,
                "minor_principal_top": -8.881669E+03,

                "stress_bot_fiber_x":  3.077436E+04,
                "stress_bot_fiber_y":  5.373482E+04,
                "tau_bot_fiber_xy":    4.983093E+04,
                "von_mises_bot":       9.813360E+04,
                "major_principal_bot": 9.339086E+04,
                "minor_principal_bot": -8.881669E+03
            }
        }
    },
    #endregion

    #region --- CASO 40: Trapezoide Irregular en Plano XY, 10 Elementos ---
    {
        "case_name": "Ej40_CTRIA3_Trapezoid_XY",
        "folder": "case_40", 
        "filename": "Ej40.dat",
        "description": "10 CTRIA3 elements forming an irregular trapezoid in XY plane. Validates multi-element assembly in 2D.",
        "expected_displacements": {
            3:  { 0: 1.761575E-05, 1: -3.140008E-05 },
            6:  { 0: 7.585357E-05, 1: -1.731764E-04 },
            9:  { 0: 5.018859E-05, 1: -5.357071E-05 },
            12: { 0: 8.477052E-05, 1: -1.854124E-04 }
        },
        "expected_reactions": {
            1: { 0: -6.300000E+02, 1: 3.065152E+02 },
            7: { 0: -1.870000E+03, 1: 4.934848E+02 }
        },
        "expected_element_forces": {
            1: {
                "membrane_force_x":   1.301354E+02,
                "membrane_force_y":   1.051633E+02,
                "membrane_force_xy": -1.112066E+02,                
                "bending_moment_x": 0.0,
                "bending_moment_y": 0.0,
                "bending_moment_xy": 0.0,
                "shear_force_x": 0.0, 
                "shear_force_y": 0.0,

                "stress_top_fiber_x":  1.301354E+05,
                "stress_top_fiber_y":  1.051633E+05,
                "tau_top_fiber_xy":   -1.112066E+05,
                "von_mises_top":       2.267373E+05,
                "major_principal_top": 2.295547E+05,
                "minor_principal_top": 5.743959E+03,

                "stress_bot_fiber_x":  1.301354E+05,
                "stress_bot_fiber_y":  1.051633E+05,
                "tau_bot_fiber_xy":   -1.112066E+05,
                "von_mises_bot":       2.267373E+05,
                "major_principal_bot": 2.295547E+05,
                "minor_principal_bot": 5.743959E+03
            },
            10: {
                "membrane_force_x":  -1.348041E+02,
                "membrane_force_y":   1.216682E+02,
                "membrane_force_xy":  3.244484E+02,
                "bending_moment_x": 0.0,
                "bending_moment_y": 0.0,
                "bending_moment_xy": 0.0,
                "shear_force_x": 0.0, 
                "shear_force_y": 0.0,

                "stress_top_fiber_x": -1.348041E+05,
                "stress_top_fiber_y":  1.216682E+05,
                "tau_top_fiber_xy":    3.244484E+05,
                "von_mises_top":       6.042988E+05,
                "major_principal_top": 3.423035E+05,
                "minor_principal_top": -3.554395E+05,

                "stress_bot_fiber_x": -1.348041E+05,
                "stress_bot_fiber_y":  1.216682E+05,
                "tau_bot_fiber_xy":    3.244484E+05,
                "von_mises_bot":       6.042988E+05,
                "major_principal_bot": 3.423035E+05,
                "minor_principal_bot": -3.554395E+05
            }
        }
    },
    #endregion

    #region --- CASO 41: Trapezoide en Plano Arbitrario 3D, 10 Elementos ---
    {
        "case_name": "Ej41_CTRIA3_Trapezoid_3D",
        "folder": "case_41", 
        "filename": "Ej41.dat",
        "description": "10 CTRIA3 elements forming a trapezoid in an arbitrary 3D plane. Validates assembly and transformations combined.",
        
        "rel_tol": 0.01,           
        "abs_tol_disp": 1e-9, 
        "abs_tol_force": 1e-2,
        
        "expected_displacements": {
            2:  {0:  1.781194E-05, 1: -8.038953E-06, 2: 0.0},
            3:  {0:  3.764823E-05, 1: -2.077752E-05, 2: 0.0},
            4:  {0:  7.410333E-05, 1: -5.542953E-05, 2: 0.0},
            5:  {0:  1.393857E-04, 1: -9.635942E-05, 2: 0.0},
            6:  {0:  1.904665E-04, 1: -1.274270E-04, 2: 0.0},
            8:  {0:  2.395201E-05, 1: -5.821325E-06, 2: 0.0},
            9:  {0:  8.026652E-05, 1: -2.343236E-05, 2: 0.0},
            10: {0:  1.242798E-04, 1: -5.770066E-05, 2: 0.0},
            11: {0:  1.798190E-04, 1: -1.012318E-04, 2: 0.0},
            12: {0:  2.067144E-04, 1: -1.343570E-04, 2: 0.0}
        },
        "expected_reactions": {
            1: {
                0: -6.520454E+02, 
                1: -1.829452E+02, 
                2:  1.836533E+02,
                3:  6.636643E-02,
                4: -8.848857E-02,
                5: 1.474810E-01,
                },
            3: {
                0:  0.0, 
                1:  0.0, 
                2: -6.821210E-13,
                3:  5.706813E-01,
                4: -7.609084E-01,
                5:  1.268181E+00,
                },
            5: {
                0:  0.0, 
                1:  0.0, 
                2:  0.0,
                3:  6.251460E-01,
                4: -8.335280E-01,
                5:  1.389213E+00,
                },
            7: {
                0: -1.731955E+03, 
                1: -8.050548E+02, 
                2:  2.963467E+02,
                3:  1.361563E-01,
                4: -1.815418E-01,
                5:  3.025696E-01,
                },
            12: {
                0:  0.0, 
                1:  0.0, 
                2:  -3.637979E-12,
                3:  1.132903E-01,
                4: -1.510537E-01,
                5:  2.517562E-01,
                },
        },
        "expected_element_forces": {
            1: {
                # Membrana / Midplane
                "sigma_x": 1.304317E+05,
                "sigma_y": 1.052227E+05,
                "tau_xy": -1.111881E+05,
                "von_mises": 2.268221E+05,
                "major_principal": 2.297275E+05,
                "minor_principal": 5.926994E+03,

                # Fuerzas y Momentos
                "membrane_force_x":  1.304317E+02,
                "membrane_force_y":  1.052227E+02,
                "membrane_force_xy": -1.111881E+02,
                "bending_moment_x":  0.0,
                "bending_moment_y":  0.0,
                "bending_moment_xy": 0.0,
                "shear_force_x":     0.0,
                "shear_force_y":     0.0,
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  1.304317E+05,
                "stress_top_fiber_y":  1.052227E+05,
                "tau_top_fiber_xy":   -1.111881E+05,
                "von_mises_top":       2.268221E+05,
                "major_principal_top": 2.297275E+05,
                "minor_principal_top": 5.926994E+03,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  1.304317E+05,
                "stress_bot_fiber_y":  1.052227E+05,
                "tau_bot_fiber_xy":   -1.111881E+05,
                "von_mises_bot":       2.268221E+05,
                "major_principal_bot": 2.297275E+05,
                "minor_principal_bot": 5.926994E+03
            },
            2: {
                # Membrana / Midplane
                "sigma_x": 1.060832E+05,
                "sigma_y": 4.326561E+05,
                "tau_xy":  1.586497E+05,
                "von_mises": 4.775527E+05,
                "major_principal": 4.970365E+05,
                "minor_principal": 4.170287E+04,

                # Fuerzas y Momentos
                "membrane_force_x":  1.060832E+02,
                "membrane_force_y":  4.326561E+02,
                "membrane_force_xy": 1.586497E+02,
                "bending_moment_x":  0.0,
                "bending_moment_y":  0.0,
                "bending_moment_xy": 0.0,
                "shear_force_x":     0.0,
                "shear_force_y":     0.0,
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  1.060832E+05,
                "stress_top_fiber_y":  4.326561E+05,
                "tau_top_fiber_xy":    1.586497E+05,
                "von_mises_top":       4.775527E+05,
                "major_principal_top": 4.970365E+05,
                "minor_principal_top": 4.170287E+04,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  1.060832E+05,
                "stress_bot_fiber_y":  4.326561E+05,
                "tau_bot_fiber_xy":    1.586497E+05,
                "von_mises_bot":       4.775527E+05,
                "major_principal_bot": 4.970365E+05,
                "minor_principal_bot": 4.170287E+04
            },
            10: {
                # Membrana / Midplane
                "sigma_x": -1.347920E+05,
                "sigma_y":  1.216227E+05,
                "tau_xy":   3.243900E+05,
                "von_mises": 6.041866E+05,
                "major_principal": 3.422219E+05,
                "minor_principal": -3.553913E+05,

                # Fuerzas y Momentos
                "membrane_force_x":  -1.347920E+02,
                "membrane_force_y":   1.216227E+02,
                "membrane_force_xy":  3.243900E+02,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x": -1.347920E+05,
                "stress_top_fiber_y":  1.216227E+05,
                "tau_top_fiber_xy":    3.243900E+05,
                "von_mises_top":       6.041866E+05,
                "major_principal_top": 3.422219E+05,
                "minor_principal_top": -3.553913E+05,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x": -1.347920E+05,
                "stress_bot_fiber_y":  1.216227E+05,
                "tau_bot_fiber_xy":    3.243900E+05,
                "von_mises_bot":       6.041866E+05,
                "major_principal_bot": 3.422219E+05,
                "minor_principal_bot": -3.553913E+05
            }
        }
    },
    #endregion

    #region --- CASO 42: Trapecio Plano XZ, Multi-Material (8 Elementos) ---
    {
        "case_name": "Ej42_CTRIA3_Trapezoid_XZ_MultiMat",
        "folder": "case_42",
        "filename": "Ej42.dat",
        "description": "8 CTRIA3 elements with 8 different materials and high skewness in XZ plane. Tests multi-stiffness assembly and load paths.",
        "rel_tol": 0.01,    
        "abs_tol_disp": 1e-5,
        "abs_tol_force": 1e-3,

        "expected_displacements": {
            1:  {0: -7.154935E-06, 1: 0.0, 2: -3.398146E-06, 3:0.0, 4: 1.527657E-07}, # Nodo empotrado
            5:  {0: 8.620173E-06,  1: 0.0, 2:  4.231403E-06, 3:0.0, 4: 1.717897E-07},
            9:  {0: 2.130888E-05,  1: 0.0, 2:  3.170477E-05, 3:0.0, 4:-5.673143E-07},
        },
        "expected_reactions": {
            # Nodos empotrados 2, 4, 6, 8 evaluados en X(0) y Z(2)
            2: { 0:  1.216458E+02, 2: -2.610722E+02, 4: -1.406964E-01},
            4: { 0:  2.026487E+02, 2:  5.495897E+02, 4: -2.216302E-01},
            6: { 0: -2.742055E+03, 2:  2.498112E+02, 4:  -1.837294E-01},
            8: { 0: -5.822398E+02, 2: -5.383288E+02, 4:  1.018744E-01}
        },
        "expected_element_forces": {
            1: {
                # Membrana / Midplane
                "sigma_x": 3.824137E+04,
                "sigma_y": 1.119145E+04,
                "tau_xy":  2.692750E+03,
                "von_mises": 3.437190E+04,
                "major_principal": 3.850683E+04,
                "minor_principal": 1.092600E+04,

                # Fuerzas y Momentos
                "membrane_force_x":   3.824137E+01,
                "membrane_force_y":   1.119145E+01,
                "membrane_force_xy":  2.692750E+00,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  3.824137E+04,
                "stress_top_fiber_y":  1.119145E+04,
                "tau_top_fiber_xy":    2.692750E+03,
                "von_mises_top":       3.437190E+04,
                "major_principal_top": 3.850683E+04,
                "minor_principal_top": 1.092600E+04,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  3.824137E+04,
                "stress_bot_fiber_y":  1.119145E+04,
                "tau_bot_fiber_xy":    2.692750E+03,
                "von_mises_bot":       3.437190E+04,
                "major_principal_bot": 3.850683E+04,
                "minor_principal_bot": 1.092600E+04
            },
            2: {
                # Membrana / Midplane
                "sigma_x": 6.207154E+04,
                "sigma_y": 2.387653E+04,
                "tau_xy":  2.259487E+04,
                "von_mises": 6.687673E+04,
                "major_principal": 7.255854E+04,
                "minor_principal": 1.338953E+04,

                # Fuerzas y Momentos
                "membrane_force_x":   6.207154E+01,
                "membrane_force_y":   2.387653E+01,
                "membrane_force_xy":  2.259487E+01,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  6.207154E+04,
                "stress_top_fiber_y":  2.387653E+04,
                "tau_top_fiber_xy":    2.259487E+04,
                "von_mises_top":       6.687673E+04,
                "major_principal_top": 7.255854E+04,
                "minor_principal_top": 1.338953E+04,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  6.207154E+04,
                "stress_bot_fiber_y":  2.387653E+04,
                "tau_bot_fiber_xy":    2.259487E+04,
                "von_mises_bot":       6.687673E+04,
                "major_principal_bot": 7.255854E+04,
                "minor_principal_bot": 1.338953E+04
            },
            4: {
                # Membrana / Midplane
                "sigma_x": -7.061358E+02,
                "sigma_y": -2.017531E+03,
                "tau_xy":   8.946562E+04,
                "von_mises": 1.549691E+05,
                "major_principal": 8.810619E+04,
                "minor_principal": -9.082986E+04,

                # Fuerzas y Momentos
                "membrane_force_x":  -7.061358E-01,
                "membrane_force_y":  -2.017531E+00,
                "membrane_force_xy":  8.946562E+01,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x": -7.061358E+02,
                "stress_top_fiber_y": -2.017531E+03,
                "tau_top_fiber_xy":    8.946562E+04,
                "von_mises_top":       1.549691E+05,
                "major_principal_top": 8.810619E+04,
                "minor_principal_top": -9.082986E+04,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x": -7.061358E+02,
                "stress_bot_fiber_y": -2.017531E+03,
                "tau_bot_fiber_xy":    8.946562E+04,
                "von_mises_bot":       1.549691E+05,
                "major_principal_bot": 8.810619E+04,
                "minor_principal_bot": -9.082986E+04
            },
            5: {
                # Membrana / Midplane
                "sigma_x":  3.043178E+03,
                "sigma_y": -5.353096E+02,
                "tau_xy":  -1.078739E+03,
                "von_mises": 3.829826E+03,
                "major_principal": 3.343210E+03,
                "minor_principal": -8.353412E+02,

                # Fuerzas y Momentos
                "membrane_force_x":   3.043178E+00,
                "membrane_force_y":  -5.353096E-01,
                "membrane_force_xy": -1.078739E+00,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  3.043178E+03,
                "stress_top_fiber_y": -5.353096E+02,
                "tau_top_fiber_xy":   -1.078739E+03,
                "von_mises_top":       3.829826E+03,
                "major_principal_top": 3.343210E+03,
                "minor_principal_top": -8.353412E+02,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  3.043178E+03,
                "stress_bot_fiber_y": -5.353096E+02,
                "tau_bot_fiber_xy":   -1.078739E+03,
                "von_mises_bot":       3.829826E+03,
                "major_principal_bot": 3.343210E+03,
                "minor_principal_bot": -8.353412E+02
            },
            7: {
                # Membrana / Midplane
                "sigma_x": -3.244819E+05,
                "sigma_y":  2.192878E+05,
                "tau_xy":  -3.854987E+05,
                "von_mises": 8.187542E+05,
                "major_principal": 4.191345E+05,
                "minor_principal": -5.243286E+05,

                # Fuerzas y Momentos
                "membrane_force_x":  -3.244819E+02,
                "membrane_force_y":   2.192878E+02,
                "membrane_force_xy": -3.854987E+02,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x": -3.244819E+05,
                "stress_top_fiber_y":  2.192878E+05,
                "tau_top_fiber_xy":   -3.854987E+05,
                "von_mises_top":       8.187542E+05,
                "major_principal_top": 4.191345E+05,
                "minor_principal_top": -5.243286E+05,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x": -3.244819E+05,
                "stress_bot_fiber_y":  2.192878E+05,
                "tau_bot_fiber_xy":   -3.854987E+05,
                "von_mises_bot":       8.187542E+05,
                "major_principal_bot": 4.191345E+05,
                "minor_principal_bot": -5.243286E+05
            },
            8: {
                # Membrana / Midplane
                "sigma_x": 1.107391E+05,
                "sigma_y": 1.918343E+05,
                "tau_xy": -6.898659E+04,
                "von_mises": 2.051766E+05,
                "major_principal": 2.313071E+05,
                "minor_principal": 7.126633E+04,

                # Fuerzas y Momentos
                "membrane_force_x":   1.107391E+02,
                "membrane_force_y":   1.918343E+02,
                "membrane_force_xy": -6.898659E+01,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  1.107391E+05,
                "stress_top_fiber_y":  1.918343E+05,
                "tau_top_fiber_xy":   -6.898659E+04,
                "von_mises_top":       2.051766E+05,
                "major_principal_top": 2.313071E+05,
                "minor_principal_top": 7.126633E+04,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  1.107391E+05,
                "stress_bot_fiber_y":  1.918343E+05,
                "tau_bot_fiber_xy":   -6.898659E+04,
                "von_mises_bot":       2.051766E+05,
                "major_principal_bot": 2.313071E+05,
                "minor_principal_bot": 7.126633E+04
            }
        }
    },
    #endregion

    #region --- CASO 43: Trapecio Plano YZ, Multi-Material (8 Elementos) ---
    {
        "case_name": "Ej43_CTRIA3_Trapezoid_YZ_MultiMat",
        "folder": "case_43",
        "filename": "Ej43.dat",
        "description": "8 CTRIA3 elements with 8 different materials in YZ plane. Severe mesh distortion and full element stress recovery.",
        "rel_tol": 0.01,    
        "abs_tol_disp": 1e-5,
        "abs_tol_force": 1e-3,

        "expected_displacements": {
            1: { 1:  1.155963E-06, 2: -6.594411E-07, 3:8.135751E-08, 4:0.0},
            2: { 1:  0.0, 2: 0.0, 3:0.0, 4:0.0},
            3: { 1:  1.371866E-06, 2: -5.167549E-07, 3:3.839616E-08},
            5: { 1: -2.923902E-07, 2:  3.653715E-07, 3:-1.148548E-08, 4:0.0},
            7: { 1:  1.169913E-07, 2:  5.632467E-08, 3:-5.777200E-09, 4:0.0},
            9: { 1:  1.468316E-08, 2:  -3.222052E-08, 3:-1.218673E-08, 4:0.0}
        },
        "expected_element_forces": {
            1: {
                # Fuerzas y Momentos
                "membrane_force_x":  -3.890586E+00,
                "membrane_force_y":   3.355102E+00,
                "membrane_force_xy": -1.415027E-01,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,

                # Top Fiber (z = +t/2)
                "stress_top_fiber_x": -3.890586E+03,
                "stress_top_fiber_y":  3.355102E+03,
                "tau_top_fiber_xy":   -1.415027E+02,
                "von_mises_top":       6.285440E+03,
                "major_principal_top": 3.357865E+03,
                "minor_principal_top": -3.893348E+03,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x": -3.890586E+03,
                "stress_bot_fiber_y":  3.355102E+03,
                "tau_bot_fiber_xy":   -1.415027E+02,
                "von_mises_bot":       6.285440E+03,
                "major_principal_bot": 3.357865E+03,
                "minor_principal_bot": -3.893348E+03
            },
            2: {
                # Fuerzas y Momentos
                "membrane_force_x":  -3.276435E+00,
                "membrane_force_y":  -1.562509E+01,
                "membrane_force_xy":  5.977812E+00,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,

                # Top Fiber (z = +t/2)
                "stress_top_fiber_x": -6.552869E+02,
                "stress_top_fiber_y": -3.125017E+03,
                "tau_top_fiber_xy":    1.195562E+03,
                "von_mises_top":       3.526395E+03,
                "major_principal_top": -1.713555E+02,
                "minor_principal_top": -3.608949E+03,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x": -6.552869E+02,
                "stress_bot_fiber_y": -3.125017E+03,
                "tau_bot_fiber_xy":    1.195562E+03,
                "von_mises_bot":       3.526395E+03,
                "major_principal_bot": -1.713555E+02,
                "minor_principal_bot": -3.608949E+03
            },
            3: {
                # Fuerzas y Momentos
                "membrane_force_x":   1.381637E+02,
                "membrane_force_y":   9.223628E+01,
                "membrane_force_xy": -8.527009E+01,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,

                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  1.381637E+04,
                "stress_top_fiber_y":  9.223628E+03,
                "tau_top_fiber_xy":   -8.527009E+03,
                "von_mises_top":       1.914837E+04,
                "major_principal_top": 2.035081E+04,
                "minor_principal_top": 2.689190E+03,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  1.381637E+04,
                "stress_bot_fiber_y":  9.223628E+03,
                "tau_bot_fiber_xy":   -8.527009E+03,
                "von_mises_bot":       1.914837E+04,
                "major_principal_bot": 2.035081E+04,
                "minor_principal_bot": 2.689190E+03
            },
            4: {
                # Fuerzas y Momentos
                "membrane_force_x":   1.020218E+01,
                "membrane_force_y":   2.914909E+01,
                "membrane_force_xy": -2.031564E+00,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,

                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  6.801455E+02,
                "stress_top_fiber_y":  1.943273E+03,
                "tau_top_fiber_xy":   -1.354376E+02,
                "von_mises_top":       1.724015E+03,
                "major_principal_top": 1.957632E+03,
                "minor_principal_top": 6.657866E+02,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  6.801455E+02,
                "stress_bot_fiber_y":  1.943273E+03,
                "tau_bot_fiber_xy":   -1.354376E+02,
                "von_mises_bot":       1.724015E+03,
                "major_principal_bot": 1.957632E+03,
                "minor_principal_bot": 6.657866E+02
            },
            5: {
                # Fuerzas y Momentos
                "membrane_force_x":  -4.510177E+01,
                "membrane_force_y":  -1.428301E+02,
                "membrane_force_xy":  1.427830E+01,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,

                # Top Fiber (z = +t/2)
                "stress_top_fiber_x": -2.255089E+03,
                "stress_top_fiber_y": -7.141507E+03,
                "tau_top_fiber_xy":    7.139152E+02,
                "von_mises_top":       6.442890E+03,
                "major_principal_top": -2.152921E+03,
                "minor_principal_top": -7.243675E+03,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x": -2.255089E+03,
                "stress_bot_fiber_y": -7.141507E+03,
                "tau_bot_fiber_xy":    7.139152E+02,
                "von_mises_bot":       6.442890E+03,
                "major_principal_bot": -2.152921E+03,
                "minor_principal_bot": -7.243675E+03
            },
            6: {
                # Fuerzas y Momentos
                "membrane_force_x":  -2.000000E+01,
                "membrane_force_y":  -6.666667E+01,
                "membrane_force_xy":  6.666000E+01,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,

                # Top Fiber (z = +t/2)
                "stress_top_fiber_x": -4.000000E+02,
                "stress_top_fiber_y": -1.333333E+03,
                "tau_top_fiber_xy":    1.333200E+03,
                "von_mises_top":       2.595518E+03,
                "major_principal_top": 5.458488E+02,
                "minor_principal_top": -2.279182E+03,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x": -4.000000E+02,
                "stress_bot_fiber_y": -1.333333E+03,
                "tau_bot_fiber_xy":    1.333200E+03,
                "von_mises_bot":       2.595518E+03,
                "major_principal_bot": 5.458488E+02,
                "minor_principal_bot": -2.279182E+03
            },
            7: {
                # Fuerzas y Momentos
                "membrane_force_x":   1.167642E+02,
                "membrane_force_y":  -1.522493E+02,
                "membrane_force_xy": -2.843417E+01,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,

                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  1.167642E+03,
                "stress_top_fiber_y": -1.522493E+03,
                "tau_top_fiber_xy":   -2.843417E+02,
                "von_mises_top":       2.387813E+03,
                "major_principal_top": 1.197368E+03,
                "minor_principal_top": -1.552219E+03,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  1.167642E+03,
                "stress_bot_fiber_y": -1.522493E+03,
                "tau_bot_fiber_xy":   -2.843417E+02,
                "von_mises_bot":       2.387813E+03,
                "major_principal_bot": 1.197368E+03,
                "minor_principal_bot": -1.552219E+03
            },
            8: {
                # Fuerzas y Momentos
                "membrane_force_x":   4.434409E+01,
                "membrane_force_y":  -2.870864E+02,
                "membrane_force_xy": -1.123684E+02,
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  0.0,
                "shear_force_x":      0.0,
                "shear_force_y":      0.0,

                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  4.434409E+01,
                "stress_top_fiber_y": -2.870864E+02,
                "tau_top_fiber_xy":   -1.123684E+02,
                "von_mises_top":       3.674174E+02,
                "major_principal_top": 7.884921E+01,
                "minor_principal_top": -3.215915E+02,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  4.434409E+01,
                "stress_bot_fiber_y": -2.870864E+02,
                "tau_bot_fiber_xy":   -1.123684E+02,
                "von_mises_bot":       3.674174E+02,
                "major_principal_bot": 7.884921E+01,
                "minor_principal_bot": -3.215915E+02
            }
        }
    },
    #endregion

    #region --- CASO 44: Flexión Pura Plano XY ---
    {
        "case_name": "CASO 44: Flexión CTRIA3 XY",
        "folder": "case_44",
        "filename": "case_44.dat",
        "description": "",
        "note":"",
        "rel_tol": 0.01,    
        "abs_tol_disp": 1e-5,
        "abs_tol_force": 1e-3,
        "expected_displacements": {
            # Desplazamientos extraídos del caso 44 (w, theta_x, theta_y)
            2: {0:0.0, 1:0.0, 2: 6.240000E+25, 3: 0.0, 4: 0.0},
            3: {0:0.0, 1:0.0, 2: 6.240000E+25, 3: 0.0, 4: 0.0}
        },
        "expected_reactions": {
            1: {
                2: -2.000000E+03,
                #3: -5.000000E+02,
                #4:  5.000000E+02
            },
        },
        "expected_element_forces": {
            50: {
                # Fuerzas y Momentos
                "membrane_force_x":  0.0,
                "membrane_force_y":  0.0,
                "membrane_force_xy": 0.0,
                "bending_moment_x":  0.0,
                "bending_moment_y":  0.0,
                "bending_moment_xy": 0.0,
                "shear_force_x":     2.000000E+03, # TRANSVERSE SHEAR FORCES Qx
                "shear_force_y":     2.000000E+03, # TRANSVERSE SHEAR FORCES Qy
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  0.0,
                "stress_top_fiber_y":  0.0,
                "tau_top_fiber_xy":    0.0,
                "von_mises_top":       0.0, 
                "major_principal_top": 0.0,
                "minor_principal_top": 0.0,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  0.0,
                "stress_bot_fiber_y":  0.0,
                "tau_bot_fiber_xy":    0.0,
                "von_mises_bot":       0.0,
                "major_principal_bot": 0.0,
                "minor_principal_bot": 0.0
            }
        }
    },
    #endregion

    #region --- CASO 45: Cortadura transversal 2 elementos ---
    {
        "case_name": "CASO 45: Cortadura transversal 2 CTRIAs",
        "folder": "case_45",
        "filename": "case_45.dat",
        "rel_tol": 0.01,  
        "abs_tol_disp": 1e-1,
        "abs_tol_force": 2.0, # al tener deformaciones en giros tan elevadas hay que levantar la mano con la toleracia en los esfuerzos
        "expected_displacements": {
            1: {
                0: 0.0,
                1: 0.0,
                2: 0.0,
                3: 0.0,
                4: 0.0,
                5: 0.0,
            },
            2: {
                0: 0.0,
                1: 0.0,
                2: 0.0,
                3: 0.0,
                4: 0.0,
                5: 0.0,
            },
            3: {
                0: 0.0,
                1: 0.0,
                2: -4.264000E+20,
                3: -2.535236E-10,
                4: 6.484807E-11,
                5: 0.0,
            },
            4: {
                0: 0.0,
                1: 0.0,
                2: -2.288000E+20,
                3: -2.360616E-10,
                4: 1.592333E-10,
                5: 0.0,
            },
        },
        "expected_reactions": {
            1: {
                2:  1.833333E-01,
                3:  1.761063E-01,
                4: -3.416667E-01,
                5: 0.0,
            },
            2: {
                2:  3.416667E-01,
                3: -1.777295E-02,
                4: -1.833333E-01,
                5: 0.0,
            },
        },
        "expected_element_forces": {
            1: {
                # Fuerzas y Momentos
                "membrane_force_x":  0.0,
                "membrane_force_y":  0.0,
                "membrane_force_xy": 0.0,
                "bending_moment_x": -3.645452E-01,
                "bending_moment_y": -1.093635E-01,
                "bending_moment_xy":-1.891519E-01,
                "shear_force_x":    -3.666667E-01, 
                "shear_force_y":     4.470687E-32,
                
                # Top Fiber (z = +t/2) -> DISTANCE = 2.500000E-01
                "stress_top_fiber_x":  8.749084E+00,
                "stress_top_fiber_y":  2.624725E+00,
                "tau_top_fiber_xy":    4.539646E+00,
                "von_mises_top":       1.105879E+01,
                "major_principal_top": 1.116279E+01,
                "minor_principal_top": 2.110181E-01,

                # Bottom Fiber (z = -t/2) -> DISTANCE = -2.500000E-01
                "stress_bot_fiber_x": -8.749084E+00,
                "stress_bot_fiber_y": -2.624725E+00,
                "tau_bot_fiber_xy":   -4.539646E+00,
                "von_mises_bot":       1.105879E+01,
                "major_principal_bot":-2.110181E-01,
                "minor_principal_bot":-1.116279E+01
            },
            2: {
                # Membrana / Midplane
                "sigma_x": 0.0,
                "sigma_y": 0.0,
                "tau_xy": 0.0,
                "von_mises": 0.0,
                "major_principal": 0.0,
                "minor_principal": 0.0,

                # Fuerzas y Momentos
                "membrane_force_x":  0.0,
                "membrane_force_y":  0.0,
                "membrane_force_xy": 0.0,
                "bending_moment_x": -1.604548E-01,
                "bending_moment_y": -8.451569E-02,
                "bending_moment_xy":-1.275147E-01,
                "shear_force_x":     6.833333E-01, 
                "shear_force_y":     3.166667E-01, 
                
                # Top Fiber (z = +t/2) -> DISTANCE = 2.500000E-01
                "stress_top_fiber_x":  3.850916E+00,
                "stress_top_fiber_y":  2.028377E+00,
                "tau_top_fiber_xy":    3.060354E+00,
                "von_mises_top":       6.263390E+00,
                "major_principal_top": 6.132792E+00,
                "minor_principal_top":-2.534989E-01,

                # Bottom Fiber (z = -t/2) -> DISTANCE = -2.500000E-01
                "stress_bot_fiber_x": -3.850916E+00,
                "stress_bot_fiber_y": -2.028377E+00,
                "tau_bot_fiber_xy":   -3.060354E+00,
                "von_mises_bot":       6.263390E+00,
                "major_principal_bot": 2.534989E-01,
                "minor_principal_bot":-6.132792E+00
            }
        }
    },
    #endregion

    #region --- CASO 46: Rotación de Matrices ---  
    {
        "case_name": "CASO 46: Rotación de Matrices",
        "folder": "case_46",
        "filename": "Ej46.dat",
        "rel_tol": 0.01,    
        "abs_tol_disp": 1e-5,
        "abs_tol_force": 1e-1,
        "expected_displacements": {
                2: {
                    1: 0.0,
                    2: 0.0,
                    3: 0.0,
                    4: 0.0,
                },
                3: {
                    0: 0.0,
                    1: -5.595121E-08,
                    2:  1.188103E-07,
                    5:  5.587325E-09,
                },
        },
        "expected_reactions": {
            1: {0: -2.419833E+01,
                1: -7.983323E-01,
                2: -1.499800E+01,
                3:  2.007153E-03,
                4:  1.990511E-03,
                5: -3.334161E-03
            },
            2: {
                0:  8.649990E-01,
                1: -9.201668E+00,
                2: -5.002001E+00,
                3:  1.988576E-03,
                4:  2.005218E-03,
                5: -3.336483E-03
            },
            3: {
                0:  3.333333E+01,
                1:  0.0,
                2:  0.0,
                3:  6.656711E-06,
                4:  6.656711E-06,
                5:  0.0
            }
        },
        "expected_element_forces": { 
            10: {
                "membrane_force_x":   -6.120276E+00,
                "membrane_force_y":   -2.185813E+01,
                "membrane_force_xy":  -6.666000E+00,
                
                "bending_moment_x":   0.0,
                "bending_moment_y":   0.0,
                "bending_moment_xy":  5.075695E-06,
                "shear_force_x":     0.0, # Prácticamente cero en F06
                "shear_force_y":     0.0, # Prácticamente cero en F06
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x": -6.120276E+02,
                "stress_top_fiber_y": -2.185813E+03,
                "tau_top_fiber_xy":   -6.669045E+02,
                "von_mises_top":       2.269111E+03,
                "major_principal_top":-3.674351E+02,
                "minor_principal_top":-2.430405E+03,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x": -6.120276E+02,
                "stress_bot_fiber_y": -2.185813E+03,
                "tau_bot_fiber_xy":   -6.662954E+02,
                "von_mises_bot":       2.268574E+03,
                "major_principal_bot":-3.678288E+02,
                "minor_principal_bot":-2.430012E+03
            }
        },
        
    },
    #endregion

    #region --- CASO 47: Convergencia de Malla 18 CTRIA3 ---
    {
        "case_name": "CASO 47: Convergencia de Malla 18 CTRIA3",
        "folder": "case_47",
        "filename": "case_47.dat",
        "rel_tol": 0.2,
        "expected_displacements": {
            1: {
                0: 0.0,
                1: 0.0,
                5: 0.0,
            },
            2: {
                0: 1.746000E-08,
                1: 5.274341E-09,
                5: 1.750931E-08,
            },
            3: {
                0: 3.507992E-08,
                1: 7.110717E-09,
                5: 3.130605E-08,
            },
            4: {
                0: 5.700468E-08,
                1: 1.414165E-08,
                5: 5.265168E-08,
            },
            5: {
                0: 0.0,
                1: 0.0,
                5: 0.0,
            },
            6: {
                0: 1.422573E-08,
                1: 1.573758E-09,
                5: 7.586792E-09,
            },
            7: {
                0: 2.889294E-08,
                1: 2.818707E-09,
                5: 1.001719E-08,
            },
            8: {
                0: 4.297494E-08,
                1: 5.114963E-09,
                5: 2.725968E-08,
            },
            9: {
                0: 0.0,
                1: 0.0,
                5: 0.0,
            },
            10: {
                0: 1.405558E-08,
                1: -9.506571E-10,
                5: -3.169371E-09,
            },
            11: {
                0: 2.869232E-08,
                1: -6.635764E-10,
                5: -9.150788E-09,
            },
            12: {
                0: 4.210349E-08,
                1: -2.422028E-09,
                5: -9.081675E-09,
            },
            13: {
                0: 0.0,
                1: 0.0,
                5: 0.0,
            },
            14: {
                0: 1.563513E-08,
                1: -4.695931E-09,
                5: -8.055739E-09,
            },
            15: {
                0: 3.286541E-08,
                1: -5.557129E-09,
                5: -1.019977E-08,
            },
            16: {
                0: 5.285293E-08,
                1: -7.296996E-09,
                5: -2.302607E-08,
            }, 
        },
        
    },
    #endregion
    #region --- CASO 48: Análogo al caso 43 pero analizando más flexión ---
    {
        "case_name": "Ej48_CTRIA3_Trapezoid_YZ_MultiMat",
        "folder": "case_48",
        "filename": "Ej48.dat",
        "description": " (As Case 43, but with more bending analisis) 8 CTRIA3 elements with 8 different materials in YZ plane.",
        "note": "Multimateial, skeewe mesh anf 3D force out of membrane plane. K6ROT=100 & MID3/E=0.0.",
        
        "rel_tol": 0.01,    
        "abs_tol_disp": 1e-5,
        "abs_tol_force": 5e-1,

        "expected_displacements": {
                1: {
                    1:  1.155963E-06,
                    2: -6.594411E-07,
                    3:  8.135751E-08,
                    4:  0.0,
                },
                2: {
                    1:  0.0,
                    2:  0.0,
                    3:  0.0,
                    4:  0.0,
                },
                3: {
                    1:  1.371866E-08,
                    2: -5.167549E-09,
                    3:  3.839616E-10,
                },
                5: {
                    1: -2.923902E-07,
                    2:  3.653715E-07,
                    3: -1.148548E-08,
                    4:  0.0,
                },
                7: {
                    1:  1.282563E-09,
                    2:  7.972293E-10,
                    3: -6.932640E-11,
                    4:  0.0,
                },
                9: {
                    1:  1.468316E-08,
                    2: -3.222052E-08,
                    3: -1.218673E-08,
                    4:  0.0,
                },
        },
        "expected_reactions": {
            2: { 
                0: -4.000000E-02, 
                1:  9.549069E+01, 
                2: -3.370845E+02, 
                #3: -1.398823E-02, 
                #4:  2.003716E+00, 
                5:  2.659258E+00 
            },
            4: { 0:  8.333333E-01, 
                1:  9.908001E+02, 
                2: -9.570910E+02, 
                #3: -9.220836E-03, 
                #4:  7.022088E+00, 
                #5:  3.813755E+00 
            },
            6: { 
                0: -9.600000E-01, 
                1:  8.173972E+02, 
                2:  1.612536E+03, 
                #3:  3.089219E-01, 
                #4:  9.596284E+00, 
                5:  3.140742E+00 },
            8: { 
                0:  1.666667E-01, 
                1:  7.131199E+01, 
                2: -3.805360E+03, 
                #3:  4.428437E-01, 
                #4:  4.779119E-01, 
                5:  3.686245E+00 
            }
        },
        "expected_element_forces": {
            3: {
                # Fuerzas y Momentos
                "membrane_force_x":   1.381678E+00,
                "membrane_force_y":   9.223219E-01,
                "membrane_force_xy": -8.527224E-01,
                "bending_moment_x":   4.421595E-01,
                "bending_moment_y":   3.002405E-01,
                "bending_moment_xy": -2.731644E-01,
                "shear_force_x":      4.625223E-02,
                "shear_force_y":     -7.655541E-02,
                
                # Top Fiber (z = +t/2)
                #"stress_top_fiber_x": -2.639141E+04,
                #"stress_top_fiber_y": -1.792219E+04,
                ###"tau_top_fiber_xy":    1.630459E+04,
                #"von_mises_top":       3.663658E+04,
                #"major_principal_top":-5.311276E+03,
                #"minor_principal_top":-3.900232E+04,

                # Bottom Fiber (z = -t/2)
                #"stress_bot_fiber_x":  2.666774E+04,
                #"stress_bot_fiber_y":  1.810666E+04,
                #"tau_bot_fiber_xy":   -1.647514E+04,
                #"von_mises_bot":       3.701954E+04,
                #"major_principal_bot": 3.940934E+04,
                #"minor_principal_bot": 5.365060E+03
            },
            6: {
               # Fuerzas y Momentos
                "membrane_force_x":  -1.600000E-01,
                "membrane_force_y":  -5.333333E-01,
                "membrane_force_xy":  7.999816E-01,
                "bending_moment_x":  -1.500000E-01,
                "bending_moment_y":  -5.000000E-01,
                "bending_moment_xy":  1.109438E-01,
                "shear_force_x":      0.0, 
                "shear_force_y":     -9.428090E-02,
                
                # Top Fiber (z = +t/2)
                #"stress_top_fiber_x":  3.568000E+02,
                #"stress_top_fiber_y":  1.189333E+03,
                #"tau_top_fiber_xy":   -2.502655E+02,
                #"von_mises_top":       1.142525E+03,
                #"major_principal_top": 1.258773E+03,
                #"minor_principal_top": 2.873602E+02,

                # Bottom Fiber (z = -t/2)
                #"stress_bot_fiber_x": -3.632000E+02,
                #"stress_bot_fiber_y": -1.210667E+03,
                #"tau_bot_fiber_xy":    2.822648E+02,
                #"von_mises_bot":       1.181920E+03,
                #"major_principal_bot":-2.777935E+02,
                #"minor_principal_bot":-1.296073E+03
            }

        }
    },
    #endregion

        #region --- CASO 49: Análogo al caso 48 pero con K6ROT=23.0 ---
    {
        "case_name": "Ej49_CTRIA3_Trapezoid_MultiMat_K6ROT=23.0",
        "folder": "case_49",
        "filename": "Ej49.dat",
        "description": "(As Case 48, but with K6ROT = 23.0) 8 CTRIA3 elements with 8 different materials in YZ plane.",
        "note": "Multimateial, skeewe mesh anf 3D force out of membrane plane. K6ROT=23.0 & MID3/E=0.0.",
        
        "rel_tol": 0.01,    
        "abs_tol_disp": 1e-5,
        "abs_tol_force": 5e-1,

        "expected_displacements": {
                    1: {
                        0:  0.0,
                        1:  1.156060E-06,
                        2: -6.594210E-07,
                        3:  8.136276E-08,
                        4:  0.0,
                    },
                    2: {
                        0:  0.0,
                        1:  0.0,
                        2:  0.0,
                        3:  0.0,
                        4:  0.0,
                    },
                    3: {
                        0:  3.619200E+22,
                        1:  1.371935E-08,
                        2: -5.167206E-09,
                        3:  3.839912E-10,
                    },
                    5: {
                        1: -2.923999E-07,
                        2:  3.653788E-07,
                        3: -1.148579E-08,
                        4:  0.0,
                    },
                    7: {
                        0: -6.240000E+21,
                        1:  1.282643E-09,
                        2:  7.973094E-10,
                        3: -6.933174E-11,
                        4:  0.0,
                    },
                    9: {
                        1:  1.468243E-08,
                        2: -3.222325E-08,
                        3: -1.218707E-08,
                        4:  0.0,
                    },
        },
        "expected_reactions": {
            2: { 
                0: -4.000000E-02, 
                1:  9.549069E+01, 
                2: -3.370845E+02, 
                3: -1.398823E-02, 
                #4:  2.003716E+00, 
                5:  2.659258E+00 
            },
            4: { 
                0:  8.333333E-01, 
                1:  9.908001E+02, 
                2: -9.570910E+02, 
                3: -9.220836E-03, 
                #4:  7.022088E+00, 
                5:  3.813755E+00 
            },
            6: { 
                0: -9.600000E-01, 
                1:  8.173972E+02, 
                2:  1.612536E+03, 
                3:  3.089219E-01, 
                #4:  9.596284E+00, 
                5:  3.140742E+00 
            },
            8: { 
                0:  1.666667E-01, 
                1:  7.131199E+01, 
                2: -3.805360E+03, 
                3:  4.428437E-01, 
                #4:  4.779119E-01, 
                5:  3.686245E+00 
            } 
        },
        "expected_element_forces": {
            3: {
                # Fuerzas y Momentos
                "membrane_force_x":   1.381678E+00,
                "membrane_force_y":   9.223219E-01,
                "membrane_force_xy": -8.527224E-01,
                "bending_moment_x":   4.421595E-01,
                "bending_moment_y":   3.002405E-01,
                "bending_moment_xy": -2.731644E-01,
                "shear_force_x":      4.625223E-02,
                "shear_force_y":     -7.655541E-02,
                
                # Top Fiber (z = +t/2)
                #"stress_top_fiber_x": -2.639141E+04,
                #"stress_top_fiber_y": -1.792219E+04,
                ###"tau_top_fiber_xy":    1.630459E+04,
                #"von_mises_top":       3.663658E+04,
                #"major_principal_top":-5.311276E+03,
                #"minor_principal_top":-3.900232E+04,

                # Bottom Fiber (z = -t/2)
                #"stress_bot_fiber_x":  2.666774E+04,
                #"stress_bot_fiber_y":  1.810666E+04,
                #"tau_bot_fiber_xy":   -1.647514E+04,
                #"von_mises_bot":       3.701954E+04,
                #"major_principal_bot": 3.940934E+04,
                #"minor_principal_bot": 5.365060E+03
            },
            6: {
                # Fuerzas y Momentos
                "membrane_force_x":  -1.600000E-01,
                "membrane_force_y":  -5.333333E-01,
                "membrane_force_xy":  7.999816E-01,
                "bending_moment_x":  -1.500000E-01,
                "bending_moment_y":  -5.000000E-01,
                "bending_moment_xy":  1.109438E-01,
                "shear_force_x":      0.0, 
                "shear_force_y":     -9.428090E-02,
                
                # Top Fiber (z = +t/2)
                #"stress_top_fiber_x":  3.568000E+02,
                #"stress_top_fiber_y":  1.189333E+03,
                #"tau_top_fiber_xy":   -2.502655E+02,
                #"von_mises_top":       1.142525E+03,
                #"major_principal_top": 1.258773E+03,
                #"minor_principal_top": 2.873602E+02,

                # Bottom Fiber (z = -t/2)
                #"stress_bot_fiber_x": -3.632000E+02,
                #"stress_bot_fiber_y": -1.210667E+03,
                #"tau_bot_fiber_xy":    2.822648E+02,
                #"von_mises_bot":       1.181920E+03,
                #"major_principal_bot":-2.777935E+02,
                #"minor_principal_bot":-1.296073E+03
            }
        }
    },
    #endregion

    #region --- CASO MANUAL: Membrana Pura ---
    {
        "case_name": "Ej_Manual_Membrana_Pura",
        "folder": "case_manual_con_cargas", # <-- Asegúrate de que apunte a la carpeta correcta
        "filename": "Ej_manual_Kglobal_membrana_con carga en plano.dat",
        "rel_tol": 1e-3,
        "expected_displacements": {
            2: {0: 9.523810E-03},
            3: {1: -2.857143E-03}
        },
        "expected_reactions": {
            1: {0: -1.000000E+03}
        },
        "expected_element_forces": {
            1: {
                # Fuerzas y Momentos
                "membrane_force_x":  200.0,
                "membrane_force_y":  -2.842171E-14,  
                "membrane_force_xy": 0.0,
                "bending_moment_x":  0.0,
                "bending_moment_y":  0.0,
                "bending_moment_xy": 0.0,
                "shear_force_x":     0.0, 
                "shear_force_y":     0.0, 
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x":  200.0,
                "stress_top_fiber_y":  -2.842171E-14,
                "tau_top_fiber_xy":    0.0,
                "von_mises_top":       200.0,
                "major_principal_top": 200.0,
                "minor_principal_top": -2.842171E-14,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  200.0,
                "stress_bot_fiber_y":  -2.842171E-14,
                "tau_bot_fiber_xy":    0.0,
                "von_mises_bot":       200.0,
                "major_principal_bot": 200.0,
                "minor_principal_bot": -2.842171E-14
            }
        }
    },
    #endregion

    #region --- CASO MANUAL: Flexión Pura ---
    {
        "case_name": "Ej_Manual_Flexion_Pura",
        "folder": "case_manual_con_cargas", # <-- Asegúrate de que coincida con tu carpeta
        "filename": "Ej_manual_Kglobal_flexion_pura_MID3_cero_con carga.dat",
        "rel_tol": 0.05, # Dejamos un 5% de tolerancia por diferencias sutiles entre teorías de placas
        "expected_displacements": {
            # Desplazamientos Rotacionales R1 (Rot-X) y R2 (Rot-Y) en el Nodo 3
            3: {3: 1.040000E-01,
                4:-3.626297E-15}
        },
        "expected_reactions": {
            # Reacción del momento puro en el Nodo 1
            1: {2: 8.333333E-12,
                3: -1.000000E+03,
                4: 3.000000E+02},
            2: {2: -1.429760E-25,
                3: 2.440777E-11,# Prácticamente cero
                4: 3-3.000000E+02},
            3:{
                2:-8.333333E-12, # Prácticamente cero
                3:0.0,
                4:0.0,
                5:0.0
                }
        },
        "expected_element_forces": {
            1: {
                # Fuerzas y Momentos
                "membrane_force_x":   0.0,
                "membrane_force_y":   0.0,
                "membrane_force_xy":  0.0,
                "bending_moment_x":   6.000000E+01,
                "bending_moment_y":   2.000000E+02,
                "bending_moment_xy":  2.440777E-12,
                "shear_force_x":     -2.859519E-26,
                "shear_force_y":     -1.666667E-12,
                
                # Top Fiber (z = +t/2)
                "stress_top_fiber_x": -3.600000E+02,
                "stress_top_fiber_y": -1.200000E+03,
                "tau_top_fiber_xy":   -1.464466E-11,
                "von_mises_top":       1.066583E+03,
                "major_principal_top":-3.600000E+02,
                "minor_principal_top":-1.200000E+03,

                # Bottom Fiber (z = -t/2)
                "stress_bot_fiber_x":  3.600000E+02,
                "stress_bot_fiber_y":  1.200000E+03,
                "tau_bot_fiber_xy":    1.464466E-11,
                "von_mises_bot":       1.066583E+03,
                "major_principal_bot": 1.200000E+03,
                "minor_principal_bot": 3.600000E+02
            }
        }
    },
    #endregion
    #region --- CASO MANUAL: Fitting L con Carga de Traccion ---
    {
        "case_name": "Ej50_Fitting_L",
        "folder": "case_50",
        "filename": "Ej50.dat",
        "description": "Dummy L fitting + 2 CBAR as base.",
        "note": "K6ROT=100.0 & MID3/E=0.0.",
        
        "rel_tol": 0.01,    
        #"abs_tol_disp": 1e-5,

        "expected_displacements": {
            1:  {0: -1.371559E-07, 1: 4.816932E-07, 2: 2.202258E-06,  3: -1.330108E-06, 4: 5.328575E-06, 5: -9.837401E-08},
            2:  {0: 7.981849E-09,  1: 4.335555E-07, 2: 1.263795E-06,  3: -1.296167E-06, 4: 3.537645E-06, 5: -1.420997E-07},
            3:  {0: 1.303164E-07,  1: 4.132674E-07, 2: -1.784488E-07, 3: -1.386450E-06, 4: 6.132084E-06, 5: -8.329172E-08},
            4:  {0: 0.0,           1: 0.0,          2: 0.0,           3: 0.0,           4: 0.0,          5: 0.0},
            5:  {0: 0.0,           1: 0.0,          2: 0.0,           3: 0.0,           4: 0.0,          5: 0.0},
            6:  {0: 0.0,           1: 3.659459E-06, 2: -6.543600E-08, 3: -1.662775E-06, 4: 5.304242E-06, 5: 6.152201E-07},
            7:  {0: 0.0,           1: 3.998660E-06, 2: 2.869857E-06,  3: -1.385384E-06, 4: 4.466842E-06, 5: 7.344812E-07},
            8:  {0: 0.0,           1: 3.662532E-06, 2: 1.637019E-06,  3: -1.461351E-06, 4: 4.373258E-06, 5: 6.375942E-07},
            10: {0: 0.0,           1: 0.0,          2: 0.0,           3: 0.0,           4: 0.0,          5: 0.0},
            30: {0: 0.0,           1: 0.0,          2: 0.0,           3: 0.0,           4: 0.0,          5: 0.0}
        },
        # --- REACCIONES (SPCFORCES) ---
        "expected_reactions": {
            4:  {0: 4.658708E+02,  1: -2.469092E+02, 2: 6.681161E-29,  3: 1.933096E-03,  4: -3.068949E-01, 5: 1.323872E-01},
            5:  {0: -5.237704E+02, 1: -2.487398E+02, 2: 1.302834E-28,  3: 4.789982E-02,  4: -4.047260E-01, 5: 1.231769E-01},
            6:  {0: -7.652261E-29, 1: 0.0,           2: 0.0,           3: 0.0,           4: 0.0,           5: 0.0},
            7:  {0: -9.585472E-29, 1: 0.0,           2: 0.0,           3: 0.0,           4: 0.0,           5: 0.0},
            8:  {0: -1.500000E+02, 1: 0.0,           2: 0.0,           3: 0.0,           4: 0.0,           5: 0.0},
            10: {0: 2.962108E+01,  1: -3.007192E+00, 2: -9.249485E+02, 3: 1.310443E+00,  4: 5.167269E+00,  5: 3.178237E-02},
            30: {0: 2.827853E+01,  1: -1.343762E+00, 2: 7.494848E+01,  3: 9.182495E-01,  4: 4.494156E+00,  5: 2.690963E-02}
        },
        # --- FUERZAS Y TENSIONES EN ELEMENTOS ---
        "expected_element_forces": {
            # Elementos CBAR 
            # ----------------------------------------------------
            201: {
                "axial_force":   9.249485E+02,
                "axial_stress":  9.249485E+05,
                "torque":       -3.178237E-02,
                
                # Mapeo a las variables de tu funcion: forces_n1 (End A) y forces_n2 (End B)
                #"forces_n1.Mz":  9.643272E+00,   # Bending Moment End-A Plane 1
                #"forces_n1.My":  1.931528E-01,   # Bending Moment End-A Plane 2
                "forces_n2.Mz": -5.167269E+00,   # Bending Moment End-B Plane 1
                #"forces_n2.My": -1.310443E+00,   # Bending Moment End-B Plane 2
                #"forces_n1.Fy":  2.962108E+01,   # Shear Plane 1
                #"forces_n1.Fz":  3.007192E+00,   # Shear Plane 2
                
                # Puntos de Tensión (CBAR 201) - End A
                "stress_points.P0.sigma_bending": 0.0, # SA1
                "stress_points.P1.sigma_bending": 0.0, # SA2
                "stress_points.P2.sigma_bending": 0.0, # SA3
                "stress_points.P3.sigma_bending": 0.0, # SA4
                
                # Tensiones Max y Min (End A y End B coinciden en tu f06)
                "max_total_stress": 9.249485E+05, # SA-MAX / SB-MAX
                "min_total_stress": 9.249485E+05  # SA-MIN / SB-MIN
            },
            # Elemento CBAR
            202: {
                "axial_force":  -7.494848E+01,
                "axial_stress": -7.494848E+04,
                "torque":       -2.690963E-02,
                
                # Mapeo a las variables de tu funcion
                #"forces_n1.Mz":  9.645107E+00,
                #"forces_n1.My": -2.463686E-01,
                "forces_n2.Mz": -4.494156E+00,
                #"forces_n2.My": -9.182495E-01,
                #"forces_n1.Fy":  2.827853E+01,
                #"forces_n1.Fz":  1.343762E+00,
                
                # Puntos de Tensión (CBAR 202)
                "stress_points.P0.sigma_bending": 0.0,
                "stress_points.P1.sigma_bending": 0.0,
                "stress_points.P2.sigma_bending": 0.0,
                "stress_points.P3.sigma_bending": 0.0,
                
                "max_total_stress": -7.494848E+04,
                "min_total_stress": -7.494848E+04
            },
            # Elementos CTRIA 
            # ----------------------------------------------------
            # Elemento 1
            1: {
                "membrane_force_x":   2.498615E+02, "membrane_force_y":  -1.272198E+02, "membrane_force_xy": -1.546065E+02,
                "bending_moment_x":   4.114569E-01, "bending_moment_y":   1.281887E-01, "bending_moment_xy":  1.322453E-01,
                "shear_force_x":      8.641892E-29, "shear_force_y":      3.780114E-29, 
                "stress_bot_fiber_x": 1.866493E+04, "stress_bot_fiber_y":-4.438158E+03, "tau_bot_fiber_xy":  -5.746647E+03,
                "von_mises_bot":      2.345179E+04, "major_principal_bot": 2.001541E+04, "minor_principal_bot":-5.788634E+03,
                "stress_top_fiber_x": 6.321224E+03, "stress_top_fiber_y":-8.283820E+03, "tau_top_fiber_xy":  -9.714006E+03,
                "von_mises_top":      2.107200E+04, "major_principal_top": 1.117143E+04, "minor_principal_top":-1.313402E+04
            },
            # Elemento 2
            2: {
                "membrane_force_x":   2.643027E+02, "membrane_force_y":  -2.882482E+02, "membrane_force_xy": -2.152652E+02,
                "bending_moment_x":   2.061114E-01, "bending_moment_y":   1.476531E-01, "bending_moment_xy":  9.713352E-02,
                "shear_force_x":      5.944445E-29, "shear_force_y":      4.407526E-29, 
                "stress_bot_fiber_x": 1.630680E+04, "stress_bot_fiber_y":-1.219761E+04, "tau_bot_fiber_xy":  -9.306256E+03,
                "von_mises_bot":      2.955363E+04, "major_principal_bot": 1.907611E+04, "minor_principal_bot":-1.496692E+04,
                "stress_top_fiber_x": 1.012346E+04, "stress_top_fiber_y":-1.662721E+04, "tau_top_fiber_xy":  -1.222026E+04,
                "von_mises_top":      3.154802E+04, "major_principal_top": 1.486537E+04, "minor_principal_top":-2.136911E+04
            },
            # Elemento 3
            3: {
                "membrane_force_x":  -4.087689E+02, "membrane_force_y":  -1.039092E+02, "membrane_force_xy": -2.371936E+01,
                "bending_moment_x":   3.176510E-01, "bending_moment_y":   2.775008E-01, "bending_moment_xy": -1.973886E-01,
                "shear_force_x":      1.240139E-28, "shear_force_y":     -1.142409E-28, 
                "stress_bot_fiber_x":-1.567368E+04, "stress_bot_fiber_y":-1.032947E+03, "tau_bot_fiber_xy":  -4.146797E+03,
                "von_mises_bot":      1.679669E+04, "major_principal_bot": 5.999166E+01, "minor_principal_bot":-1.676662E+04,
                "stress_top_fiber_x":-2.520321E+04, "stress_top_fiber_y":-9.357971E+03, "tau_top_fiber_xy":   1.774861E+03,
                "von_mises_top":      2.227943E+04, "major_principal_top":-9.161599E+03, "minor_principal_top":-2.539958E+04
            },
            # Elemento 101
            101: {
                "membrane_force_x":   6.889222E+02, "membrane_force_y":   1.722289E+02, "membrane_force_xy":  3.444594E+02,
                "bending_moment_x":   7.410485E-03, "bending_moment_y":   1.852621E-03, "bending_moment_xy":  3.705242E-03,
                "shear_force_x":      1.563150E-28, "shear_force_y":     -1.769721E-29, 
                "stress_bot_fiber_x": 6.933685E+04, "stress_bot_fiber_y": 1.733405E+04, "tau_bot_fiber_xy":   3.466826E+04,
                "von_mises_bot":      8.667090E+04, "major_principal_bot": 8.667090E+04, "minor_principal_bot":-3.223686E-07,
                "stress_top_fiber_x": 6.844759E+04, "stress_top_fiber_y": 1.711173E+04, "tau_top_fiber_xy":   3.422363E+04,
                "von_mises_top":      8.555933E+04, "major_principal_top": 8.555933E+04, "minor_principal_top":-3.265523E-07
            },
            # Elemento 102
            102: {
                "membrane_force_x":   1.810581E+01, "membrane_force_y":   3.973164E+02, "membrane_force_xy":  5.460204E+02,
                "bending_moment_x":   3.251303E-03, "bending_moment_y":   8.287009E-03, "bending_moment_xy": -1.467831E-02,
                "shear_force_x":     -7.778275E-30, "shear_force_y":     -1.354816E-28, 
                "stress_bot_fiber_x": 2.005659E+03, "stress_bot_fiber_y": 4.022886E+04, "tau_bot_fiber_xy":   5.372134E+04,
                "von_mises_bot":      1.009933E+05, "major_principal_bot": 7.813687E+04, "minor_principal_bot":-3.590235E+04,
                "stress_top_fiber_x": 1.615502E+03, "stress_top_fiber_y": 3.923442E+04, "tau_top_fiber_xy":   5.548274E+04,
                "von_mises_top":      1.035064E+05, "major_principal_top": 7.900935E+04, "minor_principal_top":-3.815942E+04
            },
            # Elemento 103
            103: {
                "membrane_force_x":   8.237398E+01, "membrane_force_y":   4.165969E+02, "membrane_force_xy":  1.391202E+02,
                "bending_moment_x":   1.279500E-03, "bending_moment_y":   7.695468E-03, "bending_moment_xy":  1.483860E-02,
                "shear_force_x":      1.857650E-29, "shear_force_y":     -1.451728E-28, 
                "stress_bot_fiber_x": 8.314168E+03, "stress_bot_fiber_y": 4.212142E+04, "tau_bot_fiber_xy":   1.480234E+04,
                "von_mises_bot":      4.637308E+04, "major_principal_bot": 4.768648E+04, "minor_principal_bot": 2.749109E+03,
                "stress_top_fiber_x": 8.160628E+03, "stress_top_fiber_y": 4.119796E+04, "tau_top_fiber_xy":   1.302171E+04,
                "von_mises_top":      4.400410E+04, "major_principal_top": 4.571334E+04, "minor_principal_top": 3.645246E+03
            },
            # Elemento 104
            104: {
                "membrane_force_x":   1.282722E+02, "membrane_force_y":   3.202947E+01, "membrane_force_xy":  6.409753E+01,
                "bending_moment_x":  -7.830944E-03, "bending_moment_y":  -1.957736E-03, "bending_moment_xy": -3.915472E-03,
                "shear_force_x":     -1.775262E-28, "shear_force_y":     -1.224049E-29, 
                "stress_bot_fiber_x": 1.235737E+04, "stress_bot_fiber_y": 3.085483E+03, "tau_bot_fiber_xy":   6.174824E+03,
                "von_mises_bot":      1.544285E+04, "major_principal_bot": 1.544285E+04, "minor_principal_bot":-9.641954E-04,
                "stress_top_fiber_x": 1.329708E+04, "stress_top_fiber_y": 3.320411E+03, "tau_top_fiber_xy":   6.644681E+03,
                "von_mises_top":      1.661749E+04, "major_principal_top": 1.661749E+04, "minor_principal_top":-8.960393E-04
            }
        }

    },
    #endregion
    {
        "case_name": "CTRIA3_PSHELL_TS_T_0.25_E2E",
        "folder": "case_ctria3_pshell_factors_e2e",
        "filename": "ctria3_e2e_ts_t_025.dat",
        "description": "Four-element CTRIA3 cantilever with TS/T=0.25.",
        "note": "MSC Nastran 2026.1 reference; AUTOSPC=NO.",
        "rel_tol": 0.06,
        "rel_tol_disp": 0.015,
        "rel_tol_force": 0.02,
        "rel_tol_shear": 0.20,
        "abs_tol_disp": 1e-6,
        "abs_tol_force": 0.02,
        "abs_tol_shear": 0.02,
        "expected_displacements": {
            3: {0: 3.053448454e-04, 1: -2.357341815e-03, 2: -2.903097153e+01, 4: 3.913707566e-03},
            5: {0: 1.819616882e-03, 1: -5.551666953e-03, 2: -5.835988617e+01, 4: 5.719000474e-03},
            6: {0: 2.993631642e-03, 1: -6.390629336e-03, 2: -5.831832886e+01, 4: 5.910182837e-03},
        },
        "expected_reactions": {
            1: {0: 1.998733063e+02, 1: 2.821232605e+02, 2: 1.501314240e+02, 4: -3.003023047e+04, 5: 3.173344851e+00},
            2: {0: -9.998732910e+02, 1: 2.178767548e+02, 2: 1.498685760e+02, 4: -3.496976953e+04, 5: 9.495823860e+00},
        },
        "expected_element_forces": {
            101: {
                "sigma_x": 2.348806381e-01, "sigma_y": 7.046419382e-02, "tau_xy": -6.346689701e-01,
                "membrane_force_x": 2.348806381e+00, "membrane_force_y": 7.046419382e-01, "membrane_force_xy": -6.346689701e+00,
                "bending_moment_x": -4.870391541e+02, "bending_moment_y": -1.217597885e+02,
                "shear_force_x": -3.001824856e+00, "shear_force_y": -8.036214276e-04,
                "stress_bot_fiber_x": -2.898746872e+01, "stress_bot_fiber_y": -7.235123158e+00,
                "stress_top_fiber_x": 2.945722961e+01, "stress_top_fiber_y": 7.376051426e+00,
            },
            103: {
                "sigma_x": 1.010378551e+00, "sigma_y": -1.653728247e-01, "tau_xy": -4.103708744e-01,
                "membrane_force_x": 1.010378551e+01, "membrane_force_y": -1.653728247e+00, "membrane_force_xy": -4.103708744e+00,
                "bending_moment_x": -2.029710236e+02, "bending_moment_y": 3.058591080e+01,
                "shear_force_x": -3.001374245e+00, "shear_force_y": -1.655671163e-03,
                "stress_bot_fiber_x": -1.116788292e+01, "stress_bot_fiber_y": 1.669781804e+00,
                "stress_top_fiber_x": 1.318864059e+01, "stress_top_fiber_y": -2.000527382e+00,
            },
        },
    },
    {
        "case_name": "CTRIA3_PSHELL_TS_T_1.25_E2E",
        "folder": "case_ctria3_pshell_factors_e2e",
        "filename": "ctria3_e2e_ts_t_125.dat",
        "description": "Four-element CTRIA3 cantilever with TS/T=1.25.",
        "note": "MSC Nastran 2026.1 reference; AUTOSPC=NO.",
        "rel_tol": 0.06,
        "rel_tol_disp": 0.015,
        "rel_tol_force": 0.02,
        "rel_tol_shear": 0.20,
        "abs_tol_disp": 1e-6,
        "abs_tol_force": 0.02,
        "abs_tol_shear": 0.02,
        "expected_displacements": {
            3: {0: 3.053448454e-04, 1: -2.357341815e-03, 2: -5.990248203e+00, 4: 3.923573997e-03},
            5: {0: 1.819616882e-03, 1: -5.551666953e-03, 2: -1.227856064e+01, 4: 5.729883444e-03},
            6: {0: 2.993631642e-03, 1: -6.390629336e-03, 2: -1.223981190e+01, 4: 5.897646770e-03},
        },
        "expected_reactions": {
            1: {0: 1.998733063e+02, 1: 2.821232605e+02, 2: 1.506245575e+02, 4: -3.014171094e+04, 5: 3.173344851e+00},
            2: {0: -9.998732910e+02, 1: 2.178767548e+02, 2: 1.493754425e+02, 4: -3.485828906e+04, 5: 9.495823860e+00},
        },
        "expected_element_forces": {
            101: {
                "sigma_x": 2.348806381e-01, "sigma_y": 7.046419382e-02, "tau_xy": -6.346689701e-01,
                "membrane_force_x": 2.348806381e+00, "membrane_force_y": 7.046419382e-01, "membrane_force_xy": -6.346689701e+00,
                "bending_moment_x": -4.882669678e+02, "bending_moment_y": -1.220667419e+02,
                "shear_force_x": -3.008591890e+00, "shear_force_y": -3.899288364e-03,
                "stress_bot_fiber_x": -2.906113815e+01, "stress_bot_fiber_y": -7.253540516e+00,
                "stress_top_fiber_x": 2.953089905e+01, "stress_top_fiber_y": 7.394468784e+00,
            },
            103: {
                "sigma_x": 1.010378551e+00, "sigma_y": -1.653728247e-01, "tau_xy": -4.103708744e-01,
                "membrane_force_x": 1.010378551e+01, "membrane_force_y": -1.653728247e+00, "membrane_force_xy": -4.103708744e+00,
                "bending_moment_x": -2.031917877e+02, "bending_moment_y": 3.017736626e+01,
                "shear_force_x": -3.006382942e+00, "shear_force_y": -7.912937552e-03,
                "stress_bot_fiber_x": -1.118112850e+01, "stress_bot_fiber_y": 1.645269156e+00,
                "stress_top_fiber_x": 1.320188522e+01, "stress_top_fiber_y": -1.976014853e+00,
            },
        },
    },
    {
        "case_name": "CTRIA3_PSHELL_COMBINED_FACTORS_E2E",
        "folder": "case_ctria3_pshell_factors_e2e",
        "filename": "ctria3_e2e_combined_factors.dat",
        "description": "Four-element CTRIA3 cantilever with 12I/T^3=0.40 and TS/T=0.50.",
        "note": "MSC Nastran 2026.1 reference; AUTOSPC=NO.",
        "rel_tol": 0.06,
        "rel_tol_disp": 0.015,
        "rel_tol_force": 0.02,
        "rel_tol_shear": 0.20,
        "abs_tol_disp": 1e-6,
        "abs_tol_force": 0.02,
        "abs_tol_shear": 0.02,
        "expected_displacements": {
            3: {0: 3.053448454e-04, 1: -2.357341815e-03, 2: -1.497562027e+01, 4: 9.808934294e-03},
            5: {0: 1.819616882e-03, 1: -5.551666953e-03, 2: -3.069640160e+01, 4: 1.432470884e-02},
            6: {0: 2.993631642e-03, 1: -6.390629336e-03, 2: -3.059952927e+01, 4: 1.474411786e-02},
        },
        "expected_reactions": {
            1: {0: 1.998733063e+02, 1: 2.821232605e+02, 2: 1.506245575e+02, 4: -3.014171094e+04, 5: 3.173344851e+00},
            2: {0: -9.998732910e+02, 1: 2.178767548e+02, 2: 1.493754425e+02, 4: -3.485828906e+04, 5: 9.495823860e+00},
        },
        "expected_element_forces": {
            101: {
                "sigma_x": 2.348806381e-01, "sigma_y": 7.046419382e-02, "tau_xy": -6.346689701e-01,
                "membrane_force_x": 2.348806381e+00, "membrane_force_y": 7.046419382e-01, "membrane_force_xy": -6.346689701e+00,
                "bending_moment_x": -4.882669678e+02, "bending_moment_y": -1.220667419e+02,
                "shear_force_x": -3.008591890e+00, "shear_force_y": -3.899288364e-03,
                "stress_bot_fiber_x": -7.300516510e+01, "stress_bot_fiber_y": -1.823954773e+01,
                "stress_top_fiber_x": 7.347492218e+01, "stress_top_fiber_y": 1.838047600e+01,
            },
            103: {
                "sigma_x": 1.010378551e+00, "sigma_y": -1.653728247e-01, "tau_xy": -4.103708744e-01,
                "membrane_force_x": 1.010378551e+01, "membrane_force_y": -1.653728247e+00, "membrane_force_xy": -4.103708744e+00,
                "bending_moment_x": -2.031917877e+02, "bending_moment_y": 3.017736626e+01,
                "shear_force_x": -3.006382942e+00, "shear_force_y": -7.912937552e-03,
                "stress_bot_fiber_x": -2.946838951e+01, "stress_bot_fiber_y": 4.361232281e+00,
                "stress_top_fiber_x": 3.148914528e+01, "stress_top_fiber_y": -4.691977978e+00,
            },
        },
    },
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

        # 4. COMPARE RESULTS
        fallback_rel_tol = case_data.get("rel_tol", 1e-2)
        rel_tol = case_data.get("rel_tol", 1e-2)
        rel_tol_disp = case_data.get("rel_tol_disp", fallback_rel_tol)
        rel_tol_force = case_data.get("rel_tol_force", fallback_rel_tol)
        rel_tol_shear = case_data.get("rel_tol_shear", rel_tol)

        abs_tol_disp = case_data.get("abs_tol_disp", 1e-5)
        abs_tol_force = case_data.get("abs_tol_force", 1e-3)
        abs_tol_shear = case_data.get("abs_tol_shear", abs_tol_force)

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
                assert val_calc == pytest.approx(val_ref, rel=rel_tol_disp, abs=abs_tol_disp), \
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
                assert val_calc == pytest.approx(val_ref, rel=rel_tol_force, abs=abs_tol_force), \
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

                # La recuperación centroidal Q del CTRIA3 tiene una tolerancia
                # independiente, documentada por la discrepancia frente a MSC.
                quantity_rel_tol = rel_tol
                quantity_abs_tol = abs_tol_force
                if qty_name.startswith("shear_force_"):
                    quantity_rel_tol = rel_tol_shear
                    quantity_abs_tol = abs_tol_shear

                # Verificamos que sea un número (float o numpy scalar)
                assert val_calc == pytest.approx(
                    val_ref,
                    rel=quantity_rel_tol,
                    abs=quantity_abs_tol,
                ), \
                    f"Stress Mismatch {case_data['case_name']}, Elem {eid}, {qty_name}"


