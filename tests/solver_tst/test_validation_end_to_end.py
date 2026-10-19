"""End-to-end CROD and CBAR validation against MSC Nastran references.
Runs full analysis cycles (Parse -> Assemble -> Solve) on external .dat files
and compares the resulting displacements against reference values extracted 
from MSC Nastran .f06 output files."""
import pytest

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
   # --- Case 1: simple rod ---
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

    # --- Case 2: rods in series with multiple materials ---
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

    # --- Case 3: planar truss in rotated coordinates ---
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

    # --- Case 4: 3D rail effect ---
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

    # --- Case 5: parser input ---
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

    # --- Case 6: oblique torsion ---
    # Validates Moment -> Torque projection on rotated element.
    {
        "case_name": "Ej6_Oblique_Torsion",
        "folder": "case_6",
        "filename": "Ej6.dat",
        "expected_displacements": { 
            2: {3: 4.062500E-01} # Rotation X (radians)
        },
        "expected_element_forces": { 
            # Element 106. Torque is amplified by 1/cos(alpha) = 1/0.8 = 1.25
            106: {
                "torque": 1.250000E+03, 
                "torsional_stress": 0.0 # Zero because C=0 in PROD entry
            } 
        }
    },

    # --- Case 7: full 3D tripod ---
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

# --- Case 8: 3D rail along basic Y ---
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

    # --- Case 9: 3D rail along basic X ---
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

    # --- Case 10: multi-bay planar classroom truss ---
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
        "DAA classroom example with inactive rotations and out-of-plane translation constrained."
    )
    },
    # --- Case 11: heterogeneous spider and torsion rod ---
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

    # --- Case 12: 3D chain with enforced displacements ---


    # --- Case 13: CBAR cantilever under combined bending and axial load ---
    {
        "case_name": "Ej13_CBAR_Cantilever_2D",
        "folder": "case_13",
        "filename": "Ej13.dat",
        "expected_displacements": {
            # Node 2: extension (x), deflection (y/z), and rotation (z).
            2: {0: 1.190476E-03, 1: 1.190481E-01, 2:0.0, 5: 3.571443E-03}
        },
        "expected_reactions": {
            # Node 1 clamp: reactions Fx, Fy, Fz, and Mz.
            1: {0: -5.000000E+02, 1: -5.000000E+02, 2:0.0, 5: -2.500000E+04}
        },
        "expected_element_forces": {
            101: {
                "axial_force": 5.000000E+02,
                "axial_stress": 5.000000E+00,
                # Dotted names address nested result dictionaries.
                # Nastran BEND-MOMENT END-A PLANE 1 maps to local Mz.
                # End A (node 1) reactions.
                "forces_n1.Mz": -2.500000E+04, # Reaction Moment (Clockwise)
                "forces_n1.Fy": -5.000000E+02, # Reaction Shear (Down)
                
                # End B (node 2), free end.
                "forces_n2.Mz": 0.0,           # Zero Moment at tip
                "forces_n2.Fy": 5.000000E+02,  # Shear transmitted
                "forces_n2.Fx": 5.000000E+02   # Axial transmitted
            }
        },
        "description": "CBAR Cantilever with combined Axial and Tip Load (Bending)."
    },

    # --- Case 14: 3D CBAR ---
    {
        "case_name": "Ej14_CBAR_3D",
        "folder": "case_14",
        "filename": "Ej14.dat",
        "expected_displacements": {
            2: {
                0: -2.225381E-05,
                1: 6.369073E-04,
                2: -3.776155E-05,
                3: -1.056354E-04,
                4: 4.194536E-06,
                5: -7.441069E-05
            }
        },
        "expected_reactions": {
            1: {
                0: -1.000000E+03, 1: -5.000000E+03, 2: 2.000000E+03,
                3: 1.100000E+04, 4: 5.000000E+02, 5: 8.000000E+03
            }
        },
        "expected_element_forces": {
            # F06 End A: Mz=1837, My=13478, Fy=492, Fz=4082, torque=452,
            # axial=3618.  With the local convention, F06 plane 1 maps to
            # -forces_n1.Mz and plane 2 maps to forces_n1.My.
            101: {
                "axial_force": 3.618136E+03,
                "axial_stress": 3.618136E+01, 
                "torque": 4.522670E+02,
                
                # Local forces at node 1 use the opposite sign from Nastran end A.
                "forces_n1.Mz": -1.837117E+03, 
                "forces_n1.My":  1.347852E+04,
                "forces_n1.Fy": -4.923660E+02,
                "forces_n1.Fz": -4.082483E+03,

                "forces_n1.Mx": -4.522670E+02  # Torque Reaction
            }
        },
        "description": "3D CBAR with biaxial bending, torsion and axial load. Validates full 12x12 matrix and 3D transformation."
    },

    # --- Case 15: CBAR elements in series with recovery points ---
    # Two CBAR elements aligned with basic Y under combined end loads. Validates
    # PBAR stress points C, D, E, and F.
    {
        "case_name": "Ej15_Series_CBAR_Stress_Points",
        "folder": "case_15",
        "filename": "Ej15.dat",
        "expected_displacements": {
            # Node 2 (central joint)
            2: {
                0:  6.031746E-08,  # T1 (x)
                1:  1.904762E-07,  # T2 (y)
                2:  1.555556E-07,  # T3 (z)
                3:  1.428571E-07,  # R1 (rx)
                4: -1.485714E-07,  # R2 (ry)
                5: -5.714286E-08   # R3 (rz)
            },
            # Node 3 (free tip)
            3: {
                0:  6.603175E-07,
                1:  7.619048E-07,
                2:  1.334127E-06,
                3:  5.000000E-07,
                4: -1.262857E-06,
                5: -3.071429E-07
            }
        },
        "expected_reactions": {
            # Node 1 clamped at the origin.
            1: {
                0: -1.000000E+00, # Fx Reaction (opposes load Fx=1)
                1: -2.000000E+00, # Fy Reaction (opposes load Fy=2)
                2: -4.000000E+00, # Fz Reaction (opposes load Fz=4)
                3: -1.900000E+01, # Mx Reaction (Moment Balance X)
                4:  3.000000E+00, # My Reaction (Moment Balance Y)
                5:  7.000000E+00  # Mz Reaction (Moment Balance Z)
            }
        },
        "expected_element_forces": {
            # Lower element 101.
            101: {
                # Internal Forces at End A (Node 1) from Nastran:
                # Axial=2.0, Torque=-3.0, M_Plane1(Mz)=7.0, M_Plane2(My)=-19.0
                # F06 plane 1 maps to -forces_n1.Mz; plane 2 maps to
                # forces_n1.My in the local CBAR convention.
                "axial_force": 2.000000E+00,
                "axial_stress": 2.000000E-02,
                "torque": -3.000000E+00,
                "forces_n1.Mz": -7.000000E+00, 
                "forces_n1.My": -1.900000E+01,
        
                # Nastran Output Page 14 (STRESSES IN BAR ELEMENTS)
                # PBAR 401 defines points C(10,0) and D(0,10).
                # SA1 = Stress at Point C, End A. Val = -7.00E-02
                # SA2 = Stress at Point D, End A. Val =  1.90E-01
                # Solver order follows the PBAR list: P0=C and P1=D.
                # P0 (point C: y=10.0, z=0.0).
                # Nastran SA1 = -7.00E-02.
                "stress_points.P0.sigma_bending": -7.000000E-02, # Point C
                "stress_points.P0.sigma_axial":    2.000000E-02, # Point C
                # Total = -0.07 + 0.02 = -0.05
                "stress_points.P0.sigma_total":   -5.000000E-02, # Point C

                # P1 (point D: y=0.0, z=10.0).
                # Nastran SA2 = 1.90E-01.
                "stress_points.P1.sigma_bending":  1.900000E-01, # Point D
                "stress_points.P1.sigma_axial":    2.000000E-02, # Point D
                # Total = 0.19 + 0.02 = 0.21
                "stress_points.P1.sigma_total":    2.100000E-01,  # Point D
                
                # MAX and MIN total stress across all recovery points (SA-MAX, SA-MIN)
                "max_total_stress":  2.100000E-01, # Matches SA-MAX
                "min_total_stress": -5.000000E-02, # Matches SA-MIN

                # SB1/SB2 are the end-B bending stresses in the F06.
                "stress_points_end_b.P0.sigma_bending": -5.000000E-02,
                "stress_points_end_b.P0.sigma_total":   -3.000000E-02,
                "stress_points_end_b.P1.sigma_bending":  1.100000E-01,
                "stress_points_end_b.P1.sigma_total":    1.300000E-01,
            },

            # --- ELEMENT 102 ---
            102: {
                # Internal Forces at End A (Node 2) from Nastran:                
                "axial_force": 2.000000E+00,
                "axial_stress": 4.000000E-02,
                "torque": -3.000000E+00,
                "forces_n1.Mz": -5.000000E+00,
                "forces_n1.My": -1.100000E+01,
 
                # --- STRESS POINTS VALIDATION ---
                # PBAR 402 defines points C(5,5) and D(-5,-5).
                # SA1 (point C, end A): 1.50E-01.
                "stress_points.P0.sigma_bending":  1.500000E-01, # Point C
                "stress_points.P0.sigma_axial":   4.000000E-02,  # Point C
                # Total = 0.15 + 0.04 = 0.19
                "stress_points.P0.sigma_total":   1.900000E-01,  # Point C
                
                # SA2 (point D, end A): -1.50E-01.
                "stress_points.P1.sigma_bending": -1.500000E-01, # Point D
                "stress_points.P1.sigma_axial":    4.000000E-02, # Point D
                # Total = -0.15 + 0.04 = -0.11
                "stress_points.P1.sigma_total":   -1.100000E-01, # Point D

                "max_total_stress":  1.900000E-01, # Matches SA-MAX
                "min_total_stress": -1.100000E-01, # Matches SA-MIN

                "stress_points_end_b.P0.sigma_bending": -7.500000E-02,
                "stress_points_end_b.P0.sigma_total":   -3.500000E-02,
                "stress_points_end_b.P1.sigma_bending":  7.500000E-02,
                "stress_points_end_b.P1.sigma_total":    1.150000E-01,
            }
        },
        "description": "Two CBARs in series (Y-axis). Validates Stress Recovery at specific points (C, D) defined in PBAR."
    },

# --- Case 16: five-member 3D frame ---
    # Statically indeterminate spatial frame with complex loads and orientations;
    # validates the complete 12x12 stiffness in arbitrary orientations.
    {
        "case_name": "Ej16_3D_Frame_5_Bars",
        "folder": "case_16",
        "filename": "Ej16.dat",
        "expected_displacements": {
            # Node 3 translations at the upper-left joint.
            3: {
                0:  1.993276E-04, 
                1:  3.672660E-05,
                2: -1.828093E-04,
                3: -1.165399E-04, # R1 (Rx)
                5: -1.681146E-04  # R3 (Rz)
            },
            # Node 5 translations at the right tip.
            5: {
                0:  3.721953E-04,
                1:  1.723244E-05,
                2: -2.991802E-04,
                4: -5.335918E-05  # R2 (Ry)
            }
        },
        "expected_reactions": {
            # Node 1 (Fixed)
            1: {
                0: -4.171473E+03, # Fx
                1: -1.052408E+03, # Fy
                2:  3.697160E+03, # Fz
                3:  8.799405E+03, # Mx
                4: -3.292138E+03, # My
                5:  9.696733E+03  # Mz
            },
            # Node 2 (Fixed)
            2: {
                0:  1.714733E+02, # Fx
                1:  4.052408E+03, # Fy
                4: -5.313543E+03  # My
            }
        },
        "expected_element_forces": {
            # Element 101 (member 1-3).
            # Nastran End A (Internal): Axial=-229.6, Torque=-364.5, Mz=4880.7, My=4690.8
            101: {
                "axial_force": -2.296478E+02,
                "torque": -3.644956E+02,
                "forces_n1.Mz": -4.880671E+03, # Reaction (-) vs Internal (+)
                "forces_n1.My": 4.690753E+03
            },
            
            # Element 102 (member 2-4).
            # Nastran End A: Axial=1770.9, Mz=-7646.0, My=-5313.5
            102: {
                "axial_force": 1.770918E+03,
                "forces_n1.Mz": 7.645998E+03,  # Reaction (+) vs Internal (-)
                "forces_n1.My": -5.313543E+03
            },

            # Element 105 (member 4-5); additional member coverage.
            # Nastran End A: Axial=600.1, Torque=-338.1, Mz=-3703.4, My=-457.49
            105: {
                "axial_force": 6.001212E+02,
                "torque": -3.381103E+02,
                "forces_n1.Mz": 3.703380E+03,
                "forces_n1.My":-4.574928E+02
            }
        },
        "description": "3D Frame with 5 CBARs. Validates full 6-DOF coupling and arbitrary orientation."
    },

    # --- Case 17: hybrid tripod (CBAR + CROD) ---
    # Three CBAR legs support a central CROD post under vertical load and
    # torsion at the joint.
    {
        "case_name": "Ej17_Hybrid_Tripod_CBAR_CROD",
        "folder": "case_17",
        "filename": "Ej17.dat",
        "expected_displacements": {
            1: {
                0: 0.0, # T1
                1: 0.0, # T2
                2: 0.0, # T3 
                3: 0.0, # R1
                4: 0.0, # R2
                5: 0.0  # R3 
            },
            # Node 4: central joint.
            4: {
                0:  1.786003E-04, # T1
                1: -9.267777E-05, # T2
                2:  5.951593E-04, # T3 (vertical translation)
                3:  3.652847E-05, # R1
                4:  1.530955E-04, # R2
                5:  1.264375E-04  # R3 (Z rotation - assembly torsion)
            },
            # Node 5: top of CROD; follows node 4 plus its own elongation.
            5: {
                0:  0.0,          # No lateral load, so Nastran reports absolute zero in T1/T2.
                2:  2.023731E-03, # T3: vertical displacement under the direct load.
                4:  0.0,          # R2 
                5:  1.264375E-04  # R3 follows node 4 as CROD torsional rigid motion.
            }
        },
        "expected_reactions": {            
            # Supports at 1, 2, 3.
            # Total Fz Reaction: -1393.1 -531.2 -75.6 ~= -2000 (Matches Load)
            1: { 
                0: -2.768665E+02,
                1:  3.374797E+02,
                2: -1.393156E+03,
                3: -2.676242E+03, # R1
                4: -2.106779E+03, # R2
                5:  4.265762E+01  # R3 
            },
            2: { 2: -5.311886E+02 },
            3: { 2: -7.565560E+01,
                 3:  1.832815E+02,
                 4: -7.837672E-01, 
                 5:  1.058018E+02
            }
        },
        "expected_element_forces": {
            # Element 104, vertical CROD post.
            # Transfers Fz=2000 as pure axial force.
            # Torque is negligible (7.10E-15 in the F06) for an ideal CROD.
            104: {
                "type": "CROD",
                "axial_force": 2.000000E+03,
                "torque": 0.0
            },
            
            # Element 101, CBAR leg 1.
            # Carries axial load and torque.
            # Nastran end A (ground): axial=1008.3, torque=298.2, Mz=-2586.0, My=-2196.8.
            101: {
                "type": "CBAR",
                "axial_force":   1.008341E+03,
                "torque":        2.982355E+02,
                "forces_n1.Mz":  2.585990E+03, # Reaction (+)
                "forces_n1.My": -2.196889E+03
            },

            # Element 103, the least-loaded CBAR leg.
            # Nastran End A: Axial=81.3, Mz=211.2, My=10.2
            103: {
                "axial_force":   8.132725E+01,
                "torque":       -7.479646E+00,
                "forces_n1.Mz": -2.112480E+02, # Reaction (-)
                "forces_n1.My":  1.025171E+01
            }
        },
        "description": "Hybrid assembly (CBAR legs + CROD post). Validates coupling of different element types and stiffness matrices."
    }, 

    # --- Case 18: hybrid tripod with all C-F stress points ---
    # Variant of case 17 validating corner stresses at all four PBAR points.
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
            # Node 4 (central joint) - values from ej18.f06 (page 11)
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
                0:  0.0,           # T1
                1:  0.0,           # T2
                2:  2.023731E-03,  # T3 (dominant axial elongation)
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
            # Element 101, CBAR leg 1.
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
                # Total = 15.64915 + 84.67422 = 100.32337
                "stress_points.P0.sigma_total":    1.003234E+02, # Total Point C (=SA-MAX)

                # Check Total (Axial + Bending) for Point E
                # Total = 15.64915 - 23.27506 = -7.62591
                "stress_points.P2.sigma_total":    -7.62591E+00, # Total Point E

                # MAX and MIN total stress across all recovery points (SA-MAX, SA-MIN)
                "max_total_stress":  1.003234E+02, # Matches SA-MAX column
                "min_total_stress": -6.902507E+01  # Matches SA-MIN column
            },
            
            # Element 103, CBAR leg 3.
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

    # --- Case 19: irregular 3D helix ---
    # Eleven elements, three PBAR properties, and two materials in a spiral
    # geometry not aligned with the basic axes. Combined tip force and moment
    # validate complex rotations and SA-MAX/MIN recovery points.
    {
        "case_name": "Ej19_Complex_Helix_Stress_Recovery_Points",
        "folder": "case_19",
        "filename": "Ej19.dat",
        "description": "Complex 3D helix with mixed materials and variable cross-sections. Validates general 3D coordinate transformations and stress recovery under asymmetric loading.",
        "expected_displacements": {

            1:{0:0.0, 1:0.0, 4:0.0, 5:0.0},

            7: {0: -5.443371E-05,
                1: -1.857526E-04,
                2: -2.951618E-04,
                3: -5.901665E-07,
                4: 2.317618E-04,
                5: -1.219883E-04
            },
            # Node 12 (Tip)
            12: {
                0:  1.347035E-02, # T1
                1:  6.321326E-03, # T2
                2: -1.308450E-02, # T3
                3: -3.590305E-04, # R1
                4:  2.037043E-04, # R2
                5: -3.121635E-04  # R3
            }
        },
        "expected_reactions": {
            # Node 1 (clamped base)
            1: {
                0: -1.000000E+02,
                1:  2.000000E+02,
                2:  5.000000E+02,
                3: -4.500000E+03,
                4: -5.500000E+03,
                5:  1.500000E+03
            }
        },
        "expected_element_forces": {
            # --- ELEMENT 101 (Base - Heavy Steel) ---
            101: {
                "axial_force": -1.341641E+02,
                "axial_stress": -5.366563E-02,
                "torque":        3.354102E+03,
                
                # Stress Points (PBAR 401: Symmetric Section)
                # SA Columns = Pure Bending in this case
                "stress_points.P0.sigma_bending":  1.031682E-01, # SA1
                "stress_points.P1.sigma_bending": -4.256780E-01, # SA2
                "stress_points.P2.sigma_bending": -1.031682E-01, # SA3
                "stress_points.P3.sigma_bending":  4.256780E-01, # SA4
                
                # MAX and MIN total stress across all recovery points (SA-MAX, SA-MIN)
                # Max = Axial(-0.053) + SA4(+0.425) = 0.372
                "max_total_stress":  3.720123E-01, # Matches SA-MAX
                "min_total_stress": -4.793436E-01  # Matches SA-MIN
            },
            
            # --- ELEMENT 106 (Mid - Light Steel - Rotated Y) ---
            # PBAR 403: Asymmetric Section
            106: {
                "axial_stress": -8.944272E-02,
                "torque":       -8.720665E+03,

                # SA Columns (Bending Only)
                "stress_points.P0.sigma_bending": -3.881966E-01, # SA1
                "stress_points.P1.sigma_bending":  6.683282E-01, # SA2
                "stress_points.P2.sigma_bending":  1.381966E-01, # SA3   
                "stress_points.P3.sigma_bending": -4.447214E-01, # SA4
                  
                
                # MAX and MIN total stress across all recovery points (SA-MAX, SA-MIN)
                # Max = Axial(-0.089) + SA2(+0.668) = 0.579
                "max_total_stress":  5.788854E-01, # Matches SA-MAX
                # Min = Axial(-0.089) + SA4(-0.444) = -0.534
                "min_total_stress": -5.341641E-01  # Matches SA-MIN
            },
            
            # --- ELEMENT 111 (Tip - Mixed Props) ---
            # High bending stresses close to load application
            111: {
                "axial_stress": -5.059644E-01,
                "torque":       -9.486833E+02,
                
                "stress_points.P0.sigma_bending":  7.869626E+00, # SA1
                "stress_points.P1.sigma_bending": -3.512897E+00, # SA2
                "stress_points.P2.sigma_bending": -4.369626E+00, # SA3
                "stress_points.P3.sigma_bending":  5.252149E+00, # SA4
                
                # MAX and MIN total stress across all recovery points (SA-MAX, SA-MIN)
                # Max = Axial(-0.505) + SA1(+7.869) = 7.363
                "max_total_stress":  7.363662E+00, # Matches SA-MAX
                # Min = Axial(-0.505) + SA3(-4.369) = -4.875
                "min_total_stress": -4.875591E+00  # Matches SA-MIN
            }
        }
    }

]


class TestNastranEndToEnd:
    """Run the parametrized CROD and CBAR reference cases."""

    @pytest.mark.parametrize("case_data", VALIDATION_DATA, ids=lambda x: x["case_name"])
    def test_run_nastran_case(self, case_data, validation_data_dir):
        """Execute and validate one Nastran reference case."""
        folder = case_data.get("folder", "")
        filename = case_data["filename"]
        
        file_path = validation_data_dir / folder / filename
        
        if not file_path.exists():
            pytest.fail(f"Test file not found: {file_path}")

        try:
            model = parse_nastran_input(str(file_path))
        except Exception as exc:
            pytest.fail(f"Parser failed on {case_data['filename']}: {exc}")

        try:
            solver = LinearStaticSolver(model)
            solver.assemble()
            results = solver.solve()
        except Exception as exc:
            pytest.fail(f"Solver failed on {case_data['filename']}: {exc}")

        # Unpack results
        displacements = results["displacements"]
        reactions = results["reactions"]
        elem_stresses = results["stresses"]
        
        rel_tol = case_data.get("rel_tol", 1e-3)
        abs_tol_disp = case_data.get("abs_tol_disp", 1e-9)
        abs_tol_force = case_data.get("abs_tol_force", 1e-4)

        expected_disp = case_data.get("expected_displacements", {})
        
        for nid, dofs_dict in expected_disp.items():
            try:
                global_indices = solver.dof_mapper.get_dof_indices(nid)
            except ValueError:
                pytest.fail(f"Node {nid} expected in results but not found in Solver Mapper.")
            
            for local_dof, val_ref in dofs_dict.items():
                idx = global_indices[local_dof]
                val_calc = displacements[idx]
                
                assert val_calc == pytest.approx(val_ref, rel=rel_tol, abs=abs_tol_disp), \
                    f"Mismatch in {case_data['case_name']}, Node {nid}, DOF {local_dof}. " \
                    f"Expected {val_ref}, Got {val_calc}"

        
        expected_react = case_data.get("expected_reactions", {})
        
        for nid, dofs_dict in expected_react.items():
            if nid not in reactions:
                # Nastran reports reactions only for constrained nodes.
                pytest.fail(f"Reaction for Node {nid} expected but not found in results.")
            
            reac_vec = reactions[nid] 
            
            for local_dof, val_ref in dofs_dict.items():
                val_calc = reac_vec[local_dof]
                
                # 0.1 % relative or 1e-4 absolute near zero.
                assert val_calc == pytest.approx(val_ref, rel=rel_tol, abs=abs_tol_force), \
                    f"Reaction Mismatch {case_data['case_name']}, Node {nid}, DOF {local_dof}"

        
        expected_forces = case_data.get("expected_element_forces", {})
        
        for eid, quantities in expected_forces.items():
            if eid not in elem_stresses:
                pytest.fail(f"Stresses for Element {eid} expected but not found in results.")
            
            calc_data = elem_stresses[eid]
            
            for qty_name, val_ref in quantities.items():
             
                val_calc = calc_data
                try:
                    if "." in qty_name:
                        keys = qty_name.split(".")
                        for k in keys:
                            # Reject a malformed scalar before nested lookup.
                            if not hasattr(val_calc, '__getitem__'):
                                raise ValueError(f"Cannot access key '{k}' in scalar value: {val_calc}")
                            val_calc = val_calc[k]
                    else:
                        val_calc = calc_data[qty_name]
                except (KeyError, IndexError, ValueError, TypeError) as exc:
                    pytest.fail(f"Error accessing '{qty_name}' in Element {eid}. Details: {exc}")

                # pytest.approx accepts Python and NumPy scalar values.
                assert val_calc == pytest.approx(val_ref, rel=rel_tol, abs=abs_tol_force), \
                    f"Stress Mismatch {case_data['case_name']}, Elem {eid}, {qty_name}"
