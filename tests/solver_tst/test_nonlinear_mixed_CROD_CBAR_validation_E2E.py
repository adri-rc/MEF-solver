"""Validación final SOL 400 de modelos no lineales mixtos CROD/CBAR.

Este módulo mantiene separada la evidencia integrada de las campañas puras:

* ``test_nonlinear_validation_end_to_end.py`` valida CROD;
* ``test_nonlinear_CBAR_validation_E2E_.py`` valida CBAR;
* este archivo valida el ensamblaje y equilibrio cuando ambos conviven.

Las referencias se transcriben del último bloque ``LOAD STEP = 1.00000E+00``
del F06/PCH. El test no lee resultados de Nastran durante su ejecución.

Claves de cada entrada de ``NONLINEAR_MIXED_VALIDATION_DATA``:

``case_name``
    Identificador único usado también como ID de pytest.
``folder`` y ``filename``
    Carpeta dentro de ``validation_cases`` y DAT exacto que se resuelve.
``description`` y ``note``
    Comportamiento validado y, opcionalmente, una observación importante.
``rel_tol`` y ``abs_tol``
    Tolerancias relativa y absoluta comunes a las magnitudes del caso.
``expected_displacements`` y ``expected_reactions``
    ``{grid_id: {dof_index: value}}``; 0=T1, 1=T2, 2=T3, 3=R1, 4=R2,
    5=R3. Los valores proceden de ``DISPLACEMENT VECTOR`` y
    ``FORCES OF SINGLE-POINT CONSTRAINT``.
``expected_element_types``
    ``{element_id: "CBAR" | "CROD"}``; demuestra que se recuperan ambos tipos.
``expected_element_results``
    ``{element_id: {result_name: value}}``. Se comparan tensión, deformación y
    fuerza axial, además de la longitud de la cuerda deformada calculada con
    el ``DISPLACEMENT VECTOR``. Para CROD, ``torsional_stress`` procede de la
    columna ``LIN. TORSIONAL STRESS``. Para CBAR, ``axial_stress`` es la
    ordenada en el centroide del plano C--F y las rutas anidadas, por ejemplo
    ``stress_points_end_a.P0.sigma_total``, corresponden a C/D/E/F = P0/P1/
    P2/P3. SOL 400 no imprime en esta tabla una tensión torsional CBAR aislada,
    por lo que no se inventa una referencia para ella.
``expected_strain_energy`` y ``expected_increments``
    Energía de ``GLOBAL WORKS & ENERGIES / TOT STRAIN ENERGY`` y número de
    incrementos convergidos solicitado por ``NLSTEP,FIXED``.

Ejemplo para añadir otro modelo mixto::

    {
        "case_name": "case_NNN_mixed_model",
        "folder": "case_NNN",
        "filename": "ejNNN_mixed_model.dat",
        "description": "Physical behavior validated by the model.",
        "rel_tol": 1.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {10: {0: 1.234, 4: 0.005}},
        "expected_reactions": {1: {0: -100000.0}},
        "expected_element_types": {1: "CBAR", 20: "CROD"},
        "expected_element_results": {
            1: {
                "axial_stress": 12.5,
                "axial_strain": 5.95e-5,
                "axial_force": 15000.0,
                "current_length": 1000.0595,
                "stress_points_end_a.P0.sigma_total": -40.0,
            },
            20: {
                "axial_stress": -30.0,
                "axial_strain": -1.43e-4,
                "axial_force": -13500.0,
                "current_length": 999.857,
                "torsional_stress": 0.0,
            },
        },
        "expected_strain_energy": 1234.5,
        "expected_increments": 40,
    }
"""

from pathlib import Path

import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.nonlinear_analysis import NonlinearStaticSolver


VALIDATION_ROOT = Path(__file__).parent / "validation_cases"

NONLINEAR_MIXED_VALIDATION_DATA = [
    {
        "case_name": "case_167_cantilever_with_crod_brace",
        "folder": "case_167",
        "filename": "ej167_mixed_cbar_crod_cantilever_brace_2_elements.dat",
        "description": "One CBAR cantilever sharing a spatial load with one CROD brace.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            2: {
                0: 0.2792798,
                1: 2.247848,
                2: 1.740790,
                3: 0.02974886,
                4: -0.005986621,
                5: 0.007268436,
            },
        },
        "expected_reactions": {
            1: {
                0: -120079.8,
                1: -7836.318,
                2: 3819.395,
                3: -7.977773e6,
                4: 9.498262e6,
                5: -2.243406e7,
            },
            3: {0: 79.76550, 1: -142163.7, 2: 86180.60, 3: 0.0, 4: 0.0, 5: 0.0},
        },
        "expected_element_types": {1: "CBAR", 2: "CROD"},
        "expected_element_results": {
            1: {
                "axial_stress": 100.08058,
                "axial_strain": 4.765741e-4,
                "axial_force": 120096.696,
                "current_length": 600.286013,
                "stress_points_end_a.P2.sigma_total": 269.8760,
                "stress_points_end_b.P0.sigma_total": -75.60275,
                "stress_points_end_b.P2.sigma_total": 275.7639,
            },
            2: {
                "axial_stress": -369.4348,
                "axial_strain": -1.759213e-3,
                "axial_force": -166245.66,
                "current_length": 582.069401,
                "torsional_stress": 0.0,
            },
        },
        "expected_strain_energy": 327358.6,
        "expected_increments": 25,
    },
    {
        "case_name": "case_168_braced_portal",
        "folder": "case_168",
        "filename": "ej168_mixed_cbar_crod_portal_4_elements.dat",
        "description": "Three-CBAR portal with a CROD diagonal and two materials.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            3: {
                0: 2.897509,
                1: 3.568657,
                2: 0.06088335,
                3: -0.006533808,
                4: 0.002397578,
                5: -0.006077024,
            },
            4: {
                0: 3.008362,
                1: -1.093099,
                2: -0.5938926,
                3: 0.002056600,
                4: 0.003733611,
                5: -0.005083175,
            },
        },
        "expected_reactions": {
            1: {
                0: -125105.1,
                1: -66129.57,
                2: -79082.60,
                3: 5.040070e7,
                4: -3.560469e7,
                5: 1.818815e6,
            },
            2: {
                0: -58894.85,
                1: 18129.57,
                2: 247082.6,
                3: -1.497400e7,
                4: -3.173664e7,
                5: 1.552748e6,
            },
        },
        "expected_element_types": {1: "CBAR", 2: "CBAR", 3: "CBAR", 4: "CROD"},
        "expected_element_results": {
            1: {
                "axial_stress": 19.44746,
                "axial_strain": 9.2607025e-5,
                "axial_force": 31115.936,
                "current_length": 800.074089,
            },
            2: {
                "axial_stress": -154.27135,
                "axial_strain": -7.34625475e-4,
                "axial_force": -246834.16,
                "current_length": 799.412515,
                "stress_points_end_a.P0.sigma_total": 101.3531,
                "stress_points_end_a.P2.sigma_total": -409.8958,
                "stress_points_end_b.P0.sigma_total": -264.3458,
                "stress_points_end_b.P2.sigma_total": -44.19694,
            },
            3: {
                "axial_stress": 32.73145,
                "axial_strain": 1.55864125e-4,
                "axial_force": 52370.32,
                "current_length": 800.124701,
            },
            4: {
                "axial_stress": 105.8428,
                "axial_strain": 1.512040e-3,
                "axial_force": 68797.82,
                "current_length": 1133.08153,
                "torsional_stress": 0.0,
            },
        },
        "expected_strain_energy": 484761.5,
        "expected_increments": 30,
    },
    {
        "case_name": "case_169_guyed_mast",
        "folder": "case_169",
        "filename": "ej169_mixed_cbar_crod_guyed_mast_5_elements.dat",
        "description": "Two-CBAR mast stabilized by three CROD guy elements.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            2: {
                0: 7.184081,
                1: -6.178676,
                2: -0.2606053,
                3: 0.009185985,
                4: 0.01029419,
                5: -0.02953908,
            },
            6: {
                0: 2.068790,
                1: -1.686877,
                2: -0.1141305,
                3: 0.005356068,
                4: 0.006470393,
                5: -0.01476831,
            },
        },
        "expected_reactions": {
            1: {
                0: -27916.48,
                1: 15179.05,
                2: 68306.42,
                3: -3.863295e7,
                4: -4.898322e7,
                5: 9.936561e6,
            },
            3: {0: -167647.0, 1: -1306.528, 2: 253694.1, 3: 0.0, 4: 0.0, 5: 0.0},
        },
        "expected_element_types": {
            1: "CBAR",
            2: "CBAR",
            3: "CROD",
            4: "CROD",
            5: "CROD",
        },
        "expected_element_results": {
            1: {
                "axial_stress": -37.8704425,
                "axial_strain": -1.80335377e-4,
                "axial_force": -68166.7965,
                "current_length": 599.891808,
                "stress_points_end_a.P1.sigma_total": 181.1859,
                "stress_points_end_a.P3.sigma_total": -256.9268,
                "stress_points_end_b.P1.sigma_total": 115.9135,
                "stress_points_end_b.P3.sigma_total": -191.6544,
            },
            2: {
                "axial_stress": -37.750165,
                "axial_strain": -1.79762725e-4,
                "axial_force": -67950.297,
                "current_length": 599.892152,
            },
            3: {
                "axial_stress": -608.1707,
                "axial_strain": -2.896051e-3,
                "axial_force": -304085.35,
                "current_length": 1438.04377,
                "torsional_stress": 0.0,
            },
            4: {
                "axial_stress": 231.7520,
                "axial_strain": 3.310743e-3,
                "axial_force": 173814.0,
                "current_length": 1450.46951,
                "torsional_stress": 0.0,
            },
            5: {
                "axial_stress": -172.8177,
                "axial_strain": -8.229413e-4,
                "axial_force": -86408.85,
                "current_length": 1444.49352,
                "torsional_stress": 0.0,
            },
        },
        "expected_strain_energy": 1576299.0,
        "expected_increments": 36,
    },
    {
        "case_name": "case_170_two_bay_braced_frame",
        "folder": "case_170",
        "filename": "ej170_mixed_cbar_crod_two_bay_frame_9_elements.dat",
        "description": "Two-bay frame with five CBARs and four CROD braces.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            4: {
                0: 0.8742845,
                1: 1.040598,
                2: -0.06962341,
                3: -0.002223749,
                4: 0.001507518,
                5: -0.0008895377,
            },
            6: {
                0: 1.593592,
                1: -0.7219370,
                2: -0.3035130,
                3: 0.001705041,
                4: 0.004046119,
                5: 0.0006903453,
            },
        },
        "expected_reactions": {
            1: {
                0: -71539.52,
                1: -34692.40,
                2: -16715.93,
                3: 2.416012e7,
                4: -1.344255e7,
                5: 400507.0,
            },
            3: {
                0: -83096.06,
                1: 16501.32,
                2: 245029.6,
                3: -1.502720e7,
                4: -1.672394e7,
                5: -301653.0,
            },
        },
        "expected_element_types": {
            1: "CBAR",
            3: "CBAR",
            5: "CBAR",
            6: "CROD",
            9: "CROD",
        },
        "expected_element_results": {
            1: {
                "axial_stress": -20.4921508,
                "axial_strain": -9.75816625e-5,
                "axial_force": -36885.8714,
                "current_length": 699.931696,
            },
            3: {
                "axial_stress": -90.417205,
                "axial_strain": -4.305580e-4,
                "axial_force": -162750.969,
                "current_length": 699.698674,
                "stress_points_end_a.P2.sigma_total": -224.1551,
                "stress_points_end_a.P3.sigma_total": -157.3707,
                "stress_points_end_b.P2.sigma_total": -143.3205,
                "stress_points_end_b.P3.sigma_total": -128.2868,
            },
            5: {
                "axial_stress": 57.929657,
                "axial_strain": 8.27566637e-4,
                "axial_force": 75308.5541,
                "current_length": 600.496746,
            },
            6: {
                "axial_stress": 128.6743,
                "axial_strain": 6.127346e-4,
                "axial_force": 70770.865,
                "current_length": 922.519359,
                "torsional_stress": 0.0,
            },
            9: {
                "axial_stress": -196.7913,
                "axial_strain": -9.371014e-4,
                "axial_force": -108235.215,
                "current_length": 921.090481,
                "torsional_stress": 0.0,
            },
        },
        "expected_strain_energy": 276615.0,
        "expected_increments": 40,
    },
    {
        "case_name": "case_171_space_platform",
        "folder": "case_171",
        "filename": "ej171_mixed_cbar_crod_space_platform_14_elements.dat",
        "description": "Spatial platform with CBAR columns/ring and CROD diagonals.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            5: {
                0: 1.080169,
                1: 0.4817405,
                2: -0.05503363,
                3: -0.0006082629,
                4: 0.001667863,
                5: 0.0005911825,
            },
            7: {
                0: 0.8907765,
                1: 0.8402616,
                2: -0.3481199,
                3: -0.0005533090,
                4: 0.002187925,
                5: -0.0005430446,
            },
        },
        "expected_reactions": {
            1: {
                0: -90733.88,
                1: -15277.24,
                2: -16027.77,
                3: 8.181880e6,
                4: -2.752697e7,
                5: -292960.7,
            },
            3: {
                0: -33235.76,
                1: -39436.74,
                2: 223209.0,
                3: 1.774214e7,
                4: -1.332804e7,
                5: 268008.7,
            },
        },
        "expected_element_types": {
            1: "CBAR",
            4: "CBAR",
            7: "CBAR",
            8: "CBAR",
            9: "CROD",
            11: "CROD",
            14: "CROD",
        },
        "expected_element_results": {
            1: {
                "axial_stress": -14.2172995,
                "axial_strain": -6.7701405e-5,
                "axial_force": -29856.3289,
                "current_length": 799.945841,
            },
            4: {
                "axial_stress": -31.3442515,
                "axial_strain": -1.49258338e-4,
                "axial_force": -65822.9282,
                "current_length": 799.880602,
            },
            7: {
                "axial_stress": 38.5316555,
                "axial_strain": 5.50452233e-4,
                "axial_force": 57797.4832,
                "current_length": 800.440483,
                "stress_points_end_a.P0.sigma_total": 72.77356,
                "stress_points_end_a.P2.sigma_total": 4.289752,
                "stress_points_end_b.P1.sigma_total": 55.54436,
                "stress_points_end_b.P3.sigma_total": 21.51895,
            },
            8: {
                "axial_stress": 5.79016937,
                "axial_strain": 8.27167305e-5,
                "axial_force": 8685.25405,
                "current_length": 600.049632,
            },
            9: {
                "axial_stress": 108.4062,
                "axial_strain": 5.162200e-4,
                "axial_force": 65043.72,
                "current_length": 1131.95489,
                "torsional_stress": 0.0,
            },
            11: {
                "axial_stress": -74.79624,
                "axial_strain": -3.561726e-4,
                "axial_force": -44877.744,
                "current_length": 1130.96789,
                "torsional_stress": 0.0,
            },
            14: {
                "axial_stress": 31.79691,
                "axial_strain": 4.542415e-4,
                "axial_force": 27027.3735,
                "current_length": 1000.45424,
                "torsional_stress": 0.0,
            },
        },
        "expected_strain_energy": 215675.7,
        "expected_increments": 40,
    },
    {
        "case_name": "case_172_two_story_frame_with_settlement",
        "folder": "case_172",
        "filename": "ej172_mixed_cbar_crod_two_story_frame_18_elements.dat",
        "description": "Two-story braced frame with a prescribed support settlement.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            3: {0: 0.0, 1: 0.0, 2: -3.0, 3: 0.0, 4: 0.0, 5: 0.0},
            7: {
                0: 4.577620,
                1: 1.881314,
                2: 0.3312835,
                3: -0.002190306,
                4: 0.004016811,
                5: -0.001588982,
            },
            9: {
                0: 5.311529,
                1: -2.443476,
                2: -3.557384,
                3: 0.002870390,
                4: 0.004677373,
                5: -0.005922480,
            },
        },
        "expected_reactions": {
            1: {
                0: -198075.1,
                1: -29089.49,
                2: -335367.9,
                3: 2.378532e7,
                4: -4.991823e7,
                5: 380846.8,
            },
            3: {
                0: -81256.22,
                1: 33586.03,
                2: 216415.0,
                3: -2.973619e7,
                4: -5.953192e7,
                5: 1.715209e6,
            },
        },
        "expected_element_types": {
            1: "CBAR",
            5: "CBAR",
            9: "CBAR",
            10: "CBAR",
            11: "CROD",
            15: "CROD",
            18: "CROD",
        },
        "expected_element_results": {
            1: {
                "axial_stress": 105.926595,
                "axial_strain": 5.04412425e-4,
                "axial_force": 233038.509,
                "current_length": 600.302724,
            },
            5: {
                "axial_stress": -136.54475,
                "axial_strain": -6.5021305e-4,
                "axial_force": -300398.45,
                "current_length": 599.609999,
            },
            9: {
                "axial_stress": 46.492645,
                "axial_strain": 6.641808e-4,
                "axial_force": 69738.9675,
                "current_length": 600.398641,
                "stress_points_end_a.P0.sigma_total": -26.05895,
                "stress_points_end_a.P3.sigma_total": 123.8873,
                "stress_points_end_b.P0.sigma_total": 143.4092,
                "stress_points_end_b.P2.sigma_total": -50.42388,
            },
            10: {
                "axial_stress": 40.9746425,
                "axial_strain": 5.85352125e-4,
                "axial_force": 61461.9638,
                "current_length": 600.351314,
            },
            11: {
                "axial_stress": 223.8548,
                "axial_strain": 1.065975e-3,
                "axial_force": 145505.62,
                "current_length": 849.432647,
                "torsional_stress": 0.0,
            },
            15: {
                "axial_stress": 295.9776,
                "axial_strain": 1.409417e-3,
                "axial_force": 192385.44,
                "current_length": 849.724068,
                "torsional_stress": 0.0,
            },
            18: {
                "axial_stress": -75.87540,
                "axial_strain": -3.613114e-4,
                "axial_force": -49319.01,
                "current_length": 848.221555,
                "torsional_stress": 0.0,
            },
        },
        "expected_strain_energy": 1159633.0,
        "expected_increments": 45,
    },
    {
        "case_name": "case_173_pyramid_roof",
        "folder": "case_173",
        "filename": "ej173_mixed_cbar_crod_pyramid_roof_22_elements.dat",
        "description": "Pyramid roof combining LC/LS CBAR sections and ten CROD braces.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            6: {
                0: 0.5666284,
                1: 0.3554892,
                2: -0.2657467,
                3: -0.0007399859,
                4: 0.0002914210,
                5: 8.068994e-6,
            },
            9: {
                0: 1.261114,
                1: 0.8338564,
                2: -1.126088,
                3: -0.0003305095,
                4: 0.002226885,
                5: -1.908803e-5,
            },
        },
        "expected_reactions": {
            1: {
                0: -61639.72,
                1: -8978.979,
                2: 8232.639,
                3: 3.856252e6,
                4: -1.849234e7,
                5: -231031.4,
            },
            3: {
                0: -51507.72,
                1: -76524.40,
                2: 305611.1,
                3: 1.475533e7,
                4: -1.329567e7,
                5: -182790.3,
            },
        },
        "expected_element_types": {
            1: "CBAR",
            6: "CBAR",
            12: "CBAR",
            13: "CROD",
            16: "CROD",
            22: "CROD",
        },
        "expected_element_results": {
            1: {
                "axial_stress": -12.0328835,
                "axial_strain": -5.7299445e-5,
                "axial_force": -28878.9204,
                "current_length": 799.954162,
            },
            6: {
                "axial_stress": 37.036851,
                "axial_strain": 1.7636598e-4,
                "axial_force": 88888.4424,
                "current_length": 900.158743,
            },
            12: {
                "axial_stress": -9.892275,
                "axial_strain": -1.413182e-4,
                "axial_force": -16816.8675,
                "current_length": 901.260445,
                "stress_points_end_a.P1.sigma_total": 64.46609,
                "stress_points_end_a.P3.sigma_total": -84.25064,
                "stress_points_end_b.P1.sigma_total": -75.61488,
                "stress_points_end_b.P3.sigma_total": 55.83033,
            },
            13: {
                "axial_stress": 47.20589,
                "axial_strain": 2.247900e-4,
                "axial_force": 33044.123,
                "current_length": 1442.54471,
                "torsional_stress": 0.0,
            },
            16: {
                "axial_stress": -77.10278,
                "axial_strain": -3.671561e-4,
                "axial_force": -53971.946,
                "current_length": 1203.71734,
                "torsional_stress": 0.0,
            },
            22: {
                "axial_stress": 10.16300,
                "axial_strain": 1.451858e-4,
                "axial_force": 9654.85,
                "current_length": 1500.21778,
                "torsional_stress": 0.0,
            },
        },
        "expected_strain_energy": 357568.3,
        "expected_increments": 50,
    },
    {
        "case_name": "case_174_two_level_space_frame",
        "folder": "case_174",
        "filename": "ej174_mixed_cbar_crod_two_level_space_frame_24_elements.dat",
        "description": "Two-level space frame with rotated loads and prescribed top sway.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            10: {
                0: 14.70629,
                1: 6.953951,
                2: -0.8135936,
                3: -0.004083529,
                4: 0.01227148,
                5: 0.008460244,
            },
            12: {
                0: 12.0,
                1: 2.577555,
                2: -0.1608321,
                3: -0.001178493,
                4: 0.01000640,
                5: 0.003010375,
            },
        },
        "expected_reactions": {
            1: {
                0: -249614.5,
                1: -16916.42,
                2: -235603.8,
                3: 1.418426e7,
                4: -1.154050e8,
                5: -745721.5,
            },
            3: {
                0: -198763.7,
                1: -35017.57,
                2: 577396.2,
                3: 3.512429e7,
                4: -8.972543e7,
                5: -429238.6,
            },
            12: {0: 298579.5, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
        },
        "expected_element_types": {
            1: "CBAR",
            8: "CBAR",
            13: "CBAR",
            16: "CBAR",
            17: "CROD",
            22: "CROD",
            24: "CROD",
        },
        "expected_element_results": {
            1: {
                "axial_stress": 57.6666925,
                "axial_strain": 2.7460325e-4,
                "axial_force": 149933.40,
                "current_length": 800.219713,
            },
            8: {
                "axial_stress": -31.698595,
                "axial_strain": -1.509456e-4,
                "axial_force": -82416.347,
                "current_length": 799.879253,
            },
            13: {
                "axial_stress": 46.96415,
                "axial_strain": 4.083840e-4,
                "axial_force": 70446.225,
                "current_length": 1200.49017,
                "stress_points_end_a.P1.sigma_total": 222.9963,
                "stress_points_end_a.P2.sigma_total": 244.9143,
                "stress_points_end_b.P0.sigma_total": 321.7996,
                "stress_points_end_b.P2.sigma_total": -227.8713,
            },
            16: {
                "axial_stress": -19.4808742,
                "axial_strain": -1.69398907e-4,
                "axial_force": -29221.3114,
                "current_length": 899.847554,
            },
            17: {
                "axial_stress": -0.3653495,
                "axial_strain": -5.219279e-6,
                "axial_force": -328.81455,
                "current_length": 1499.99217,
                "torsional_stress": 0.0,
            },
            22: {
                "axial_stress": 228.3971,
                "axial_strain": 1.087605e-3,
                "axial_force": 148458.115,
                "current_length": 1205.46911,
                "torsional_stress": 0.0,
            },
            24: {
                "axial_stress": -95.82736,
                "axial_strain": -4.563208e-4,
                "axial_force": -62287.784,
                "current_length": 1203.60997,
                "torsional_stress": 0.0,
            },
        },
        "expected_strain_energy": 4876912.0,
        "expected_increments": 55,
    },
    {
        "case_name": "case_175_warehouse_with_offsets",
        "folder": "case_175",
        "filename": "ej175_mixed_cbar_crod_warehouse_28_elements.dat",
        "description": "Warehouse frame with offsets, three materials and prescribed sway.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            7: {
                0: 5.664712,
                1: 1.179148,
                2: -0.05399930,
                3: -0.001221629,
                4: 0.006530191,
                5: -0.001594645,
            },
            8: {
                0: 10.0,
                1: -2.124265,
                2: -0.7629680,
                3: 0.002032168,
                4: 0.007228253,
                5: -0.003033269,
            },
            10: {
                0: 8.208095,
                1: 0.3616815,
                2: -4.276338,
                3: -0.002434928,
                4: -0.004285634,
                5: -0.001719637,
            },
        },
        "expected_reactions": {
            1: {
                0: -166894.0,
                1: -118324.1,
                2: -26125.97,
                3: 1.308923e7,
                4: -2.679996e7,
                5: 687845.3,
            },
            3: {
                0: -313577.9,
                1: -70331.80,
                2: -93826.48,
                3: 1.146119e7,
                4: -8.547304e7,
                5: 781072.2,
            },
            8: {0: 1.097951e6, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
        },
        "expected_element_types": {
            1: "CBAR",
            7: "CBAR",
            12: "CBAR",
            13: "CBAR",
            14: "CROD",
            17: "CROD",
            21: "CROD",
            28: "CROD",
        },
        "expected_element_results": {
            1: {
                "axial_stress": -23.04735,
                "axial_strain": -1.0974926e-4,
                "axial_force": -69142.05,
                "current_length": 1199.86831,
            },
            7: {
                "axial_stress": -38.36784,
                "axial_strain": -1.8270405e-4,
                "axial_force": -84409.248,
            },
            12: {
                "axial_stress": 32.5890207,
                "axial_strain": 2.8338281e-4,
                "axial_force": 55401.3353,
                "current_length": 2000.56685,
            },
            13: {
                "axial_stress": 249.17265,
                "axial_strain": 2.1667185e-3,
                "axial_force": 423593.505,
                "current_length": 2004.33814,
                "stress_points_end_a.P1.sigma_total": 331.8495,
                "stress_points_end_a.P2.sigma_total": 316.6895,
                "stress_points_end_b.P0.sigma_total": 310.8384,
                "stress_points_end_b.P3.sigma_total": 343.1747,
            },
            14: {
                "axial_stress": 139.1984,
                "axial_strain": 6.628497e-4,
                "axial_force": 104398.8,
                "current_length": 2333.92678,
                "torsional_stress": 0.0,
            },
            17: {
                "axial_stress": -439.6645,
                "axial_strain": -2.093641e-3,
                "axial_force": -329748.375,
                "current_length": 2327.49759,
                "torsional_stress": 0.0,
            },
            21: {
                "axial_stress": 45.3790,
                "axial_strain": 3.946000e-4,
                "axial_force": 40841.10,
                "current_length": 3232.37388,
                "torsional_stress": 0.0,
            },
            28: {
                "axial_stress": 101.6090,
                "axial_strain": 8.835563e-4,
                "axial_force": 91448.10,
                "current_length": 3803.35751,
                "torsional_stress": 0.0,
            },
        },
        "expected_strain_energy": 5893006.0,
        "expected_increments": 60,
    },
    {
        "case_name": "case_176_braced_tower_with_settlement",
        "folder": "case_176",
        "filename": "ej176_mixed_cbar_crod_braced_tower_30_elements.dat",
        "description": "Thirty-element braced tower with several loads and support settlement.",
        "rel_tol": 3.0e-2,
        "abs_tol": 1.0e-5,
        "expected_displacements": {
            4: {0: 0.0, 1: 0.0, 2: -4.0, 3: 0.0, 4: 0.0, 5: 0.0},
            9: {
                0: 1.227783,
                1: 2.421820,
                2: -1.246257,
                3: -0.001919806,
                4: 0.0008868527,
                5: 0.004023843,
            },
            10: {
                0: 1.261553,
                1: 4.417283,
                2: -0.9836593,
                3: -0.003038120,
                4: 0.001840893,
                5: 0.007258938,
            },
            12: {
                0: -1.526421,
                1: 2.213418,
                2: -2.743663,
                3: -0.001739466,
                4: -0.001364104,
                5: 0.006345933,
            },
        },
        "expected_reactions": {
            1: {
                0: -45263.52,
                1: 8529.238,
                2: 603600.5,
                3: 1.361331e8,
                4: -2.780826e7,
                5: -712105.6,
            },
            3: {
                0: -192449.6,
                1: -228995.4,
                2: 1.256365e6,
                3: 5.765227e7,
                4: 1.267082e8,
                5: 6.856129e6,
            },
            4: {
                0: -56330.53,
                1: 199333.3,
                2: -1.177743e6,
                3: 1.325127e8,
                4: 1.338707e8,
                5: -1.327627e6,
            },
        },
        "expected_element_types": {
            1: "CBAR",
            9: "CBAR",
            18: "CBAR",
            19: "CBAR",
            20: "CBAR",
            21: "CROD",
            24: "CROD",
            30: "CROD",
        },
        "expected_element_results": {
            1: {
                "axial_stress": 0.0,
                "axial_strain": 0.0,
                "axial_force": 0.0,
                "current_length": 800.0,
            },
            9: {
                "axial_stress": 10.0287925,
                "axial_strain": 8.720690e-5,
                "axial_force": 15043.1887,
                "current_length": 440.038373,
            },
            18: {
                "axial_stress": -106.573005,
                "axial_strain": -5.074906e-4,
                "axial_force": -277089.813,
                "current_length": 713.780514,
            },
            19: {
                "axial_stress": -197.333003,
                "axial_strain": -9.39680925e-4,
                "axial_force": -513065.807,
                "current_length": 713.472091,
                "stress_points_end_a.P0.sigma_total": -321.6536,
                "stress_points_end_a.P1.sigma_total": -283.0127,
                "stress_points_end_b.P2.sigma_total": -322.9977,
                "stress_points_end_b.P3.sigma_total": -297.6326,
            },
            20: {
                "axial_stress": 45.461745,
                "axial_strain": 2.1648455e-4,
                "axial_force": 118200.537,
                "current_length": 714.297460,
            },
            21: {
                "axial_stress": 51.30231,
                "axial_strain": 2.442967e-4,
                "axial_force": 35911.617,
                "current_length": 1007.61892,
                "torsional_stress": 0.0,
            },
            24: {
                "axial_stress": 326.8244,
                "axial_strain": 1.556306e-3,
                "axial_force": 228777.08,
                "current_length": 1008.94060,
                "torsional_stress": 0.0,
            },
            30: {
                "axial_stress": 78.64115,
                "axial_strain": 6.838361e-4,
                "axial_force": 66844.9775,
                "current_length": 622.679487,
                "torsional_stress": 0.0,
            },
        },
        "expected_strain_energy": 3618518.0,
        "expected_increments": 60,
    },
]


def _nested_result_value(result, path):
    value = result
    for key in path.split("."):
        value = value[key]
    return value


@pytest.mark.parametrize(
    "reference",
    [
        pytest.param(reference, id=reference["case_name"])
        for reference in NONLINEAR_MIXED_VALIDATION_DATA
    ],
)
def test_nonlinear_mixed_crod_cbar_final_results(reference):
    deck = VALIDATION_ROOT / reference["folder"] / reference["filename"]
    assert deck.is_file(), f"Missing DAT for {reference['case_name']}: {deck}"

    solver = NonlinearStaticSolver(parse_nastran_input(str(deck)))
    result = solver.solve()
    rel_tol = reference["rel_tol"]
    abs_tol = reference["abs_tol"]

    assert result["load_factor"] == pytest.approx(1.0)
    assert len(result["increments"]) == reference["expected_increments"]

    for node_id, dofs in reference.get("expected_displacements", {}).items():
        computed = result["displacements"][solver.dof_mapper.get_dof_indices(node_id)]
        for dof, expected in dofs.items():
            assert computed[dof] == pytest.approx(
                expected, rel=rel_tol, abs=abs_tol
            ), f"GRID {node_id} displacement DOF {dof + 1}"

    for node_id, dofs in reference.get("expected_reactions", {}).items():
        computed = result["reactions"][node_id]
        for dof, expected in dofs.items():
            assert computed[dof] == pytest.approx(
                expected, rel=rel_tol, abs=abs_tol
            ), f"GRID {node_id} reaction DOF {dof + 1}"

    for element_id, expected_type in reference["expected_element_types"].items():
        assert result["stresses"][element_id]["type"] == expected_type

    for element_id, expected_results in reference["expected_element_results"].items():
        computed = result["stresses"][element_id]
        for result_name, expected in expected_results.items():
            assert _nested_result_value(computed, result_name) == pytest.approx(
                expected, rel=rel_tol, abs=abs_tol
            ), f"Element {element_id} result {result_name}"

    if "expected_strain_energy" in reference:
        assert result["strain_energy"] == pytest.approx(
            reference["expected_strain_energy"], rel=rel_tol, abs=abs_tol
        )
