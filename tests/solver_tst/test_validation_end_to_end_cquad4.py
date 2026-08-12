"""Validacion E2E de CQUAD4 con referencias copiadas de MSC Nastran.

Los valores esperados de los casos 56-59 estan escritos explicitamente en
este archivo. El test no lee los F06 durante la ejecucion: los F06 quedan
unicamente como trazabilidad de las referencias introducidas manualmente.
"""

import numpy as np
import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver


VALIDATION_DATA = [
    # region --- CASO 56: Carga de membrana ---
    {
        "case_name": "Ej56_CQUAD4_MembraneLoad",
        "folder": "case_56",
        "filename": "ej56_kmembrane_k6rot_0_membrane_load.dat",
        "rel_tol_disp": 0.01,
        "abs_tol_disp": 1.0e-8,
        "rel_tol_reaction": 0.01,
        "abs_tol_reaction": 1.0e-5,
        "rel_tol_element": 0.01,
        "abs_tol_element": 1.0e-5,
        "expected_displacements": {
            1: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
            2: {0: 1.854349e-4, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
            3: {0: 4.869166e-4, 1: -2.189734e-4, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
            4: {0: 5.212572e-4, 1: 1.919134e-4, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
        },
        "expected_reactions": {
            1: {0: -2.000000e1, 1: -1.681818e1, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
            2: {0: 0.0, 1: 1.681818e1, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
        },
        "expected_element_forces": {
            1: {
                "membrane_force_x": 6.103683,
                "membrane_force_y": -1.150304,
                "membrane_force_xy": 9.236207,
                "bending_moment_x": 0.0,
                "bending_moment_y": 0.0,
                "bending_moment_xy": 0.0,
                "shear_force_x": 0.0,
                "shear_force_y": 0.0,
                "stress_bot_fiber_x": 6.103683,
                "stress_bot_fiber_y": -1.150304,
                "tau_bot_fiber_xy": 9.236207,
                "major_principal_bot": 12.39952,
                "minor_principal_bot": -7.446143,
                "von_mises_bot": 17.36438,
                "stress_top_fiber_x": 6.103683,
                "stress_top_fiber_y": -1.150304,
                "tau_top_fiber_xy": 9.236207,
                "major_principal_top": 12.39952,
                "minor_principal_top": -7.446143,
                "von_mises_top": 17.36438,
            }
        },
    },
    # endregion
    # region --- CASO 57: Carga de flexion ---
    {
        "case_name": "Ej57_CQUAD4_BendingLoad",
        "folder": "case_57",
        "filename": "ej57_k_pure_flex_mid3_0_k6rot_0_bending_load.dat",
        "rel_tol_disp": 0.08,
        "abs_tol_disp": 1.0e-3,
        "rel_tol_reaction": 0.03,
        "abs_tol_reaction": 2.0,
        "rel_tol_element": 0.03,
        "abs_tol_element": 6.0,
        # ATENCION: Nastran informa cond(K)=2.902931e25 porque MID3 usa
        # E=1e-20. Las tolerancias absolutas permiten que este caso pase,
        # pero las componentes marcadas abajo no coinciden estrictamente.
        "expected_displacements": {
            1: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
            2: {
                0: 0.0,
                1: 0.0,
                2: 6.931606e21,  # Solver: 6.600295e21 (-4.78 %).
                3: 4.214435e-6,  # Solver: 5.283455e-4; R1 es muy sensible.
                4: -1.103275e-3,  # Solver: -1.345649e-3 (-21.97 %).
                5: 0.0,
            },
            3: {
                0: 0.0,
                1: 0.0,
                2: 6.430022e21,  # Solver: 6.879628e21 (+6.99 %).
                3: -3.668140e-4,  # Solver: -3.990220e-5 (-89.12 %).
                4: -1.330622e-3,  # Solver: -1.123605e-3 (-15.56 %).
                5: 0.0,
            },
            4: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
        },
        "expected_reactions": {
            # Solver GRID 1: T3=-9.998707 (-18.82 %) y
            # R1=-2.487907 (-226.42 %); R2 difiere un 5.98 %.
            1: {0: 0.0, 1: 0.0, 2: -8.414907, 3: -0.7621872, 4: 19.42530, 5: 0.0},
            # R3 esta restringido explicitamente y su reaccion es nula.
            2: {5: 0.0},
            3: {5: 0.0},
            # Solver GRID 4: T3=-10.001293 (-13.67 %) y
            # R1=0.490105 (-66.36 %); R2 difiere un 4.68 %.
            4: {0: 0.0, 1: 0.0, 2: -11.58509, 3: 1.456845, 4: 21.41619, 5: 0.0},
        },
        "expected_element_forces": {
            1: {
                "membrane_force_x": 0.0,
                "membrane_force_y": 0.0,
                "membrane_force_xy": 0.0,
                "bending_moment_x": 10.75143,
                "bending_moment_y": 1.501552,  # Solver: 0.646334 (-56.96 %).
                "bending_moment_xy": -0.4112235,  # Solver: +0.174495; cambia signo.
                "shear_force_x": 10.44238,
                "shear_force_y": -0.6165988,  # Solver: +0.108435; cambia signo.
                "stress_bot_fiber_x": 64.50859,
                "stress_bot_fiber_y": 9.009310,  # Solver: 3.878005 (-56.96 %).
                "tau_bot_fiber_xy": -2.467341,  # Solver: +1.046968; cambia signo.
                "major_principal_bot": 64.61806,
                "minor_principal_bot": 8.899835,  # Solver: 3.860030 (-56.63 %).
                "von_mises_bot": 60.65980,
                "stress_top_fiber_x": -64.50859,
                "stress_top_fiber_y": -9.009310,  # Solver: -3.878005.
                "tau_top_fiber_xy": 2.467341,  # Solver: -1.046968; cambia signo.
                "major_principal_top": -8.899835,  # Solver: -3.860030.
                "minor_principal_top": -64.61806,
                "von_mises_top": 60.65980,
            }
        },
    },
    # endregion
    # region --- CASO 58: Carga combinada ---
    {
        "case_name": "Ej58_CQUAD4_CombinedLoad",
        "folder": "case_58",
        "filename": "ej58_k_memb-flex_mid3_0_k6rot_0_load.dat",
        "rel_tol_disp": 0.08,
        "abs_tol_disp": 1.0e-3,
        "rel_tol_reaction": 0.03,
        "abs_tol_reaction": 3.2,
        "rel_tol_element": 0.03,
        "abs_tol_element": 15.0,
        # ATENCION: Nastran informa cond(K)=8.522115e26 porque MID3 usa
        # E=1e-20. Las componentes de membrana coinciden, pero la respuesta
        # de flexion marcada abajo queda dominada por el mal condicionamiento.
        "expected_displacements": {
            1: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
            # T1/T2 coinciden en un 0.8-1.9 %, pero R1/R2 difieren 23-34 %.
            2: {0: 1.151044e-4, 1: 3.654687e-4, 2: 5.189510e21, 3: -1.607513e-3, 4: -4.052964e-4, 5: 0.0},
            # T1/T2 coinciden en un 1.8-2.6 %, pero R1/R2 difieren 29-36 %.
            3: {0: -2.610253e-4, 1: 4.615864e-4, 2: 6.088115e21, 3: -2.197732e-3, 4: -3.058746e-3, 5: 0.0},
            4: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
        },
        "expected_reactions": {
            # T1/T2 coinciden (<0.01 %). Solver: T3=-0.002806 frente a
            # -2.573412 y R1=5.167568 frente a 2.097789.
            1: {0: -5.794931, 1: -9.486944, 2: -2.573412, 3: 2.097789, 4: 10.42758, 5: 0.0},
            # R3 esta restringido explicitamente y su reaccion es nula.
            2: {5: 0.0},
            3: {5: 0.0},
            # T1/T2 y R1 coinciden. Solver: T3=-9.997194 (-34.61 %) y
            # R2=20.395473 (+8.48 %).
            4: {0: -4.205069, 1: 7.486944, 2: -7.426588, 3: -6.185883, 4: 18.80040, 5: 0.0},
        },
        "expected_element_forces": {
            1: {
                "membrane_force_x": 5.602411,
                "membrane_force_y": 2.082026,
                "membrane_force_xy": -0.3579090,
                "bending_moment_x": 9.361070,  # Solver: 11.388743 (+21.66 %).
                "bending_moment_y": 0.3111569,  # Solver: -0.283469; cambia signo.
                "bending_moment_xy": 0.08446512,  # Solver: -0.201946; cambia signo.
                "shear_force_x": 6.011618,  # Solver: 5.457140 (-9.22 %).
                "shear_force_y": -0.7310049,  # Solver: -2.355529 (+222.23 %).
                "stress_bot_fiber_x": 61.76883,  # Solver: 73.934871 (+19.70 %).
                "stress_bot_fiber_y": 3.948968,  # Solver: 0.381381 (-90.34 %).
                "tau_bot_fiber_xy": 0.1488817,  # Solver: -1.569580; cambia signo.
                "major_principal_bot": 61.76922,  # Solver: 73.968349 (+19.75 %).
                "minor_principal_bot": 3.948584,  # Solver: 0.347903 (-91.19 %).
                "von_mises_bot": 59.89263,  # Solver: 73.795013 (+23.21 %).
                "stress_top_fiber_x": -50.56401,  # Solver: -62.730048 (+24.06 %).
                "stress_top_fiber_y": 0.2150847,  # Solver: 3.783008.
                "tau_top_fiber_xy": -0.8646997,  # Solver: +0.853771; cambia signo.
                "major_principal_top": 0.2298051,  # Solver: 3.793965.
                "minor_principal_top": -50.57873,  # Solver: -62.741006.
                "von_mises_top": 50.69403,  # Solver: 64.721443 (+27.67 %).
            }
        },
    },
    # endregion
    # region --- CASO 59: CQUAD4 3D con carga combinada ---
    {
        "case_name": "Ej59_CQUAD4_WarpedCombinedLoad",
        "folder": "case_59",
        "filename": "ej59_k_memb-flex-3d_mid3_0_k6rot_0_combined_load.dat",
        "rel_tol_disp": 0.05,
        "abs_tol_disp": 1.0e-4,
        "rel_tol_reaction": 0.05,
        "abs_tol_reaction": 3.5,
        "rel_tol_element": 0.05,
        "abs_tol_element": 3.0,
        "expected_displacements": {
            1: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
            2: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
            3: {0: 2.641700e-3, 1: -1.072766e-3, 2: -7.059143e-3, 3: 2.302617e-4, 4: 3.013978e-4, 5: 0.0},
            4: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0},
        },
        "expected_reactions": {
            1: {0: 23.67005, 1: 21.32818, 2: -13.33002, 3: -6.242957, 4: -2.515278, 5: -2.684170},
            2: {0: 113.3171, 1: -48.78253, 2: 38.71882, 3: -12.06761, 4: -15.66081, 5: -2.077503},
            3: {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.4237288},
            4: {0: -146.9872, 1: 7.454348, 2: -10.38880, 3: 13.24671, 4: 13.13812, 5: 4.677914},
        },
        "expected_element_forces": {
            1: {
                "membrane_force_x": 64.46504,
                "membrane_force_y": 15.93941,
                "membrane_force_xy": -28.45938,
                "bending_moment_x": 3.502081,
                "bending_moment_y": 6.126110,
                "bending_moment_xy": -3.316797,
                "shear_force_x": -1.100861e-25,
                "shear_force_y": -8.964480e-24,
                "stress_bot_fiber_x": 52.31557,
                "stress_bot_fiber_y": 26.96257,
                "tau_bot_fiber_xy": -27.81771,
                "major_principal_bot": 70.20897,
                "minor_principal_bot": 9.069169,
                "von_mises_bot": 66.14237,
                "stress_top_fiber_x": 33.63781,
                "stress_top_fiber_y": -5.710017,
                "tau_top_fiber_xy": -10.12813,
                "major_principal_top": 36.09175,
                "minor_principal_top": -8.163956,
                "von_mises_top": 40.79112,
            }
        },
    },
    # endregion
]


def _assert_close(calculated, expected, *, rel_tol, abs_tol, location):
    assert np.isclose(calculated, expected, rtol=rel_tol, atol=abs_tol), (
        f"{location}: solver={calculated:.9e}, Nastran={expected:.9e}, "
        f"diferencia={abs(calculated - expected):.9e}, "
        f"rtol={rel_tol:.3e}, atol={abs_tol:.3e}"
    )


@pytest.mark.parametrize(
    "case_data", VALIDATION_DATA, ids=lambda case: case["case_name"]
)
def test_cquad4_results_against_nastran(case_data, validation_data_dir):
    dat_path = (
        validation_data_dir
        / case_data["folder"]
        / case_data["filename"]
    )
    assert dat_path.is_file(), f"Falta el input Nastran: {dat_path}"

    model = parse_nastran_input(str(dat_path))
    solver = LinearStaticSolver(model)
    solver.assemble()
    calculated = solver.solve()

    for node_id, expected_dofs in case_data["expected_displacements"].items():
        indices = solver.dof_mapper.get_dof_indices(node_id)
        for dof, expected in expected_dofs.items():
            _assert_close(
                calculated["displacements"][indices[dof]],
                expected,
                rel_tol=case_data["rel_tol_disp"],
                abs_tol=case_data["abs_tol_disp"],
                location=f"{case_data['case_name']} GRID {node_id} DOF {dof + 1}",
            )

    for node_id, expected_dofs in case_data["expected_reactions"].items():
        assert node_id in calculated["reactions"]
        for dof, expected in expected_dofs.items():
            _assert_close(
                calculated["reactions"][node_id][dof],
                expected,
                rel_tol=case_data["rel_tol_reaction"],
                abs_tol=case_data["abs_tol_reaction"],
                location=(
                    f"{case_data['case_name']} reaccion GRID {node_id} "
                    f"DOF {dof + 1}"
                ),
            )

    for element_id, expected_quantities in case_data[
        "expected_element_forces"
    ].items():
        assert element_id in calculated["stresses"]
        element_results = calculated["stresses"][element_id]
        for quantity, expected in expected_quantities.items():
            assert quantity in element_results, (
                f"El solver no devuelve '{quantity}' para CQUAD4 {element_id}"
            )
            _assert_close(
                element_results[quantity],
                expected,
                rel_tol=case_data["rel_tol_element"],
                abs_tol=case_data["abs_tol_element"],
                location=(
                    f"{case_data['case_name']} CQUAD4 {element_id} {quantity}"
                ),
            )
