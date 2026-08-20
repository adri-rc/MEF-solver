"""Final CROD SOL 400 validation against explicit MSC Nastran values.

No F06/PCH file is read by these tests.  Every reference below was
transcribed from the converged 2026.1 outputs in ``validation_cases``.
"""

from pathlib import Path

import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.nonlinear_analysis import NonlinearStaticSolver


VALIDATION_ROOT = Path(__file__).parent / "validation_cases"

# DOF indices: 0=T1, 1=T2, 2=T3, 3=R1, 4=R2, 5=R3.
NONLINEAR_VALIDATION_DATA = [
    {"case": 92, "u": {2: {0: 1.0}}, "r": {1: {0: -2.1e5}}},
    {"case": 93, "u": {2: {0: 1.0}}, "s": {1: {"axial_stress": 2100.0, "axial_strain": 0.01}}},
    {"case": 94, "u": {2: {0: 20.0}}, "r": {1: {0: -4.2e6}}, "s": {1: {"axial_force": 4.2e6, "axial_stress": 42000.0, "axial_strain": 0.2}}, "energy": 4.2e7, "increments": 20},
    {"case": 95, "u": {2: {0: 20.0}}, "s": {1: {"axial_stress": 42000.0}}, "increments": 9},
    {"case": 96, "u": {2: {0: -10.0}}, "r": {1: {0: 2.1e6}}, "s": {1: {"axial_stress": -21000.0, "axial_strain": -0.1}}},
    {"case": 97, "u": {2: {0: 12.24745}}, "r": {1: {0: -2.1e6}}, "s": {1: {"axial_stress": 21000.0, "axial_strain": 0.1}}},
    {"case": 98, "u": {2: {0: 10.0}, 3: {0: 25.0}}, "s": {1: {"axial_stress": 21000.0}, 2: {"axial_stress": 21000.0}}},
    {"case": 99, "u": {2: {0: 2.666667}, 3: {0: 20.66667}}, "r": {1: {0: -8.4e5}}, "s": {1: {"axial_stress": 7000.0}, 2: {"axial_stress": 10500.0}}},
    {"case": 100, "u": {2: {0: 10.0}}, "r": {1: {0: -2.45e6}}, "s": {1: {"axial_stress": 21000.0}, 2: {"axial_stress": 7000.0}}},
    {"case": 101, "u": {2: {0: 8.095238}}, "r": {1: {0: -1.7e6}}, "s": {1: {"axial_stress": 17000.0, "axial_strain": 0.08095238}}},
    {"case": 102, "u": {3: {1: -223.3105}}, "r": {1: {0: -1.621922e6, 1: 2.0e6}, 2: {0: 1.621922e6, 1: 2.0e6}}, "s": {1: {"axial_stress": 25750.01, "axial_strain": 0.1226191}, 2: {"axial_stress": 25750.01}}},
    {"case": 103, "u": {3: {0: 9.454499, 1: -39.27897}}, "r": {1: {0: 1.608725e6, 1: 7.545773e5}, 2: {0: -2.608725e6, 1: 2.245423e6}}, "s": {1: {"axial_stress": -17769.03}, 2: {"axial_stress": -57366.68}}},
    {"case": 104, "u": {3: {1: -6.221705}}, "r": {1: {0: 266264.9, 1: 50000.0}}, "s": {1: {"axial_stress": -2709.188, "axial_strain": -0.01290089}}, "increments": 14},
    {"case": 105, "u": {2: {0: -13.27239}, 3: {0: -6.636195, 1: -233.9893}}, "r": {1: {1: 2.0e6}, 2: {1: 2.0e6}}, "s": {1: {"axial_stress": -13936.01}, 2: {"axial_stress": 32501.96}, 3: {"axial_stress": 32501.96}}},
    {"case": 106, "u": {3: {0: 22.52493, 1: -12.27022}, 4: {0: 18.92084, 1: -4.190009}}, "r": {1: {0: -576573.9}, 2: {0: -423426.1, 1: 2.998984e6}}, "s": {3: {"axial_stress": 8229.271}, 6: {"axial_stress": -23623.47}}},
    {"case": 107, "u": {4: {0: 10.33237, 1: -5.384856, 2: -9.4136}}, "r": {2: {0: -797209.0, 1: 475349.0, 2: 1.453991e6}}, "s": {2: {"axial_stress": -24642.71, "axial_strain": -0.1173462}}},
    {"case": 108, "u": {5: {0: 10.92065, 1: 1.751269, 2: -21.62411}}, "r": {3: {0: -1.318034e6, 1: -1.492986e6, 2: 1.877012e6}}, "s": {3: {"axial_stress": -27366.76, "axial_strain": -0.1303179}}},
    {"case": 109, "u": {3: {1: -6.221705}}, "s": {1: {"axial_stress": -2709.188}}, "increments": 5},
    {"case": 110, "u": {3: {1: -6.221705}}, "s": {1: {"axial_stress": -2709.188}}, "increments": 40},
    {"case": 111, "u": {5: {0: 2.536503, 1: 0.03771908, 2: -5.326263}, 6: {0: -7.533348, 1: 44.10162, 2: -17.23535}}, "r": {3: {0: -870512.1, 1: -459667.5, 2: 1.442809e6}}, "s": {5: {"axial_stress": 120.2893}, 9: {"axial_stress": -22619.5}}, "increments": 15},
    {"case": 112, "u": {3: {1: -15.38383}}, "r": {1: {0: 1.181807e6, 1: 1.0e6}}, "s": {1: {"axial_force": -1.548118e6, "axial_stress": -15481.18, "axial_strain": -0.07371989}}},
    {"case": 113, "u": {3: {0: 0.0, 1: -15.38383, 2: 0.0}}, "r": {1: {0: 585445.7, 1: 1.189084e6, 2: 8.0e5}}, "s": {1: {"axial_stress": -15481.18}}},
    {"case": 114, "u": {3: {0: 0.006127186, 1: -15.38383}}, "s": {1: {"axial_stress": -15474.23}, 2: {"axial_stress": -15488.12}}},
    {"case": 115, "u": {3: {0: -0.006127186, 1: -15.38383}}, "s": {1: {"axial_stress": -15488.12}, 2: {"axial_stress": -15474.23}}},
    {"case": 116, "u": {2: {3: 0.9904762}}, "r": {1: {3: -40000.0}}, "s": {1: {"torque": 40000.0, "torsional_stress": 2000.0}}, "energy": 19809.52},
    {"case": 117, "u": {2: {0: 20.0, 3: 0.9904762}}, "r": {1: {0: -4.2e6, 3: -40000.0}}, "s": {1: {"axial_stress": 42000.0, "torque": 40000.0, "torsional_stress": 2000.0}}, "energy": 4.201981e7},
    {"case": 118, "u": {3: {1: -15.38383}}, "r": {2: {0: -1.181807e6, 1: 1.0e6}}, "s": {2: {"axial_stress": -15481.18}}},
    {"case": 119, "u": {3: {1: -40.32496}}, "r": {1: {0: 3.183911e6, 1: 1.9e6}}, "s": {1: {"axial_stress": -37077.34, "axial_strain": -0.1765587}}},
    {"case": 121, "u": {3: {1: -0.269374}}, "r": {1: {0: 20000.0, 1: 20000.0}}},
    {"case": 122, "u": {3: {1: -0.2699212}}, "r": {1: {0: 20054.13, 1: 20000.0}}, "s": {1: {"axial_stress": -283.2257, "axial_strain": -0.001348694}}},
]


@pytest.mark.parametrize(
    "reference", NONLINEAR_VALIDATION_DATA, ids=lambda item: f"case_{item['case']}"
)
def test_nonlinear_crod_final_results(reference):
    case_id = reference["case"]
    deck = next((VALIDATION_ROOT / f"case_{case_id}").glob("*.dat"))
    model = parse_nastran_input(str(deck))
    solver = NonlinearStaticSolver(model)
    result = solver.solve()

    assert result["load_factor"] == pytest.approx(1.0)
    if "increments" in reference:
        assert len(result["increments"]) == reference["increments"]

    for node_id, dofs in reference.get("u", {}).items():
        nodal = result["displacements"][solver.dof_mapper.get_dof_indices(node_id)]
        for dof, expected in dofs.items():
            assert nodal[dof] == pytest.approx(expected, rel=3e-5, abs=3e-7)

    for node_id, dofs in reference.get("r", {}).items():
        nodal = result["reactions"][node_id]
        for dof, expected in dofs.items():
            assert nodal[dof] == pytest.approx(expected, rel=3e-5, abs=5e-2)

    for element_id, expected_results in reference.get("s", {}).items():
        computed = result["stresses"][element_id]
        for key, expected in expected_results.items():
            assert computed[key] == pytest.approx(expected, rel=3e-5, abs=3e-5)

    if "energy" in reference:
        assert result["strain_energy"] == pytest.approx(
            reference["energy"], rel=3e-7, abs=5e-2
        )
