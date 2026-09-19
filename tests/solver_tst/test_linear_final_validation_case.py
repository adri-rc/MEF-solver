"""Structural and MSC Nastran checks for the final linear defence model.

The external values are transcribed from ``tfm_final_airframe_frame.f06``,
generated with MSC Nastran Student Edition 2026.1 on 2026-08-21.  The test does
not read F06/OP2 at runtime, so every reference remains explicit and auditable.
"""

from pathlib import Path

import numpy as np
import pytest

from src.preprocessing.data_model import (
    BarElement,
    CQUAD4Element,
    CTRIA3Element,
    RodElement,
    SinglePointLoad,
    SinglePointMoment,
)
from src.preprocessing.nastran_reader import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver


CASE_FILE = (
    Path(__file__).parent
    / "validation_cases"
    / "case_linear_test_final"
    / "tfm_final_airframe_frame.dat"
)


MSC_DISPLACEMENTS = {
    1: (2.614711e-2, 1.071225e-1, 0.0, 0.0, 0.0, -1.687735e-4),
    165: (5.151528e-2, 1.213636e-1, 0.0, 0.0, 0.0, 2.438577e-3),
    376: (3.176561e-2, 7.440119e-1, 0.0, 0.0, 0.0, 0.0),
    377: (5.569107e-2, 8.000000e-1, 0.0, 0.0, 0.0, 0.0),
}

MSC_SPC_REACTIONS = {
    181: (-1.209305e3, -6.517070e3, 0.0, 0.0, 0.0, -6.553100e3),
    184: (-1.397134e3, -1.674250e3, 0.0, 0.0, 0.0, -3.024559e3),
    187: (-2.937250e3, -5.394952e3, 0.0, 0.0, 0.0, -1.873448e4),
    366: (1.797159e3, -4.596829e3, 0.0, 0.0, 0.0, 2.002545e4),
    369: (3.720260e2, -1.609596e3, 0.0, 0.0, 0.0, -6.530562e2),
    372: (8.745043e2, -5.201126e3, 0.0, 0.0, 0.0, 7.713286e3),
    377: (0.0, 1.299382e4, 0.0, 0.0, 0.0, 0.0),
}

# axial force, torque, axial stress and torsional stress
MSC_CROD_RESULTS = {
    8001: (7.063809e3, 0.0, 1.412762e2, 0.0),
    8003: (7.648824e3, 0.0, 1.529765e2, 0.0),
}

# CBAR force rows follow the F06 order:
# MA1, MA2, MB1, MB2, shear 1, shear 2, axial force and torque.
MSC_CBAR_FORCES = {
    7001: (0.0, -1.150988e4, 0.0, 8.875711e3, 0.0, -2.047848e2, 2.314552e2, 0.0),
    7025: (0.0, 8.117794e4, 0.0, -3.406619e4, 0.0, 1.157693e3, -7.562862e2, 0.0),
    7056: (0.0, 2.011640e4, 0.0, -8.650859e4, 0.0, 1.071108e3, 3.153551e2, 0.0),
}

# Four bending fibres and the axial/max/min stresses at end A, followed by the
# four bending fibres and max/min stresses at end B.
MSC_CBAR_STRESSES = {
    7001: (
        -5.137369, -5.137369, 17.33378, 17.33378, 0.9849159, 18.31869,
        -4.152453, 3.961624, 3.961624, -13.36675, -13.36675, 4.946540,
        -12.38183,
    ),
    7025: (
        36.23332, 36.23332, -122.2533, -122.2533, -3.218239, 33.01508,
        -125.4715, -15.20526, -15.20526, 51.30340, 51.30340, 48.08517,
        -18.42349,
    ),
    7056: (
        8.978845, 8.978845, -30.29514, -30.29514, 1.341937, 10.32078,
        -28.95320, -38.61263, -38.61263, 130.2812, 130.2812, 131.6232,
        -37.27069,
    ),
}

# Membrane resultants Nx, Ny, Nxy and top-fibre sx, sy, txy, von Mises.
MSC_SHELL_RESULTS = {
    1001: (10.31059, 36.98256, 5.675639, 5.155296, 18.49128, 2.837819, 17.24344),
    1102: (41.60005, 20.26746, 38.83020, 20.80002, 10.13373, 19.41510, 38.14958),
    1256: (1.578388, 14.14789, -9.883197, 0.7891940, 7.073944, -4.941598, 10.87837),
    5001: (1.275846, 40.61016, 4.265709, 0.6379229, 20.30508, 2.132855, 20.33217),
    5030: (47.98331, 54.58718, 39.82791, 23.99166, 27.29359, 19.91396, 43.07457),
    5128: (29.98747, 51.72194, 29.83225, 14.99373, 25.86097, 14.91612, 34.25378),
}


@pytest.fixture(scope="module")
def solved_case():
    model = parse_nastran_input(str(CASE_FILE))
    solver = LinearStaticSolver(model)
    solver.assemble()
    return model, solver, solver.solve()


def test_final_case_contains_the_intended_supported_physics(solved_case):
    model, _solver, _results = solved_case
    type_counts = {
        RodElement: 0,
        BarElement: 0,
        CTRIA3Element: 0,
        CQUAD4Element: 0,
    }
    for element in model.elements.values():
        type_counts[type(element)] += 1

    assert len(model.nodes) == 377
    assert type_counts == {
        RodElement: 4,
        BarElement: 76,
        CTRIA3Element: 128,
        CQUAD4Element: 256,
    }
    assert sum(isinstance(load, SinglePointLoad) for load in model.loads) == 2
    assert sum(isinstance(load, SinglePointMoment) for load in model.loads) == 1
    assert any(
        constraint.node_id == 377
        and str(constraint.components) == "2"
        and constraint.value == pytest.approx(0.8)
        for constraint in model.constraints
    )


def test_final_case_solves_with_prescribed_motion_and_global_equilibrium(solved_case):
    model, solver, results = solved_case
    displacements = results["displacements"]
    assert np.all(np.isfinite(displacements))
    right_clevis_t2 = solver.dof_mapper.get_dof_indices(377)[1]
    assert displacements[right_clevis_t2] == pytest.approx(0.8, abs=1.0e-12)
    assert np.max(np.abs(displacements)) < 5.0

    applied_force = np.zeros(3)
    applied_moment = np.zeros(3)
    for load in model.loads:
        vector = load.module * np.asarray(load.vector, dtype=float)
        position = model.nodes[load.appnode_id].get_coords()
        if isinstance(load, SinglePointLoad):
            applied_force += vector
            applied_moment += np.cross(position, vector)
        else:
            applied_moment += vector

    reaction_force = np.zeros(3)
    reaction_moment = np.zeros(3)
    for node_id, reaction in results["reactions"].items():
        position = model.nodes[node_id].get_coords()
        reaction_force += reaction[:3]
        reaction_moment += reaction[3:] + np.cross(position, reaction[:3])

    assert reaction_force + applied_force == pytest.approx(
        np.zeros(3), rel=1.0e-9, abs=1.0e-6
    )
    assert reaction_moment + applied_moment == pytest.approx(
        np.zeros(3), rel=1.0e-9, abs=1.0e-3
    )


def test_final_case_matches_transcribed_msc_nastran_reference(solved_case):
    """Correlate all supported element families with explicit F06 values."""
    _model, solver, results = solved_case
    displacements = results["displacements"]

    for node_id, expected in MSC_DISPLACEMENTS.items():
        dofs = solver.dof_mapper.get_dof_indices(node_id)
        assert displacements[dofs] == pytest.approx(expected, rel=2.0e-5, abs=2.0e-8)

    for node_id, expected in MSC_SPC_REACTIONS.items():
        assert results["reactions"][node_id] == pytest.approx(
            expected, rel=1.0e-4, abs=5.0e-4
        )

    for element_id, expected in MSC_CROD_RESULTS.items():
        data = results["stresses"][element_id]
        actual = (
            data["axial_force"],
            data["torque"],
            data["axial_stress"],
            data["torsional_stress"],
        )
        assert actual == pytest.approx(expected, rel=1.0e-6, abs=1.0e-5)

    for element_id, expected in MSC_CBAR_FORCES.items():
        data = results["stresses"][element_id]
        end_a = data["forces_n1"]
        end_b = data["forces_n2"]
        actual = (
            -end_a["Mz"],
            end_a["My"],
            end_b["Mz"],
            -end_b["My"],
            -end_a["Fy"],
            -end_a["Fz"],
            data["axial_force"],
            data["torque"],
        )
        assert actual == pytest.approx(expected, rel=5.0e-4, abs=5.0e-4)

    for element_id, expected in MSC_CBAR_STRESSES.items():
        data = results["stresses"][element_id]
        end_a = data["stress_points_end_a"]
        end_b = data["stress_points_end_b"]
        actual = tuple(end_a[f"P{i}"]["sigma_bending"] for i in range(4))
        actual += (
            data["axial_stress"],
            max(end_a[f"P{i}"]["sigma_total"] for i in range(4)),
            min(end_a[f"P{i}"]["sigma_total"] for i in range(4)),
        )
        actual += tuple(end_b[f"P{i}"]["sigma_bending"] for i in range(4))
        actual += (
            max(end_b[f"P{i}"]["sigma_total"] for i in range(4)),
            min(end_b[f"P{i}"]["sigma_total"] for i in range(4)),
        )
        assert actual == pytest.approx(expected, rel=5.0e-4, abs=5.0e-4)

    for element_id, expected in MSC_SHELL_RESULTS.items():
        data = results["stresses"][element_id]
        actual = (
            data["membrane_force_x"],
            data["membrane_force_y"],
            data["membrane_force_xy"],
            data["stress_top_fiber_x"],
            data["stress_top_fiber_y"],
            data["tau_top_fiber_xy"],
            data["von_mises_top"],
        )
        assert actual == pytest.approx(expected, rel=1.0e-3, abs=5.0e-4)
