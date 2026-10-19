"""Structural checks for the final mixed CROD/CBAR nonlinear bridge model.

The MSC numerical checkpoints belong to the shared mixed E2E module; this file
keeps the deck-integrity, equilibrium and nonlinear-signal checks together.
"""

from pathlib import Path

import numpy as np
import pytest

from src.preprocessing.data_model import (
    AnalysisType,
    BarElement,
    RodElement,
    SinglePointLoad,
    SinglePointMoment,
)
from src.preprocessing.nastran_reader import parse_nastran_input
from src.solver.nonlinear_analysis import NonlinearStaticSolver


CASE_FILE = (
    Path(__file__).parent
    / "validation_cases"
    / "case_nonlinear_test_final_updated"
    / "tfm_final_nonlinear_bridge_updated.dat"
)


@pytest.fixture(scope="module")
def solved_bridge():
    model = parse_nastran_input(str(CASE_FILE))
    solver = NonlinearStaticSolver(model)
    return model, solver, solver.solve()


def test_bridge_contains_the_intended_supported_physics(solved_bridge):
    model, _solver, _results = solved_bridge
    type_counts = {BarElement: 0, RodElement: 0}
    for element in model.elements.values():
        type_counts[type(element)] += 1

    assert len(model.nodes) == 202
    assert type_counts == {BarElement: 144, RodElement: 428}
    assert sum(isinstance(load, SinglePointLoad) for load in model.loads) == 44
    assert sum(isinstance(load, SinglePointMoment) for load in model.loads) == 3
    assert model.analysis.solution_sequence == 400
    assert model.analysis.analysis_type is AnalysisType.NONLINEAR_STATIC
    assert model.analysis.large_displacements
    assert not model.analysis.follower_force_stiffness
    assert model.analysis.nonlinear_step.increment_method == "FIXED"
    assert model.analysis.nonlinear_step.number_of_increments == 20
    assert "PARAM,AUTOSPC,NO" in CASE_FILE.read_text(
        encoding="ascii"
    ).replace(" ", "")

    expected_settlements = {4011: -2.0, 4022: -6.0, 4034: -10.0}
    for grid_id, expected in expected_settlements.items():
        assert any(
            constraint.node_id == grid_id
            and str(constraint.components) == "3"
            and constraint.value == pytest.approx(expected)
            for constraint in model.constraints
        )

    inner_grid_ids = set(range(1101, 1132)) | set(range(1201, 1232))
    constrained_inner_rotations = {
        constraint.node_id
        for constraint in model.constraints
        if str(constraint.components) == "456"
    }
    assert constrained_inner_rotations == inner_grid_ids


def test_bridge_converges_at_full_load_and_satisfies_equilibrium(solved_bridge):
    _model, solver, results = solved_bridge
    displacements = results["displacements"]

    assert results["load_factor"] == pytest.approx(1.0, abs=1.0e-14)
    assert len(results["increments"]) == 20
    assert results["cutbacks"] == 0
    assert np.all(np.isfinite(displacements))

    translations = np.array(
        [
            displacements[solver.dof_mapper.get_dof_indices(node_id)[:3]]
            for node_id in sorted(solver.model.nodes)
        ]
    )
    assert np.max(np.linalg.norm(translations, axis=1)) < 20.0
    for grid_id, expected in ((4011, -2.0), (4022, -6.0), (4034, -10.0)):
        assert displacements[solver.dof_mapper.get_dof_indices(grid_id)[2]] == (
            pytest.approx(expected, abs=1.0e-12)
        )

    final_state = solver.assemble_state(displacements)
    free_dofs, _fixed_dofs, _prescribed = solver.get_dof_partition()
    residual = solver.external_load - final_state.internal_force
    reference = max(1.0, np.max(np.abs(solver.external_load[free_dofs])))
    assert np.max(np.abs(residual[free_dofs])) / reference < 3.0e-7


def test_bridge_response_is_materially_geometrically_nonlinear(solved_bridge):
    """Compare the thesis bridge with the same stage evaluated at LGDISP=-1."""
    _model, nonlinear_solver, nonlinear_results = solved_bridge
    linear_control = parse_nastran_input(str(CASE_FILE))
    linear_control.analysis.large_displacement_parameter = -1
    linear_control.analysis.large_displacements = False
    linear_solver = NonlinearStaticSolver(linear_control)
    linear_results = linear_solver.solve()

    nonlinear_t2 = nonlinear_results["displacements"][
        nonlinear_solver.dof_mapper.get_dof_indices(1015)[1]
    ]
    linear_t2 = linear_results["displacements"][
        linear_solver.dof_mapper.get_dof_indices(1015)[1]
    ]
    relative_change = abs(nonlinear_t2 - linear_t2) / abs(linear_t2)
    assert relative_change > 0.015
