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
from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.nonlinear_analysis import NonlinearStaticSolver


CASE_FILE = (
    Path(__file__).parent
    / "validation_cases"
    / "case_nonlinear_test_final"
    / "tfm_final_nonlinear_bridge.dat"
)


@pytest.fixture(scope="module")
def solved_bridge():
    model = parse_nastran_input(str(CASE_FILE))
    solver = NonlinearStaticSolver(model)
    return model, solver, solver.solve()


def test_final_bridge_contains_the_intended_supported_physics(solved_bridge):
    model, _solver, _results = solved_bridge
    type_counts = {BarElement: 0, RodElement: 0}
    for element in model.elements.values():
        type_counts[type(element)] += 1

    assert len(model.nodes) == 32
    assert type_counts == {BarElement: 48, RodElement: 30}
    assert sum(isinstance(load, SinglePointLoad) for load in model.loads) == 15
    assert sum(isinstance(load, SinglePointMoment) for load in model.loads) == 3
    assert model.analysis.solution_sequence == 400
    assert model.analysis.analysis_type is AnalysisType.NONLINEAR_STATIC
    assert model.analysis.large_displacements
    assert not model.analysis.follower_force_stiffness
    assert model.analysis.nonlinear_step.increment_method == "FIXED"
    assert model.analysis.nonlinear_step.number_of_increments == 50
    assert any(
        constraint.node_id == 109
        and str(constraint.components) == "3"
        and constraint.value == pytest.approx(-12.0)
        for constraint in model.constraints
    )
    assert any(
        constraint.node_id == 209
        and str(constraint.components) == "3"
        and constraint.value == pytest.approx(-6.0)
        for constraint in model.constraints
    )


def test_final_bridge_converges_at_full_load_and_satisfies_equilibrium(solved_bridge):
    _model, solver, results = solved_bridge
    displacements = results["displacements"]

    assert results["load_factor"] == pytest.approx(1.0, abs=1.0e-14)
    assert len(results["increments"]) == 50
    assert results["cutbacks"] == 0
    assert np.all(np.isfinite(displacements))
    assert np.max(np.abs(displacements)) < 50.0
    assert displacements[solver.dof_mapper.get_dof_indices(109)[2]] == pytest.approx(
        -12.0, abs=1.0e-12
    )
    assert displacements[solver.dof_mapper.get_dof_indices(209)[2]] == pytest.approx(
        -6.0, abs=1.0e-12
    )

    final_state = solver.assemble_state(displacements)
    free_dofs, _fixed_dofs, _prescribed = solver.get_dof_partition()
    residual = solver.external_load - final_state.internal_force
    reference = max(1.0, np.max(np.abs(solver.external_load[free_dofs])))
    assert np.max(np.abs(residual[free_dofs])) / reference < 3.0e-7


def test_final_bridge_response_is_materially_geometrically_nonlinear(solved_bridge):
    """Compare with the same SOL 400 stage evaluated with LGDISP=-1."""
    _model, nonlinear_solver, nonlinear_results = solved_bridge
    linear_control = parse_nastran_input(str(CASE_FILE))
    linear_control.analysis.large_displacement_parameter = -1
    linear_control.analysis.large_displacements = False
    linear_solver = NonlinearStaticSolver(linear_control)
    linear_results = linear_solver.solve()

    nonlinear_t1 = nonlinear_results["displacements"][
        nonlinear_solver.dof_mapper.get_dof_indices(105)[0]
    ]
    linear_t1 = linear_results["displacements"][
        linear_solver.dof_mapper.get_dof_indices(105)[0]
    ]
    relative_change = abs(nonlinear_t1 - linear_t1) / abs(linear_t1)
    assert relative_change > 0.015
