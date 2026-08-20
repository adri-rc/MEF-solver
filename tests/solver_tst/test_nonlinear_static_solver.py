from pathlib import Path

import numpy as np
import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.nonlinear_analysis import (
    NonlinearConvergenceError,
    NonlinearStaticSolver,
)


VALIDATION_ROOT = Path(__file__).parent / "validation_cases"


def _solver(case_id):
    path = next((VALIDATION_ROOT / f"case_{case_id}").glob("*.dat"))
    return NonlinearStaticSolver(parse_nastran_input(str(path)))


def _node_result(solver, vector, node_id):
    return vector[solver.dof_mapper.get_dof_indices(node_id)]


def test_large_axial_extension_converges_incrementally():
    solver = _solver(94)
    result = solver.solve()

    assert result["load_factor"] == pytest.approx(1.0)
    assert len(result["increments"]) == 20
    assert _node_result(solver, result["displacements"], 2)[0] == pytest.approx(20.0)
    assert result["reactions"][1][0] == pytest.approx(-4.2e6)
    assert result["stresses"][1]["axial_strain"] == pytest.approx(0.2)
    assert result["stresses"][1]["axial_stress"] == pytest.approx(42000.0)


def test_adaptive_control_changes_step_without_changing_final_equilibrium():
    solver = _solver(95)
    result = solver.solve()

    factors = np.array([increment.load_factor for increment in result["increments"]])
    step_sizes = np.diff(np.concatenate(([0.0], factors)))
    assert factors[-1] == pytest.approx(1.0)
    assert len(factors) < 20
    assert len(np.unique(np.round(step_sizes, 10))) > 1
    assert _node_result(solver, result["displacements"], 2)[0] == pytest.approx(20.0)


def test_torsion_and_axial_response_remain_uncoupled():
    solver = _solver(117)
    result = solver.solve()
    node_2 = _node_result(solver, result["displacements"], 2)

    assert node_2[0] == pytest.approx(20.0)
    assert node_2[3] == pytest.approx(0.9904762, rel=1e-7)
    assert result["stresses"][1]["torque"] == pytest.approx(40000.0)
    assert result["stresses"][1]["torsional_stress"] == pytest.approx(2000.0)


def test_force_control_limit_point_reports_last_converged_state():
    solver = _solver(120)

    with pytest.raises(NonlinearConvergenceError) as error:
        solver.solve()

    assert 0.95 < error.value.last_converged_load_factor < 1.0
    assert error.value.cutbacks > 0
    assert solver.increment_history


def test_assembled_tangent_and_internal_force_are_symmetric_and_equilibrated():
    solver = _solver(112)
    result = solver.solve()
    state = solver.assemble_state(result["displacements"])
    free_dofs, _, _ = solver.get_dof_partition()
    residual = solver.external_load - state.internal_force

    assert state.tangent_stiffness.toarray() == pytest.approx(
        state.tangent_stiffness.toarray().T
    )
    assert np.linalg.norm(residual[free_dofs], ord=np.inf) < 1.0e-5
    assert state.strain_energy == pytest.approx(result["strain_energy"])


def test_unsupported_convergence_without_force_residual_is_rejected():
    path = next((VALIDATION_ROOT / "case_94").glob("*.dat"))
    model = parse_nastran_input(str(path))
    model.analysis.nonlinear_step.convergence = "U"

    with pytest.raises(NotImplementedError, match="requires force-residual"):
        NonlinearStaticSolver(model)
