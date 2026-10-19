"""Intermediate SOL 400 references written explicitly from MSC Nastran."""

from pathlib import Path

import numpy as np
import pytest

from src.preprocessing.nastran_reader import parse_nastran_input
from src.solver.nonlinear_analysis import NonlinearStaticSolver


VALIDATION_ROOT = Path(__file__).parent / "validation_cases"

# CASE 112 PCH: converged states from the 40 fixed increments.  The test never
# reads the PCH; these values are deliberately auditable reference data.
CASE_112_INTERMEDIATE = {
    0.25: {"uy": -3.458443, "axial_stress": -3.599418e3},
    0.50: {"uy": -7.129332, "axial_stress": -7.347492e3},
    0.75: {"uy": -11.07308, "axial_stress": -1.128630e4},
    1.00: {"uy": -15.38383, "axial_stress": -1.548118e4},
}

# CASES 114/115: (Ux(+1000)-Ux(-1000))/2000 around CASE 112.
CASE_112_NASTRAN_TANGENT_COMPLIANCE_XX = 6.127186e-6


def _case_112_solver():
    path = next((VALIDATION_ROOT / "case_112").glob("*.dat"))
    return NonlinearStaticSolver(parse_nastran_input(str(path)))


def test_converged_increment_history_matches_nastran_load_path():
    solver = _case_112_solver()
    result = solver.solve()
    by_factor = {
        round(increment.load_factor, 8): increment
        for increment in result["increments"]
    }

    for load_factor, expected in CASE_112_INTERMEDIATE.items():
        increment = by_factor[load_factor]
        node_3 = increment.displacements[solver.dof_mapper.get_dof_indices(3)]
        element_state = solver.element_solvers[1].evaluate(
            solver.get_element_displacements_basic(
                solver.model.elements[1], increment.displacements
            )
        )
        assert node_3[1] == pytest.approx(expected["uy"], rel=2e-6, abs=2e-6)
        assert element_state.axial_stress == pytest.approx(
            expected["axial_stress"], rel=2e-6, abs=2e-3
        )


def test_deformed_global_tangent_matches_nastran_central_probe():
    solver = _case_112_solver()
    result = solver.solve()
    state = solver.assemble_state(result["displacements"])
    node_3_dofs = solver.dof_mapper.get_dof_indices(3)
    active = node_3_dofs[:2]
    reduced_tangent = state.tangent_stiffness[active, :][:, active].toarray()
    compliance = np.linalg.inv(reduced_tangent)

    assert compliance[0, 0] == pytest.approx(
        CASE_112_NASTRAN_TANGENT_COMPLIANCE_XX,
        rel=2e-6,
        abs=1e-12,
    )
