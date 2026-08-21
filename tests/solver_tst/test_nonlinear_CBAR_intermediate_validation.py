"""Intermediate SOL 400 CBAR references transcribed from MSC Nastran.

The test deliberately does not read PCH/OP2 files.  Displacements come from
the converged ``DISPLACEMENTS`` blocks at the stated ``LOAD FACTOR``.  The
tangent reference is the central probe formed by cases 138 and 139 around
case 137: ``(T2(+100)-T2(-100))/(200*load_factor)``.
"""

from pathlib import Path

import numpy as np
import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.nonlinear_analysis import NonlinearStaticSolver


VALIDATION_ROOT = Path(__file__).parent / "validation_cases"

CASE_132_INTERMEDIATE_DISPLACEMENTS = {
    0.2: [-0.06745054, 11.9037247, 0.0, 0.0, 0.0, 0.01785595],
    0.4: [-0.26968592, 23.8012714, 0.0, 0.0, 0.0, 0.03570480],
    0.6: [-0.60635805, 35.6864853, 0.0, 0.0, 0.0, 0.05353947],
    0.8: [-1.07688892, 47.5532455, 0.0, 0.0, 0.0, 0.07135291],
    1.0: [-1.68047369, 59.3954926, 0.0, 0.0, 0.0, 0.08913815],
}

CASE_137_INTERMEDIATE_DISPLACEMENTS = {
    0.2: [0.47787645, 11.319180, -8.0521393, 0.04838159, 0.01007500, 0.01855162],
    0.4: [0.77574509, 22.283243, -16.147793, 0.09447157, 0.02091959, 0.03615898],
    0.6: [0.90810114, 32.913235, -24.083239, 0.13832106, 0.03215692, 0.05289228],
    0.8: [0.89248760, 43.220272, -31.732460, 0.18003432, 0.04353360, 0.06881820],
    1.0: [0.74746597, 53.210960, -39.022457, 0.21974313, 0.05488243, 0.08399875],
}

CASE_137_TANGENT_COMPLIANCE_T2 = {
    0.2: 1.9423961639404298e-4,
    0.4: 1.9044876098632811e-4,
    0.6: 1.8688837687174480e-4,
    0.8: 1.8343925476074218e-4,
    1.0: 1.8007278442382812e-4,
}


def _solve(case_id):
    deck = next((VALIDATION_ROOT / f"case_{case_id}").glob("*.dat"))
    solver = NonlinearStaticSolver(parse_nastran_input(str(deck)))
    result = solver.solve()
    increments = {
        round(increment.load_factor, 8): increment
        for increment in result["increments"]
    }
    return solver, increments


def test_cbar_planar_bending_converged_path_matches_nastran():
    solver, increments = _solve(132)
    node_dofs = solver.dof_mapper.get_dof_indices(2)

    for load_factor, expected in CASE_132_INTERMEDIATE_DISPLACEMENTS.items():
        computed = increments[load_factor].displacements[node_dofs]
        assert computed == pytest.approx(expected, rel=5.0e-4, abs=2.0e-5)


def test_cbar_spatial_combined_path_and_tangent_match_nastran_probes():
    solver, increments = _solve(137)
    node_dofs = solver.dof_mapper.get_dof_indices(2)

    for load_factor, expected in CASE_137_INTERMEDIATE_DISPLACEMENTS.items():
        increment = increments[load_factor]
        computed = increment.displacements[node_dofs]
        assert computed == pytest.approx(expected, rel=3.0e-2, abs=2.0e-5)

        assembly = solver.assemble_state(increment.displacements)
        reduced_tangent = assembly.tangent_stiffness[node_dofs, :][:, node_dofs]
        compliance_t2 = np.linalg.inv(reduced_tangent.toarray())[1, 1]
        assert compliance_t2 == pytest.approx(
            CASE_137_TANGENT_COMPLIANCE_T2[load_factor],
            rel=8.5e-3,
            abs=2.0e-9,
        )
