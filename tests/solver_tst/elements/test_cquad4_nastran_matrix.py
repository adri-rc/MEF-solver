"""Regression tests for CQUAD4 stiffness matrices exported by Nastran.

Each validation case contains a Nastran input deck and the corresponding KAAX
DMIG converted to CSV.  Only physically active degrees of freedom are compared
for planar cases; the inclined three-dimensional case compares all 24 DOFs.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver


VALIDATION_ROOT = Path(__file__).parents[1] / "validation_cases"


@pytest.mark.parametrize(
    "folder, deck, reference, active_dofs, tolerance",
    [
        pytest.param(
            "case_51",
            "ej51_kmembrane_k6rot_0.dat",
            "ej51_kmembrane_k6rot_0_Nastran.csv",
            [0, 1, 6, 7, 12, 13, 18, 19],
            0.01,  # Maximum relative Frobenius error: 1% for flat membrane.
            id="case_51_membrane",
        ),
        pytest.param(
            "case_52",
            "ej52_k_pure_flex_mid3_0_k6rot_0.dat",
            "ej52_k_pure_flex_mid3_0_k6rot_0_Nastran.csv",
            [2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22],
            0.01,  # Maximum relative Frobenius error: 1% for flat bending.
            id="case_52_bending",
        ),
        pytest.param(
            "case_53",
            "ej53_k_memb-flex_mid3_0_k6rot_0.dat",
            "ej53_k_memb-flex_mid3_0_k6rot_0_Nastran.csv",
            [
                0, 1, 2, 3, 4,       # Node 1: drilling rotation is inactive.
                6, 7, 8, 9, 10,      # Node 2.
                12, 13, 14, 15, 16,  # Node 3.
                18, 19, 20, 21, 22,  # Node 4.
            ],
            0.01,  # Maximum relative Frobenius error: 1% for flat combined K.
            id="case_53_membrane_bending",
        ),
        pytest.param(
            "case_54",
            "ej54_k_memb-flex-3d_mid3_0_k6rot_0.dat",
            "ej54_k_memb-flex-3d_mid3_0_k6rot_0_Nastran.csv",
            # The element is inclined in 3D, so local membrane/bending terms
            # contribute to every translational and rotational global DOF.
            list(range(24)),
            0.03,  # Maximum relative Frobenius error: 3% for warped full K.
            id="case_54_warped_full_matrix",
        ),
    ],
)
def test_cquad4_active_stiffness_block_against_nastran(
    folder, deck, reference, active_dofs, tolerance
):
    """Compare the assembled active stiffness block with Nastran's KAAX DMIG."""
    case_path = VALIDATION_ROOT / folder
    model = parse_nastran_input(str(case_path / deck))
    solver = LinearStaticSolver(model)
    solver.build_stiffness_matrix()

    calculated = solver.K_global.toarray()[np.ix_(active_dofs, active_dofs)]
    nastran_full = pd.read_csv(case_path / reference).iloc[:, 2:].to_numpy(float)
    nastran = nastran_full[np.ix_(active_dofs, active_dofs)]

    relative_frobenius_error = np.linalg.norm(
        calculated - nastran
    ) / np.linalg.norm(nastran)

    assert relative_frobenius_error < tolerance
