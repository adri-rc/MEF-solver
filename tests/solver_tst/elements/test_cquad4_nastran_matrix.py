from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver


VALIDATION_ROOT = Path(__file__).parents[1] / "validation_cases"


@pytest.mark.parametrize(
    "folder, deck, reference, active_dofs",
    [
        (
            "case_51",
            "ej51_kmembrane_k6rot_0.dat",
            "ej51_kmembrane_k6rot_0_Nastran.csv",
            [0, 1, 6, 7, 12, 13, 18, 19],
        ),
        (
            "case_52",
            "ej52_k_pure_flex_mid3_0_k6rot_0.dat",
            "ej52_k_pure_flex_mid3_0_k6rot_0_Nastran.csv",
            [2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22],
        ),
    ],
)
def test_cquad4_active_stiffness_block_against_nastran(
    folder, deck, reference, active_dofs
):
    case_path = VALIDATION_ROOT / folder
    model = parse_nastran_input(str(case_path / deck))
    solver = LinearStaticSolver(model)
    solver.build_stiffness_matrix()

    calculated = solver.K_global.toarray()[np.ix_(active_dofs, active_dofs)]
    nastran_full = pd.read_csv(case_path / reference).iloc[:, 2:].to_numpy(float)
    nastran = nastran_full[np.ix_(active_dofs, active_dofs)]

    relative_frobenius_error = np.linalg.norm(calculated - nastran) / np.linalg.norm(nastran)
    assert relative_frobenius_error < 0.03
