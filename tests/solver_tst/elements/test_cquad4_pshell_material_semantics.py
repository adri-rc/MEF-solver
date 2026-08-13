import numpy as np

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver


def test_blank_mid2_produces_no_cquad4_plate_stiffness(tmp_path):
    """A PSHELL with only MID1 is membrane-only, including its inactive DOFs."""
    deck = tmp_path / "blank_mid2_cquad4.dat"
    deck.write_text(
        """SOL 101
CEND
BEGIN BULK
PARAM,K6ROT,0.0
MAT1,1,2.1E5,,0.3
PSHELL,1,1,1.0
GRID,1,,0.,0.,0.
GRID,2,,2.,0.,0.
GRID,3,,2.,1.,0.
GRID,4,,0.,1.,0.
CQUAD4,1,1,1,2,3,4
ENDDATA
"""
    )

    model = parse_nastran_input(str(deck))
    assert model.properties[1].mat_id_bending is None

    solver = LinearStaticSolver(model)
    solver.build_stiffness_matrix()
    stiffness = solver.K_global.toarray()
    plate_dofs = [2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22]
    np.testing.assert_allclose(
        stiffness[np.ix_(plate_dofs, plate_dofs)], 0.0, atol=1.0e-12
    )
