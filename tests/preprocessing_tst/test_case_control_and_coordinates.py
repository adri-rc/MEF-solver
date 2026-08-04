import numpy as np
import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver


def test_subcase_selection_load_combination_and_spcadd(tmp_path):
    deck = """SOL 101
CEND
SUBCASE 1
  SPC = 100
  LOAD = 10
SUBCASE 2
  SPC = 201
  LOAD = 30
BEGIN BULK
GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
MAT1,1,1000.,,0.3
PROD,1,1,1.0
CROD,1,1,1,2
SPC1,100,123456,1
SPC1,200,123456,2
SPCADD,201,200
FORCE,10,2,,10.,1.,0.,0.
FORCE,20,2,,10.,1.,0.,0.
LOAD,30,2.0,3.0,20
ENDDATA
"""
    path = tmp_path / "subcases.dat"
    path.write_text(deck)

    with pytest.raises(ValueError, match="multiple analysis subcases"):
        parse_nastran_input(str(path))

    model = parse_nastran_input(str(path), subcase_id=2)
    assert model.active_subcase_id == 2
    assert model.active_load_sid == 30
    assert model.active_spc_sid == 201
    assert len(model.loads) == 1
    assert model.loads[0].module == pytest.approx(60.0)
    assert {(constraint.node_id, constraint.components) for constraint in model.constraints} == {
        (2, "123456")
    }


def test_grid_cd_force_cid_and_permanent_ps(tmp_path):
    deck = """CORD2R,10,0,0.,0.,0.,0.,0.,1.
,0.,1.,0.
GRID,1,,0.,0.,0.,10,123456
GRID,2,,0.,1.,0.,10
MAT1,1,1000.,,0.3
PROD,1,1,1.0
CROD,1,1,1,2
SPC1,1,23456,2
FORCE,1,2,10,10.,1.,0.,0.
"""
    path = tmp_path / "coordinates.bdf"
    path.write_text(deck)

    model = parse_nastran_input(str(path))
    np.testing.assert_allclose(model.loads[0].vector, [0.0, 10.0 / 10.0, 0.0], atol=1e-12)
    assert any(c.sid == 0 and c.node_id == 1 for c in model.constraints)

    solver = LinearStaticSolver(model)
    solver.assemble()
    results = solver.solve()
    node_2 = solver.dof_mapper.get_dof_indices(2)
    # Results are expressed in GRID.CD. Local T1 is aligned with basic +Y.
    assert results["displacements"][node_2[0]] == pytest.approx(0.01)
    assert results["reactions"][1][0] == pytest.approx(-10.0)
