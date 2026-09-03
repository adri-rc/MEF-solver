from pathlib import Path

import numpy as np
import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.elements.ctria3_element import CTRIA3ElementSolver
from src.solver.linear_static_solver import LinearStaticSolver
from tests.solver_tst.elements.CTRIA3_pch_reader import leer_dmig_pch


CASE_DIR = (
    Path(__file__).resolve().parents[1]
    / "validation_cases"
    / "case_ctria3_pshell_material_semantics"
)
MEMBRANE_DOFS = [0, 1, 6, 7, 12, 13]
PLATE_DOFS = [2, 3, 4, 8, 9, 10, 14, 15, 16]


def _read_model(case_name: str):
    return parse_nastran_input(str(CASE_DIR / f"{case_name}.dat"))


def _build_global_stiffness(case_name: str) -> tuple[object, np.ndarray]:
    model = _read_model(case_name)
    solver = LinearStaticSolver(model)
    solver.build_stiffness_matrix()
    return model, solver.K_global.toarray()


def _make_element_solver(case_name: str) -> CTRIA3ElementSolver:
    model = _read_model(case_name)
    element = model.elements[100]
    prop = model.properties[element.prop_id]
    nodes = [model.nodes[nid] for nid in element.nodes]
    return CTRIA3ElementSolver(
        *nodes,
        prop,
        model.materials.get(prop.mat_id),
        mat_bending=model.materials.get(prop.mat_id_bending),
        mat_shear=model.materials.get(prop.mat_id_shear),
        k6rot=0.0,
        theta=element.theta,
    )


def _block(matrix: np.ndarray, dofs: list[int]) -> np.ndarray:
    return matrix[np.ix_(dofs, dofs)]


def test_nastran_reference_files_record_expected_outcomes():
    valid_cases = (
        "ctria3_01_mid1_only",
        "ctria3_02_mid1_mid2",
        "ctria3_03_mid1_mid2_mid3",
        "ctria3_04_mid2_only",
        "ctria3_07_mid2_mid3",
    )
    invalid_cases = (
        "ctria3_05_mid3_only_invalid",
        "ctria3_06_mid1_mid3_invalid",
    )

    for case_name in valid_cases:
        f06_text = (CASE_DIR / f"{case_name}.f06").read_text(
            encoding="latin-1"
        )
        assert "USER FATAL MESSAGE" not in f06_text
        assert "* * * END OF JOB * * *" in f06_text
        assert (CASE_DIR / f"{case_name}.pch").is_file()

    for case_name in invalid_cases:
        f06_text = (CASE_DIR / f"{case_name}.f06").read_text(
            encoding="latin-1"
        )
        assert "USER FATAL MESSAGE 9994 (BULAEM)" in f06_text
        assert "MID3 = 0 when MID2 =  0" in f06_text
        assert not (CASE_DIR / f"{case_name}.pch").exists()


def test_mid1_only_is_membrane_only_and_matches_nastran_membrane():
    _, stiffness = _build_global_stiffness("ctria3_01_mid1_only")
    nastran = leer_dmig_pch(
        CASE_DIR / "ctria3_01_mid1_only.pch", num_nudos=3
    )

    np.testing.assert_allclose(_block(nastran, PLATE_DOFS), 0.0, atol=0.0)
    np.testing.assert_allclose(_block(stiffness, PLATE_DOFS), 0.0, atol=1.0e-12)
    np.testing.assert_allclose(
        _block(stiffness, MEMBRANE_DOFS),
        _block(nastran, MEMBRANE_DOFS),
        rtol=1.0e-9,
        atol=1.0e-5,
    )


def test_mid2_only_is_accepted_without_creating_membrane_stiffness():
    model, stiffness = _build_global_stiffness("ctria3_04_mid2_only")
    prop = model.properties[10]
    nastran = leer_dmig_pch(
        CASE_DIR / "ctria3_04_mid2_only.pch", num_nudos=3
    )

    assert prop.mat_id is None
    assert prop.mat_id_bending == 2
    assert prop.mat_id_shear is None
    np.testing.assert_allclose(_block(stiffness, MEMBRANE_DOFS), 0.0, atol=1.0e-12)
    np.testing.assert_allclose(_block(nastran, MEMBRANE_DOFS), 0.0, atol=0.0)
    assert np.linalg.norm(_block(stiffness, PLATE_DOFS)) > 0.0
    assert np.linalg.norm(_block(nastran, PLATE_DOFS)) > 0.0


def test_blank_mid3_is_not_replaced_with_mid2():
    element_solver = _make_element_solver("ctria3_02_mid1_mid2")

    assert element_solver.mat is not None
    assert element_solver.mat_bending is not None
    assert element_solver.mat_shear is None
    assert element_solver.mat.E == pytest.approx(70000.0)
    assert element_solver.mat_bending.E == pytest.approx(140000.0)
    np.testing.assert_allclose(element_solver._get_D_shear(), 0.0, atol=0.0)


def test_all_three_materials_remain_independent():
    element_solver = _make_element_solver("ctria3_03_mid1_mid2_mid3")

    assert element_solver.mat.E == pytest.approx(70000.0)
    assert element_solver.mat_bending.E == pytest.approx(140000.0)
    assert element_solver.mat_shear.E == pytest.approx(28000.0)
    assert np.linalg.norm(element_solver._get_D_bending()) > 0.0
    assert np.linalg.norm(element_solver._get_D_shear()) > 0.0


def test_mid2_mid3_without_mid1_is_accepted_as_plate_only():
    model, stiffness = _build_global_stiffness("ctria3_07_mid2_mid3")
    prop = model.properties[10]

    assert prop.mat_id is None
    assert prop.mat_id_bending == 2
    assert prop.mat_id_shear == 3
    np.testing.assert_allclose(_block(stiffness, MEMBRANE_DOFS), 0.0, atol=1.0e-12)
    assert np.linalg.norm(_block(stiffness, PLATE_DOFS)) > 0.0


@pytest.mark.parametrize(
    "case_name",
    ("ctria3_05_mid3_only_invalid", "ctria3_06_mid1_mid3_invalid"),
)
def test_mid3_without_mid2_is_rejected_like_nastran(case_name):
    with pytest.raises(ValueError, match="MID3.*MID2"):
        _read_model(case_name)
