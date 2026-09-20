import numpy as np
import pytest

from src.preprocessing.data_model import IsotropicMaterial, Node, RodProperty
from src.solver.elements.crod_nonlinear import CRODNonlinear


def _node(nid, xyz):
    return Node(
        nid=nid,
        coord_id=0,
        x1=xyz[0],
        x2=xyz[1],
        x3=xyz[2],
        cd_id=0,
        ps="",
        seid=0,
        x=xyz[0],
        y=xyz[1],
        z=xyz[2],
    )


def _element(*, j=0.0, c=0.0, large_displacements=True):
    return CRODNonlinear(
        _node(1, (0.0, 0.0, 0.0)),
        _node(2, (100.0, 0.0, 0.0)),
        RodProperty(pid=1, mat_id=1, A=100.0, J=j, C=c, NSM=0.0),
        IsotropicMaterial(mid=1, E=210000.0, nu=0.3),
        large_displacements=large_displacements,
    )


def test_finite_axial_extension_uses_nastran_engineering_strain():
    element = _element()
    displacement = np.zeros(12)
    displacement[6] = 20.0

    state = element.evaluate(displacement)

    assert state.current_length == pytest.approx(120.0)
    assert state.axial_strain == pytest.approx(0.2)
    assert state.axial_stress == pytest.approx(42000.0)
    assert state.axial_force == pytest.approx(4.2e6)
    assert state.internal_force[[0, 6]] == pytest.approx([-4.2e6, 4.2e6])

    # Material stiffness acts axially; the current axial force supplies the
    # differential/geometric stiffness in both transverse directions.
    current_block = state.tangent_stiffness[:3, :3]
    assert np.diag(current_block) == pytest.approx([210000.0, 35000.0, 35000.0])


def test_rigid_rotation_changes_direction_without_creating_axial_force():
    element = _element()
    displacement = np.zeros(12)
    displacement[6:9] = [-100.0, 100.0, 0.0]

    state = element.evaluate(displacement)

    assert state.current_direction == pytest.approx([0.0, 1.0, 0.0])
    assert state.axial_strain == pytest.approx(0.0, abs=1e-14)
    assert state.internal_force == pytest.approx(np.zeros(12), abs=1e-8)
    assert state.tangent_stiffness[1, 1] == pytest.approx(210000.0)
    assert state.tangent_stiffness[0, 0] == pytest.approx(0.0, abs=1e-10)


def test_consistent_tangent_matches_central_difference_of_internal_force():
    element = _element(j=50.0, c=2.5)
    displacement = np.zeros(12)
    displacement[0:3] = [1.2, -0.7, 0.4]
    displacement[6:9] = [8.0, 12.0, -4.0]
    displacement[3] = -0.1
    displacement[9] = 0.3
    state = element.evaluate(displacement)

    step = 1.0e-5
    numerical = np.zeros((12, 12))
    for column in range(12):
        plus = displacement.copy()
        minus = displacement.copy()
        plus[column] += step
        minus[column] -= step
        numerical[:, column] = (
            element.evaluate(plus).internal_force
            - element.evaluate(minus).internal_force
        ) / (2.0 * step)

    assert state.tangent_stiffness == pytest.approx(numerical, rel=2e-8, abs=2e-5)
    assert state.tangent_stiffness == pytest.approx(state.tangent_stiffness.T)


def test_reference_axis_torsion_matches_cases_116_and_117():
    element = _element(j=50.0, c=2.5)
    displacement = np.zeros(12)
    displacement[9] = 0.9904762

    state = element.evaluate(displacement)

    assert state.twist_rate == pytest.approx(0.009904762)
    assert state.torque == pytest.approx(40000.0, rel=2e-8)
    assert state.torsional_stress == pytest.approx(2000.0, rel=2e-8)
    assert state.internal_force[[3, 9]] == pytest.approx([-40000.0, 40000.0], rel=2e-8)


def test_lgdisp_disabled_keeps_reference_direction_and_linear_stiffness():
    element = _element(large_displacements=False)
    displacement = np.zeros(12)
    displacement[7] = 25.0

    state = element.evaluate(displacement)

    assert state.axial_strain == pytest.approx(0.0)
    assert state.axial_force == pytest.approx(0.0)
    assert state.tangent_stiffness[0, 0] == pytest.approx(210000.0)
    assert state.tangent_stiffness[1, 1] == pytest.approx(0.0)


def test_collapsed_current_geometry_is_rejected():
    element = _element()
    displacement = np.zeros(12)
    displacement[6] = -100.0

    with pytest.raises(ValueError, match="collapsed current length"):
        element.evaluate(displacement)
