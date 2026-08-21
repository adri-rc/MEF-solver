"""Unit tests for the corotational nonlinear CBAR formulation."""

import numpy as np
import pytest

from src.preprocessing.data_model import BarProperty, IsotropicMaterial, Node
from src.solver.elements.bar_element import BarElementSolver
from src.solver.elements.cbar_non_linear import CBARNonlinear


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


def _element(*, pa=0, pb=0, wa=None, wb=None, large_displacements=True):
    node_a = _node(1, (0.0, 0.0, 0.0))
    node_b = _node(2, (1000.0, 0.0, 0.0))
    prop = BarProperty(
        pid=1,
        mat_id=1,
        A=1000.0,
        I1=8.0e6,
        I2=2.0e6,
        J=1.5e6,
        stress_points=[(25.0, 40.0), (-25.0, 40.0), (-25.0, -40.0), (25.0, -40.0)],
    )
    material = IsotropicMaterial(mid=1, E=210000.0, nu=0.3)
    element = CBARNonlinear(
        node_a,
        node_b,
        prop,
        material,
        orientation=[0.0, 1.0, 0.0],
        pa=pa,
        pb=pb,
        wa=wa,
        wb=wb,
        large_displacements=large_displacements,
    )
    return element, node_a, node_b, prop, material


def test_undeformed_tangent_matches_the_established_linear_cbar():
    element, node_a, node_b, prop, material = _element()
    state = element.evaluate(np.zeros(12))
    linear = BarElementSolver(
        node_a, node_b, prop, material, orientation=[0.0, 1.0, 0.0]
    ).get_stiffness_matrix()

    assert state.internal_force == pytest.approx(np.zeros(12), abs=1.0e-8)
    assert state.tangent_stiffness == pytest.approx(linear, rel=2.0e-5, abs=2.0e-2)


def test_finite_rigid_rotation_is_objective():
    element, *_ = _element()
    displacement = np.zeros(12)
    displacement[6:9] = [-1000.0, 1000.0, 0.0]
    displacement[3:6] = [0.0, 0.0, np.pi / 2.0]
    displacement[9:12] = [0.0, 0.0, np.pi / 2.0]

    state = element.evaluate(displacement)

    assert state.current_length == pytest.approx(1000.0)
    assert state.current_axes[:, 0] == pytest.approx([0.0, 1.0, 0.0], abs=1.0e-12)
    assert state.basic_deformations == pytest.approx(np.zeros(6), abs=2.0e-12)
    assert state.internal_force == pytest.approx(np.zeros(12), abs=1.0e-5)


def test_finite_axial_extension_uses_nastran_logarithmic_strain():
    element, *_ = _element()
    displacement = np.zeros(12)
    displacement[6] = 1000.0 * (np.exp(0.02) - 1.0)

    state = element.evaluate(displacement)

    assert state.current_length == pytest.approx(1000.0 * np.exp(0.02))
    assert state.axial_strain == pytest.approx(0.02)
    assert state.axial_stress == pytest.approx(4200.0)
    assert state.axial_force == pytest.approx(4.2e6)
    assert state.internal_force[[0, 6]] == pytest.approx([-4.2e6, 4.2e6])


def test_pure_torsion_uses_gj_over_reference_length():
    element, *_ = _element()
    displacement = np.zeros(12)
    displacement[9] = 1.0

    state = element.evaluate(displacement)

    expected_torque = element.mat.G * element.prop.J / element.reference_length
    assert state.twist_rate == pytest.approx(1.0 / 1000.0)
    assert state.torque == pytest.approx(expected_torque)
    assert state.internal_force[[3, 9]] == pytest.approx(
        [-expected_torque, expected_torque]
    )


def test_tangent_matches_spatial_internal_force_derivative():
    element, *_ = _element()
    displacement = np.array(
        [
            1.2, -0.7, 0.4, 0.03, -0.02, 0.01,
            8.0, 32.0, -14.0, 0.16, 0.05, 0.09,
        ]
    )
    state = element.evaluate(displacement)

    numerical = np.zeros((12, 12))
    for column in range(12):
        step = 1.0e-5 if column % 6 < 3 else 1.0e-7
        plus = displacement.copy()
        minus = displacement.copy()
        plus[column] += step
        minus[column] -= step
        numerical[:, column] = (
            element.evaluate_internal_force(plus)
            - element.evaluate_internal_force(minus)
        ) / (2.0 * step)

    assert state.tangent_stiffness == pytest.approx(
        numerical, rel=2.0e-4, abs=5.0e2
    )


def test_released_end_rotation_carries_no_spurious_action():
    element, *_ = _element(pb=6)
    displacement = np.zeros(12)
    displacement[11] = 0.2

    state = element.evaluate(displacement)

    assert state.internal_force == pytest.approx(np.zeros(12), abs=1.0e-7)


def test_rigid_rotation_with_large_rotation_offsets_is_objective():
    offset = [0.0, 40.0, 0.0]
    element, *_ = _element(wa=offset, wb=offset)
    displacement = np.zeros(12)
    displacement[6:9] = [-1000.0, 1000.0, 0.0]
    displacement[3:6] = [0.0, 0.0, np.pi / 2.0]
    displacement[9:12] = [0.0, 0.0, np.pi / 2.0]

    state = element.evaluate(displacement)

    assert state.basic_deformations == pytest.approx(np.zeros(6), abs=2.0e-12)
    assert state.internal_force == pytest.approx(np.zeros(12), abs=1.0e-5)
