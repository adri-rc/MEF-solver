"""Tests for linear CROD stiffness and result recovery."""

import numpy as np
import pytest

from src.solver.elements.rod_element import RodElementSolver


class MockNode:
    def __init__(self, nid, coords):
        self.nid = nid
        self._coords = np.array(coords, dtype=float)
    def get_coords(self):
        return self._coords

class MockProp:
    def __init__(self, A, J, C=0.0):
        self.A = A
        self.J = J
        self.C = C

class MockMat:
    def __init__(self, E, G):
        self.E = E
        self.G = G

class MockDofMapper:
    """Provide the contiguous six-DOF layout used by these two-node models."""

    def get_dof_indices(self, node_id):
        start = (node_id - 1) * 6
        return list(range(start, start + 6))


class TestRodElementSolver:
    def test_stiffness_axial_x_aligned(self):
        """Check EA/L stiffness and decoupling for a rod along global X."""
        n1 = MockNode(1, [0.0, 0.0, 0.0])
        n2 = MockNode(2, [2.0, 0.0, 0.0])

        prop = MockProp(A=0.5, J=0.0)
        mat = MockMat(E=1000.0, G=0.0)
        expected_k = 250.0

        solver = RodElementSolver(n1, n2, prop, mat, n_gaussp=1)
        K = solver.get_stiffness_matrix()

        assert K[0, 0] == pytest.approx(expected_k)
        assert K[0, 6] == pytest.approx(-expected_k)
        assert K[6, 0] == pytest.approx(-expected_k)
        assert K[6, 6] == pytest.approx(expected_k)

        axial_indices = [0, 6]

        for row in axial_indices:
            for col in range(12):
                if col not in axial_indices:
                    assert K[row, col] == pytest.approx(0.0, abs=1e-9), \
                        f"Unexpected coupling between Axial (idx {row}) and DOF {col}"

    def test_stiffness_torsional_x_aligned(self):
        """Check GJ/L stiffness and decoupling for torsion about global X."""
        n1 = MockNode(1, [0.0, 0.0, 0.0])
        n2 = MockNode(2, [1.0, 0.0, 0.0])

        prop = MockProp(A=0.0, J=0.1)
        mat = MockMat(E=0.0, G=800.0)
        expected_k = 80.0

        solver = RodElementSolver(n1, n2, prop, mat)
        K = solver.get_stiffness_matrix()

        assert K[3, 3] == pytest.approx(expected_k)
        assert K[3, 9] == pytest.approx(-expected_k)
        assert K[9, 9] == pytest.approx(expected_k)

        torsion_indices = [3, 9]

        for row in torsion_indices:
            for col in range(12):
                if col not in torsion_indices:
                    assert K[row, col] == pytest.approx(0.0, abs=1e-9), \
                        f"Unexpected coupling between Torsion (idx {row}) and DOF {col}"

    def test_transformation_45_degrees(self):
        """
        A 45-degree rod projects half of its axial stiffness into each global
        X/Y direct and coupling term.
        """
        n1 = MockNode(1, [0.0, 0.0, 0.0])
        n2 = MockNode(2, [1.0, 1.0, 0.0])

        L = np.sqrt(2.0)
        E_val = 100.0
        A_val = 1.0
        k_local = (E_val * A_val) / L

        prop = MockProp(A=A_val, J=0.0)
        mat = MockMat(E=E_val, G=0.0)

        solver = RodElementSolver(n1, n2, prop, mat, n_gaussp=1)
        K = solver.get_stiffness_matrix()

        k_projected = k_local * 0.5

        assert K[0, 0] == pytest.approx(k_projected)
        assert K[1, 1] == pytest.approx(k_projected)
        assert K[0, 1] == pytest.approx(k_projected)
        assert K[0, 6] == pytest.approx(-k_projected)
        assert K[0, 7] == pytest.approx(-k_projected)

    def test_zero_length_error(self):
        """Reject coincident rod nodes with a diagnostic message."""
        n1 = MockNode(1, [1.0, 1.0, 1.0])
        n2 = MockNode(2, [1.0, 1.0, 1.0])
        prop = MockProp(A=1.0, J=1.0)
        mat = MockMat(E=1.0, G=1.0)

        with pytest.raises(ValueError) as excinfo:
            RodElementSolver(n1, n2, prop, mat)
        error_msg = str(excinfo.value)
        assert "Zero length" in error_msg
        assert "nodes 1 and 2" in error_msg

    def test_calculate_stresses_axial_tension(self):
        """Recover axial stress and force for an X-aligned rod in tension."""
        n1 = MockNode(1, [0., 0., 0.])
        n2 = MockNode(2, [2., 0., 0.])
        prop = MockProp(A=10.0, J=1.0)
        mat = MockMat(E=2000.0, G=1000.0)

        solver = RodElementSolver(n1, n2, prop, mat)
        mapper = MockDofMapper()

        u_global = np.zeros(12)

        u_global[6] = 0.1

        results = solver.calculate_stresses(u_global, mapper)

        assert results['axial_stress'] == pytest.approx(100.0)
        assert results['axial_force'] == pytest.approx(1000.0)
        assert results['torque'] == pytest.approx(0.0)

    def test_calculate_stresses_rotated_rod(self):
        """
        A global Y displacement is axial for a Y-aligned rod; global X is
        transverse and produces no stress in the CROD formulation.
        """
        n1 = MockNode(1, [0., 0., 0.])
        n2 = MockNode(2, [0., 1., 0.])
        prop = MockProp(A=1.0, J=1.0)
        mat = MockMat(E=1000.0, G=1000.0)

        solver = RodElementSolver(n1, n2, prop, mat)
        mapper = MockDofMapper()

        u_global = np.zeros(12)

        u_global[7] = 0.1

        results = solver.calculate_stresses(u_global, mapper)

        assert results['axial_stress'] == pytest.approx(100.0)

        u_global_shear = np.zeros(12)
        u_global_shear[6] = 0.1

        results_shear = solver.calculate_stresses(u_global_shear, mapper)
        assert results_shear['axial_stress'] == pytest.approx(0.0, abs=1e-9)

    def test_calculate_stresses_rotated_rod_torsion(self):
        """
        For a Y-aligned rod, global Ry is local twist; global Rx is transverse
        and produces no stress because CROD has no bending stiffness.
        """
        n1 = MockNode(1, [0., 0., 0.])
        n2 = MockNode(2, [0., 1., 0.])

        prop = MockProp(A=10.0, J=2.0, C=0.5)
        mat = MockMat(E=1000.0, G=1000.0)

        solver = RodElementSolver(n1, n2, prop, mat)
        mapper = MockDofMapper()

        u_global_tor = np.zeros(12)
        u_global_tor[10] = 0.1

        results = solver.calculate_stresses(u_global_tor, mapper)

        assert results['torque'] == pytest.approx(200.0)
        assert results['torsional_stress'] == pytest.approx(50.0)

        u_global_bend = np.zeros(12)
        u_global_bend[9] = 0.1

        results_bend = solver.calculate_stresses(u_global_bend, mapper)

        assert results_bend['torque'] == pytest.approx(0.0, abs=1e-9)
        assert results_bend['torsional_stress'] == pytest.approx(0.0, abs=1e-9)

    def test_stress_recovery_axial_and_strain(self):
        """Recover strain, stress, and force from a 0.01 m extension."""
        n1, n2 = MockNode(1, [0,0,0]), MockNode(2, [1,0,0])
        prop = MockProp(A=0.1, J=0.0)
        mat = MockMat(E=200e9, G=0)

        solver = RodElementSolver(n1, n2, prop, mat)
        mapper = MockDofMapper()

        u_global = np.zeros(12)
        u_global[6] = 0.01

        results = solver.calculate_stresses(u_global, mapper)

        expected_strain = 0.01 / 1.0
        expected_stress = expected_strain * 200e9

        assert results['axial_strain'] == pytest.approx(expected_strain)
        assert results['axial_stress'] == pytest.approx(expected_stress)
        assert results['axial_force'] == pytest.approx(expected_stress * 0.1)

    def test_stress_recovery_torsion_and_twist(self):
        """Recover twist rate and torque for a Y-aligned rod."""
        n1, n2 = MockNode(1, [0,0,0]), MockNode(2, [0,2.0,0])
        prop = MockProp(A=0.1, J=0.5)
        mat = MockMat(E=0, G=80e9)

        solver = RodElementSolver(n1, n2, prop, mat)
        mapper = MockDofMapper()

        # Global Ry is local Rx for this rod.
        u_global = np.zeros(12)
        u_global[10] = 0.1

        results = solver.calculate_stresses(u_global, mapper)

        expected_twist = 0.1 / 2.0
        expected_torque = 80e9 * 0.5 * expected_twist

        assert results['twist_rate'] == pytest.approx(expected_twist)
        assert results['torque'] == pytest.approx(expected_torque)
        assert results['axial_force'] == pytest.approx(0.0)
