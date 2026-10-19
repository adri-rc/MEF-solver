"""Tests for linear CBAR stiffness and result recovery."""

import numpy as np
import pytest

from src.solver.elements.bar_element import BarElementSolver


class MockNode:
    def __init__(self, nid, coords):
        self.nid = nid
        self._coords = np.array(coords, dtype=float)
    def get_coords(self):
        return self._coords

class MockMat:
    def __init__(self, E, G):
        self.E = E
        self.G = G

class MockProp:
    def __init__(self, A, I1, I2, J, stress_points=None):
        self.pid = 1
        self.A, self.I1, self.I2, self.J = A, I1, I2, J
        self.stress_points = stress_points if stress_points is not None else []

class MockDofMapper:
    """Provide the contiguous six-DOF layout used by these two-node models."""

    def get_dof_indices(self, nid):
        start = (nid - 1) * 6
        return list(range(start, start + 6))


class TestBarElementSolver:
    def test_stiffness_matrix_bending_xy_analytical(self):
        """Compare XY bending terms with Euler-Bernoulli beam coefficients."""
        L = 2.0
        E = 1000.0
        I1 = 0.001

        n1 = MockNode(1, [0,0,0])
        n2 = MockNode(2, [L,0,0])
        prop = MockProp(A=0.1, I1=I1, I2=0.001, J=0.0)
        mat = MockMat(E=E, G=1000.0)
        orient = [0.0, 1.0, 0.0]

        solver = BarElementSolver(n1, n2, prop, mat, orient)
        K = solver.get_stiffness_matrix()

        # In the XY plane, Ty uses indices 1 and 7; Rz uses 5 and 11.
        k_force = (12 * E * I1) / L**3
        k_moment = (4 * E * I1) / L
        k_couple = (6 * E * I1) / L**2

        assert K[1, 1] == pytest.approx(k_force)
        assert K[7, 7] == pytest.approx(k_force)
        assert K[5, 5] == pytest.approx(k_moment)
        assert K[11, 11] == pytest.approx(k_moment)

        assert K[1, 5] == pytest.approx(k_couple)
        assert K[1, 7] == pytest.approx(-k_force)

    def test_stiffness_matrix_bending_xz_sign_convention(self):
        """Check the Tz-Ry sign convention for bending in the XZ plane."""
        L = 1.0
        E = 1.0
        I2 = 1.0

        n1, n2 = MockNode(1, [0,0,0]), MockNode(2, [1,0,0])
        prop = MockProp(A=0, I1=0, I2=I2, J=0)
        mat = MockMat(E=E, G=0)

        solver = BarElementSolver(n1, n2, prop, mat, [0,1,0])
        K = solver.get_stiffness_matrix()

        assert K[4, 4] == pytest.approx(4.0)

        # Tz-Ry has the opposite sign to Ty-Rz because theta_y = -dw/dx.
        assert K[2, 4] == pytest.approx(-6.0)

    def test_geometric_transformation_90_degrees(self):
        """A beam along global Y maps its axial stiffness to global Ty."""
        n1 = MockNode(1, [0,0,0])
        n2 = MockNode(2, [0, 10.0, 0])
        prop = MockProp(A=1.0, I1=0, I2=0, J=0)
        mat = MockMat(E=100.0, G=0)

        # The orientation vector cannot be parallel to the beam axis.
        orient = [1.0, 0.0, 0.0]
        solver = BarElementSolver(n1, n2, prop, mat, orient)
        K = solver.get_stiffness_matrix()

        assert K[1, 1] == pytest.approx(10.0)
        assert K[0, 0] == pytest.approx(0.0)

    def test_arbitrary_rotation_3d_kinematics(self):
        """
        Check that an axial global displacement on vector [1, 2, 3] produces
        the expected local strain and stress.
        """
        vec = np.array([1.0, 2.0, 3.0])
        L = np.linalg.norm(vec)
        unit_vec = vec / L

        n1 = MockNode(1, [0,0,0])
        n2 = MockNode(2, vec)
        prop = MockProp(A=0.1, I1=1e-4, I2=1e-4, J=1e-4)
        mat = MockMat(E=200e9, G=80e9)

        orient = [0.0, 1.0, 0.0]

        solver = BarElementSolver(n1, n2, prop, mat, orient)

        delta = 0.01
        disp_vec = unit_vec * delta

        u_global = np.zeros(12)
        u_global[6:9] = disp_vec

        mapper = MockDofMapper()
        res = solver.calculate_stresses(u_global, mapper)

        expected_strain = delta / L
        expected_stress = expected_strain * 200e9

        assert res['axial_strain'] == pytest.approx(expected_strain, rel=1e-5)
        assert res['axial_stress'] == pytest.approx(expected_stress, rel=1e-5)

        forces_n1 = res['forces_n1']
        assert abs(forces_n1['My']) < 1.0
        assert abs(res['torque']) < 1.0

    def test_stress_recovery_points_torsion(self):
        """Verify that torsional shear stress scales with section radius."""
        n1, n2 = MockNode(1, [0,0,0]), MockNode(2, [2,0,0])

        points = [(0.0, 0.0), (0.0, 0.05), (0.1, 0.0)]
        prop = MockProp(A=0.1, I1=1, I2=1, J=0.5, stress_points=points)

        mat = MockMat(E=0, G=1000.0)

        solver = BarElementSolver(n1, n2, prop, mat, [0,1,0])

        u_global = np.zeros(12)
        u_global[9] = 0.1

        mapper = MockDofMapper()
        results = solver.calculate_stresses(u_global, mapper)

        expected_twist = 0.05
        expected_torque = 1000.0 * 0.5 * expected_twist

        assert results['twist_rate'] == pytest.approx(expected_twist)
        assert results['torque'] == pytest.approx(expected_torque)

        stress_data = results['stress_points']

        assert stress_data['P0']['tau_xy'] == pytest.approx(0.0)
        assert stress_data['P1']['tau_xy'] == pytest.approx(2.5)
        assert stress_data['P2']['tau_xy'] == pytest.approx(5.0)

    def test_stress_recovery_bending_curvature(self):
        """Check bending-stress symmetry and axial/bending superposition."""
        n1, n2 = MockNode(1, [0,0,0]), MockNode(2, [1,0,0])

        points = [(0.0, 0.0), (0.1, 0.0), (-0.1, 0.0)]
        prop = MockProp(A=0.1, I1=1e-4, I2=1e-4, J=1e-4, stress_points=points)

        mat = MockMat(E=200e9, G=0)
        solver = BarElementSolver(n1, n2, prop, mat, [0,1,0])

        # Opposite end rotations produce constant positive curvature.
        theta = 0.001
        u_global = np.zeros(12)
        u_global[5] = -theta
        u_global[11] = theta

        mapper = MockDofMapper()
        res = solver.calculate_stresses(u_global, mapper)

        r_neutral = res['stress_points']['P0']
        r_top = res['stress_points']['P1']
        r_bot = res['stress_points']['P2']

        assert r_neutral['sigma_axial'] == pytest.approx(0.0, abs=1e-9)
        assert r_top['sigma_axial'] == pytest.approx(0.0, abs=1e-9)

        assert r_neutral['sigma_bending'] == pytest.approx(0.0, abs=1e-5)

        # Positive curvature compresses the top fiber and tensions the bottom.
        assert r_top['sigma_bending'] < 0
        assert r_bot['sigma_bending'] > 0
        assert abs(r_top['sigma_bending']) == pytest.approx(abs(r_bot['sigma_bending']))

        val_total_top = r_top['sigma_bending'] + r_top['sigma_axial']
        val_total_bot = r_bot['sigma_bending'] + r_bot['sigma_axial']

        assert r_top['sigma_total'] == pytest.approx(val_total_top)
        assert r_bot['sigma_total'] == pytest.approx(val_total_bot)

    def test_zero_length_error(self):
        """Reject coincident beam nodes."""
        n1 = MockNode(1, [0,0,0])
        n2 = MockNode(2, [0,0,0])
        prop = MockProp(A=1, I1=1, I2=1, J=1)
        mat = MockMat(E=1, G=1)

        with pytest.raises(ValueError, match="Zero length"):
            BarElementSolver(n1, n2, prop, mat, [0,1,0])
