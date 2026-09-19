"""Tests for element coordinate systems and transformation matrices."""

import numpy as np
import pytest

from src.solver.transformations import CoordinateTransformation


class TestCoordinateTransformation:
    @pytest.mark.parametrize("n1_coords, n2_coords", [
        ([0.0,0.0,0.0], [1.0,1.0,1.0]),         # Positive 3D diagonal
        ([0.0,0.0,0.0], [-1.0,-1.0,-1.0]),      # Negative 3D diagonal
        ([10.0,0.0,30.0], [50.0,0.0,30.0]),     # X displacement
        ([0.0,10.0,-10.0], [0.0,20.0,-10.0]),   # Y displacement (critical case for the auxiliary vector)
        ([5.0,5.0,5.0], [5.0,6.0,5.0]),         # Y displacement (critical case for the auxiliary vector)
        ([-10.0,20.0,30.0], [-10.0,20.0,40.0]), # Z displacement
        ([0.0,0.0,0.0], [1.0, 0.1, 0.0])           # Almost parallel to X (numerical stability test)
    ])

    def test_orthogonality_parametrized(self, n1_coords, n2_coords):
        n1 = np.array(n1_coords)
        n2 = np.array(n2_coords)

        vx, vy, vz = CoordinateTransformation.compute_local_axes(n1, n2)

        norms = [np.linalg.norm(v) for v in [vx, vy, vz]]
        np.testing.assert_allclose(norms, 1.0, err_msg="Vectors are not unitary")

        assert np.dot(vx, vy) == pytest.approx(0.0, abs=1e-9)
        assert np.dot(vx, vz) == pytest.approx(0.0, abs=1e-9)
        assert np.dot(vy, vz) == pytest.approx(0.0, abs=1e-9)

    @pytest.mark.parametrize("n1_coords, n2_coords", [
        ([0.0,0.0,0.0], [1.0,1.0,1.0]),         # Positive 3D diagonal
        ([0.0,0.0,0.0], [-1.0,-1.0,-1.0]),      # Negative 3D diagonal
        ([10.0,0.0,30.0], [50.0,0.0,30.0]),     # X displacement
        ([0.0,10.0,-10.0], [0.0,20.0,-10.0]),   # Y displacement (critical case for the auxiliary vector)
        ([5.0,5.0,5.0], [5.0,6.0,5.0]),         # Y displacement (critical case for the auxiliary vector)
        ([-10.0,20.0,30.0], [-10.0,20.0,40.0]), # Z displacement
        ([0.0,0.0,0.0], [1.0, 0.1, 0.0])           # Almost parallel to X (numerical stability test)
    ])

    def test_rotation_matrix_properties(self, n1_coords, n2_coords):
        """The local basis must be a proper orthogonal rotation matrix."""
        n1 = np.array(n1_coords)
        n2 = np.array(n2_coords)

        vx, vy, vz = CoordinateTransformation.compute_local_axes(n1, n2)
        R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)

        det = np.linalg.det(R)
        assert det == pytest.approx(1.0), f"Det(R) is '{det}', value must be 1.0 (pure rotation)"

        R_transpose = R.T
        R_inverse = np.linalg.inv(R)
        np.testing.assert_allclose(R_inverse, R_transpose, atol=1e-9,
                                   err_msg="Matrix is not orthogonal: Inverse != Transpose")

        expected_vz = np.cross(vx, vy)
        np.testing.assert_allclose(vz, expected_vz, atol=1e-9,
                                   err_msg="Does not comply with Right-hand rule.")

    def test_random_orientations_stress(self):
        """Exercise proper rotations over 100 deterministic random bars."""
        np.random.seed(42)
        for _ in range(100):
            n1 = np.random.uniform(-100, 100, 3)
            diff = np.random.uniform(-100, 100, 3)
            if np.linalg.norm(diff) < 1e-3:
                continue
            n2 = n1 + diff

            vx, vy, vz = CoordinateTransformation.compute_local_axes(n1, n2)
            R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)

            assert np.linalg.det(R) == pytest.approx(1.0), f"Det(R) is '{np.linalg.det(R)}', value must be 1.0 (pure rotation)"

            identity_check = R @ R.T
            np.testing.assert_allclose(
                identity_check, np.eye(3),
                atol=1e-9,
                err_msg=f"Orthogonality failure with n1={n1} n2={n2}"
            )

    def test_zero_length_vector(self):
        """Reject coincident nodes with a diagnostic message."""
        n1 = np.array([10.0, 10.0, 10.0])
        n2 = n1

        with pytest.raises(ValueError) as excinfo:
            CoordinateTransformation.compute_local_axes(n1, n2)

        assert "Zero length" in str(excinfo.value) or "zero length" in str(excinfo.value)

    def test_full_matrix_structure_standard_rod(self):
        """Copy R into the translational and rotational blocks of both nodes."""
        R_dummy = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0]])

        T = CoordinateTransformation.build_full_transformation_matrix(
            R_3x3=R_dummy,
            nodes_count=2,
            dofs_per_node=6
        )

        assert T.shape == (12, 12)

        np.testing.assert_array_equal(T[0:3, 0:3], R_dummy, err_msg="Block node 1 failure - Translation")
        np.testing.assert_array_equal(T[3:6, 3:6], R_dummy, err_msg="Block node 1 failure - Rotation")
        np.testing.assert_array_equal(T[6:9, 6:9], R_dummy, err_msg="Block node 2 failure - Translation")
        np.testing.assert_array_equal(T[9:12, 9:12], R_dummy, err_msg="Block node 2 failure - Rotation")

    def test_full_matrix_zeros_off_diagonal(self):
        """Keep every off-diagonal transformation block zero."""
        R_dummy = np.ones((3,3))
        T = CoordinateTransformation.build_full_transformation_matrix(R_dummy, 2, 6)

        block_off_1 = T[0:3, 6:9]
        block_off_2 = T[3:6, 0:3]
        block_off_3 = T[6:9, 9:]
        block_off_4 = T[9:, 0:3]
        block_off_5 = T[9:, 6:9]
        expected_zeros = np.zeros((3,3))

        np.testing.assert_array_equal(block_off_1, expected_zeros, err_msg="Transformation T matrix is not 0 block diagonal.")
        np.testing.assert_array_equal(block_off_2, expected_zeros, err_msg="Transformation T matrix is not 0 block diagonal.")
        np.testing.assert_array_equal(block_off_3, expected_zeros, err_msg="Transformation T matrix is not 0 block diagonal.")
        np.testing.assert_array_equal(block_off_4, expected_zeros, err_msg="Transformation T matrix is not 0 block diagonal.")
        np.testing.assert_array_equal(block_off_4, expected_zeros, err_msg="Transformation T matrix is not 0 block diagonal.")
        np.testing.assert_array_equal(block_off_5, expected_zeros, err_msg="Transformation T matrix is not 0 block diagonal.")

    def test_compute_local_axes_triangle_standard(self):
        """For an XY triangle, X follows N1-N2 and Z is the plane normal."""
        n1 = np.array([0.0, 0.0, 0.0])
        n2 = np.array([3.0, 0.0, 0.0])
        n3 = np.array([0.0, 4.0, 0.0])
        
        vx, vy, vz = CoordinateTransformation.compute_local_axes_triangle(n1, n2, n3)
        
        np.testing.assert_allclose(vx, [1,0,0], atol=1e-9)
        np.testing.assert_allclose(vz, [0,0,1], atol=1e-9)
        np.testing.assert_allclose(vy, [0,1,0], atol=1e-9)

    def test_compute_local_axes_triangle_orthogonality_random(self):
        """Build orthonormal bases for deterministic random triangles."""
        np.random.seed(99)
        for _ in range(20):
            n1 = np.random.rand(3)
            n2 = np.random.rand(3)
            n3 = np.random.rand(3)
            
            if np.linalg.norm(np.cross(n2-n1, n3-n1)) < 1e-4:
                continue

            vx, vy, vz = CoordinateTransformation.compute_local_axes_triangle(n1, n2, n3)

            assert np.linalg.norm(vx) == pytest.approx(1.0)
            assert np.linalg.norm(vy) == pytest.approx(1.0)
            assert np.linalg.norm(vz) == pytest.approx(1.0)
            
            assert np.dot(vx, vy) == pytest.approx(0.0, abs=1e-9)
            assert np.dot(vx, vz) == pytest.approx(0.0, abs=1e-9)
            assert np.dot(vy, vz) == pytest.approx(0.0, abs=1e-9)

    def test_compute_local_axes_triangle_collinear_error(self):
        """Reject collinear triangle nodes."""
        n1 = np.array([0,0,0])
        n2 = np.array([1,0,0])
        n3 = np.array([2,0,0])

        with pytest.raises(ValueError, match="collinear"):
            CoordinateTransformation.compute_local_axes_triangle(n1, n2, n3)
    
    def test_transformation_matrix_energy_conservation(self):
        """An orthogonal T preserves energy because T.T @ T equals identity."""
        n1 = np.array([1.0, 1.0, 1.0])
        n2 = np.array([4.0, 1.0, 1.0])
        n3 = np.array([2.0, 4.0, 1.0])

        vx, vy, vz = CoordinateTransformation.compute_local_axes_triangle(n1, n2, n3)
        
        R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)
        T = CoordinateTransformation.build_full_transformation_matrix(R, nodes_count=3, dofs_per_node=6)

        assert T.shape == (18, 18)

        I_calc = T.T @ T
        I_expected = np.eye(18)

        np.testing.assert_allclose(
            I_calc, 
            I_expected, 
            atol=1e-12, 
            err_msg="Transformation matrix T is not orthogonal. Energy will leak!"
        )

    def test_right_hand_rule(self):
        """The local Z axis follows the right-hand rule."""
        vx, vy, vz = CoordinateTransformation.compute_local_axes_triangle(
            np.array([0,0,0]), np.array([1,0,0]), np.array([0,1,0])
        )
        
        expected_z = np.array([0.0, 0.0, 1.0])
        np.testing.assert_allclose(vz, expected_z, atol=1e-12)
