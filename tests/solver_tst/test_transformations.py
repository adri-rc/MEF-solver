import pytest
import numpy as np
from src.solver.transformations import CoordinateTransformation

"""
@file test_transformations.py
@brief Verification of geometric transformations and linear algebra logic.
@details
This suite tests the robustness of coordinate system transformations under various
spatial orientations.

Test Coverage:
    - Vector Orthogonality: Ensures local axes (x, y, z) are unitary and perpendicular.
    - SO(3) Properties: Validates that Rotation matrices have Determinant=1.0 and
      that their Inverse equals their Transpose (R^-1 = R^T).
    - Right-Hand Rule: Checks consistency of cross products.
    - Matrix Assembly: Verifies the correct construction of the 12x12 transformation
      matrix (block-diagonal structure).
    - Edge Cases: Tests stability with random vectors and parallel alignment cases.
"""

class TestCoordinateTransformation:
    """
    @class TestCoordinateTransformation
    @brief Check linear algebra logic for local coordinate systems.
    """

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
        """
        @test Orthogonality in multiple directions.
        @brief Check that for various orientations, the local corrd system remains orthonormal.
        """
        n1 = np.array(n1_coords)
        n2 = np.array(n2_coords)
        
        vx, vy, vz = CoordinateTransformation.compute_local_axes(n1, n2)            
        
        norms = [np.linalg.norm(v) for v in [vx, vy, vz]]
        np.testing.assert_allclose(norms, 1.0, err_msg="Vectors are not unitary")

        # Orthogonality across components (dot product)            
        assert np.dot(vx, vy) == pytest.approx(0.0, abs=1e-9) # vx . vy            
        assert np.dot(vx, vz) == pytest.approx(0.0, abs=1e-9) # vx . vz            
        assert np.dot(vy, vz) == pytest.approx(0.0, abs=1e-9) # vy . vz
    
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
        """
        @test Properties of the Rotation Matrix.
        @brief
        @details
        -> Det(R)=1.0
        -> R^-1=R_transpose
        -> Right-hand rule (Cross Product)
        """
        n1 = np.array(n1_coords)
        n2 = np.array(n2_coords)
        
        vx, vy, vz = CoordinateTransformation.compute_local_axes(n1, n2)
        R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)
        
        #Det(R)=1.0
        det = np.linalg.det(R) 
        assert det == pytest.approx(1.0), f"Det(R) is '{det}', value must be 1.0 (pure rotation)"

        #R^-1=R_transpose
        R_transpose = R.T
        R_inverse = np.linalg.inv(R)
        np.testing.assert_allclose(R_inverse,R_transpose,atol=1e-9,
                                   err_msg="Matrix is not orthogonal: Inverse != Transpose")
        
        # Right-hand rule (Cross Product)        
        expected_vz = np.cross(vx, vy) # Z must be exactly X cross Y
        np.testing.assert_allclose(vz, expected_vz, atol=1e-9,
                                   err_msg="Does not comply with Right-hand rule.")
    
    def test_random_orientations_stress(self):
        """
        @test Stress Test Aleatorio.
        @brief Genera 100 barras con coordenadas aleatorias para detectar inestabilidad numérica.
        """
        np.random.seed(42)        
        for _ in range(100):
            
            n1 = np.random.uniform(-100, 100, 3) # Random point
            diff = np.random.uniform(-100, 100, 3) #random displacement vector
            if np.linalg.norm(diff) < 1e-3: continue # Ensure they are not coincident (Length > 0)
            n2 = n1 + diff

            vx, vy, vz = CoordinateTransformation.compute_local_axes(n1, n2)
            R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)
            
            assert np.linalg.det(R) == pytest.approx(1.0), f"Det(R) is '{np.linalg.det(R)}', value must be 1.0 (pure rotation)"
        
            identity_check = R @ R.T
            np.testing.assert_allclose(
                identity_check,np.eye(3),
                atol=1e-9, 
                err_msg=f"Orthogonality failure with n1={n1} n2={n2}"
            )
    
    def test_zero_length_vector(self):
        """
        @test Validation of 'zero length' messages.
        @brief Raise an error if element nodes n1 and n2 are coincident.
        """
        n1 = np.array([10.0, 10.0, 10.0])
        n2 = n1
        
        with pytest.raises(ValueError) as excinfo:
            CoordinateTransformation.compute_local_axes(n1, n2)
        
        assert "Zero length" in str(excinfo.value) or "zero length" in str(excinfo.value)

    def test_full_matrix_structure_standard_rod(self):
        """
        @test T-matrix structure for a standard bar (2 nodes, 6 DOFs).
        @brief Check that R is copied 4 times on the diagonal (Translation1, Rotation1, Translation2, Rotation2).
        """
        # Easy-to-identify fictitious “R” matrix
        R_dummy = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0]])
        
        
        T = CoordinateTransformation.build_full_transformation_matrix(
            R_3x3=R_dummy, 
            nodes_count=2, 
            dofs_per_node=6
        )

        assert T.shape == (12, 12) # check shape (2 nodes * 6 dofs= 12x12)
        
        # Check the 4 blocks of 3x3 on the diagonal
        np.testing.assert_array_equal(T[0:3, 0:3], R_dummy, err_msg="Block node 1 failure - Translation")        
        np.testing.assert_array_equal(T[3:6, 3:6], R_dummy, err_msg="Block node 1 failure - Rotation")        
        np.testing.assert_array_equal(T[6:9, 6:9], R_dummy, err_msg="Block node 2 failure - Translation")        
        np.testing.assert_array_equal(T[9:12, 9:12], R_dummy, err_msg="Block node 2 failure - Rotation")
    
    def test_full_matrix_zeros_off_diagonal(self):
        """
        @test Zeros outside the T matrix diagonal.
        @brief Checks that there is no coupling between nodes in the transformation matrix.
        """
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

    # --- BLOCK 2: TRIANGULAR ELEMENT AXES ---

    def test_compute_local_axes_triangle_standard(self):
        """
        @test Local Axes for Triangle (Standard 3-4-5).
        @brief Validates X along N1-N2, Z normal to plane.
        """
        n1 = np.array([0.0, 0.0, 0.0])
        n2 = np.array([3.0, 0.0, 0.0])
        n3 = np.array([0.0, 4.0, 0.0])
        
        vx, vy, vz = CoordinateTransformation.compute_local_axes_triangle(n1, n2, n3)
        
        # X: Along N1->N2 = (1,0,0)
        np.testing.assert_allclose(vx, [1,0,0], atol=1e-9)
        
        # Z: Normal to plane XY = (0,0,1)
        np.testing.assert_allclose(vz, [0,0,1], atol=1e-9)
        
        # Y: Cross(Z, X) = (0,1,0)
        np.testing.assert_allclose(vy, [0,1,0], atol=1e-9)

    def test_compute_local_axes_triangle_orthogonality_random(self):
        """
        @test Orthogonality for Random Triangle.
        """
        np.random.seed(99)
        for _ in range(20):
            n1 = np.random.rand(3)
            n2 = np.random.rand(3)
            n3 = np.random.rand(3)
            
            # Avoid collinearity crudely
            if np.linalg.norm(np.cross(n2-n1, n3-n1)) < 1e-4: continue
            
            vx, vy, vz = CoordinateTransformation.compute_local_axes_triangle(n1, n2, n3)
            
            # Checks norms
            assert np.linalg.norm(vx) == pytest.approx(1.0)
            assert np.linalg.norm(vy) == pytest.approx(1.0)
            assert np.linalg.norm(vz) == pytest.approx(1.0)
            
            # Check orthogonality
            assert np.dot(vx, vy) == pytest.approx(0.0, abs=1e-9)
            assert np.dot(vx, vz) == pytest.approx(0.0, abs=1e-9)
            assert np.dot(vy, vz) == pytest.approx(0.0, abs=1e-9)

    def test_compute_local_axes_triangle_collinear_error(self):
        """@test Error on collinear nodes."""
        n1 = np.array([0,0,0])
        n2 = np.array([1,0,0])
        n3 = np.array([2,0,0]) # Collinear along X
        
        with pytest.raises(ValueError, match="collinear"):
            CoordinateTransformation.compute_local_axes_triangle(n1, n2, n3)


