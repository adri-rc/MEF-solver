"""
@file test_rod_element.py
@brief Unit tests for the Rod Element Solver (Physics & Assembly).
@details
Validates:
1. Axial Stiffness calculation (EA/L).
2. Torsional Stiffness calculation (GJ/L).
3. Zero Coupling verification (Axial vs Torsion vs Bending).
4. Coordinate rotation logic (45 degree projection check).
"""
import pytest
import numpy as np
from src.solver.elements.rod_element_v2 import RodElementSolverV2

# --- MOCKS (Data Objects Simulation) ---

class MockNode:
    def __init__(self, nid, coords):
        self.nid = nid
        self._coords = np.array(coords, dtype=float)
    def get_coords(self):
        return self._coords

class MockProp:
    def __init__(self, A, J):
        self.A = A
        self.J = J

class MockMat:
    def __init__(self, E, G):
        self.E = E
        self.G = G

# --- TESTS
class TestRodElementSolver:

    def test_stiffness_axial_x_aligned(self):
        """
        @test Axial Stiffness (Aligned with Global X).
        @brief Validates K = (EA/L) * [1 -1; -1 1] for Tx DOFs, , and ensures no axial coupling.
        """
        n1 = MockNode(1, [0.0, 0.0, 0.0])
        n2 = MockNode(2, [2.0, 0.0, 0.0])
        
        # Physics: E=1000, A=0.5 -> EA/L = 500 / 2 = 250 (expected)
        prop = MockProp(A=0.5, J=0.0)
        mat = MockMat(E=1000.0, G=0.0)
        expected_k = 250.0

        solver = RodElementSolverV2(n1, n2, prop, mat, n_gaussp=1)
        K = solver.get_stiffness_matrix()
        
        # Verify Axial Indices (Tx for Node 1 and Node 2 -> Indices 0 and 6)
        assert K[0, 0] == pytest.approx(expected_k)
        assert K[0, 6] == pytest.approx(-expected_k)
        assert K[6, 0] == pytest.approx(-expected_k)
        assert K[6, 6] == pytest.approx(expected_k)
        
        # Verify decoupling.
        axial_indices = [0, 6]
        
        for row in axial_indices:
            for col in range(12):
                if col not in axial_indices:
                    assert K[row, col] == pytest.approx(0.0, abs=1e-9), \
                        f"Unexpected coupling between Axial (idx {row}) and DOF {col}"
    
    def test_stiffness_torsional_x_aligned(self):
        """
        @test Torsional Stiffness.
        @brief Validates K_tor = (GJ/L) * [1 -1; -1 1] for Rx DOFs, and ensures no torsion coupling.
        """
        n1 = MockNode(1, [0.0, 0.0, 0.0])
        n2 = MockNode(2, [1.0, 0.0, 0.0]) # L = 1.0
        
        # Physics: G=800, J=0.1 -> GJ/L = 80 / 1 = 80
        prop = MockProp(A=0.0, J=0.1)
        mat = MockMat(E=0.0, G=800.0)
        expected_k = 80.0

        solver = RodElementSolverV2(n1, n2, prop, mat)
        K = solver.get_stiffness_matrix()

        # Verify Torsion Indices (Rx for Node 1 and 2 -> Indices 3 and 9)
        assert K[3, 3] == pytest.approx(expected_k)
        assert K[3, 9] == pytest.approx(-expected_k)
        assert K[9, 9] == pytest.approx(expected_k)

        # Verify decoupling.
        torsion_indices = [3, 9]

        for row in torsion_indices:
            for col in range(12):
                if col not in torsion_indices:
                    assert K[row, col] == pytest.approx(0.0, abs=1e-9), \
                        f"Unexpected coupling between Torsion (idx {row}) and DOF {col}"
    
    def test_transformation_45_degrees(self):
        """
        @test Transformation Matrix Correctness.
        @brief Verifies that a bar at 45 degrees distributes stiffness correctly.
        @details
        If a bar is at 45 deg in XY plane:
        - cos(a) = sin(a) = 1/sqrt(2) approx 0.707
        - Local Axial Stiffness 'k' projects to Global as:
            K_xx = k * cos^2(a) = k * 0.5
            K_yy = k * sin^2(a) = k * 0.5
            K_xy = k * cos(a)*sin(a) = k * 0.5
        """
        # Setup: 45 degree bar
        # L = sqrt(1^2 + 1^2) = sqrt(2) approx 1.4142
        n1 = MockNode(1, [0.0, 0.0, 0.0])
        n2 = MockNode(2, [1.0, 1.0, 0.0]) 
        
        L = np.sqrt(2.0)
        E_val = 100.0
        A_val = 1.0
        k_local = (E_val * A_val) / L  # ~70.71
        
        prop = MockProp(A=A_val, J=0.0)
        mat = MockMat(E=E_val, G=0.0)
        
        solver = RodElementSolverV2(n1, n2, prop, mat, n_gaussp=1)
        K = solver.get_stiffness_matrix()
        
        # Expected projection factor: cos(45)^2 = 0.5
        k_projected = k_local * 0.5
        
        # 1. Verify Diagonal terms (Stiffness in X and Y should be equal)        
        assert K[0, 0] == pytest.approx(k_projected) # K_xx (0,0)
        assert K[1, 1] == pytest.approx(k_projected) # K_yy (1,1)
        
        # 2. Verify Coupling term (Force in X causes displacement in Y)        
        assert K[0, 1] == pytest.approx(k_projected) # K_xy (0,1)
        
        # 3. Verify Node 2 interaction (Negative values)        
        assert K[0, 6] == pytest.approx(-k_projected) # K_xx_n1n2 (0, 6)        
        assert K[0, 7] == pytest.approx(-k_projected) # K_xy_n1n2 (0, 7)
    

    def test_zero_length_error(self):
        """
        @test Error Handling for Geometry.
        @brief Verifies that the solver raises a ValueError with a specific message 
               when nodes are coincident.
        """
        n1 = MockNode(1, [1.0, 1.0, 1.0])
        n2 = MockNode(2, [1.0, 1.0, 1.0])
        prop = MockProp(A=1.0, J=1.0)
        mat = MockMat(E=1.0, G=1.0)
        
        with pytest.raises(ValueError)as excinfo:
            RodElementSolverV2(n1, n2, prop, mat)
        error_msg = str(excinfo.value)
        assert "Zero length" in error_msg
        assert "Nodes 1-2" in error_msg or "id nodes" in error_msg






