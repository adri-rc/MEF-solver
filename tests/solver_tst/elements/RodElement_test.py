"""
@file RodElement_test.py
@brief Unit tests for the Rod Element Solver (Physics & Assembly).
@details
Validates:
1. Axial Stiffness calculation (EA/L).
2. Torsional Stiffness calculation (GJ/L).
3. Zero Coupling verification (Axial vs Torsion vs Bending).
4. Coordinate rotation logic (45 degree projection check).
5. Axial strain and twist rate values.
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
    def __init__(self, A, J, C=0.0):
        self.A = A
        self.J = J
        self.C = C

class MockMat:
    def __init__(self, E, G):
        self.E = E
        self.G = G

class MockDofMapper:
    """Simulates the DofMapper behavior for stress calculation tests."""
    def get_dof_indices(self, node_id):
        # Logic: Node 1 -> [0..5], Node 2 -> [6..11]
        start = (node_id - 1) * 6
        return list(range(start, start + 6))

# --- TESTS ---
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

    def test_calculate_stresses_axial_tension(self):
        """
        @test Stress Recovery (Axial Tension).
        @brief Verifies stress calculation for a bar pulled in X.
        """
        # Setup: L=2.0
        n1 = MockNode(1, [0., 0., 0.])
        n2 = MockNode(2, [2., 0., 0.])
        prop = MockProp(A=10.0, J=1.0)
        mat = MockMat(E=2000.0, G=1000.0)
        
        solver = RodElementSolverV2(n1, n2, prop, mat)
        mapper = MockDofMapper()
        
        # Create a Global Displacement Vector (size 12 for 2 nodes)
        u_global = np.zeros(12)
        
        # Apply displacement: Node 2 moves +0.1 in X (Global Index 6)
        # Node 1 is fixed (0.0).
        u_global[6] = 0.1
        
        # Expected Physics:
        # Strain = delta_L / L = 0.1 / 2.0 = 0.05
        # Stress = E * Strain = 2000 * 0.05 = 100.0
        # Force = Stress * A = 100.0 * 10.0 = 1000.0
        
        results = solver.calculate_stresses(u_global, mapper)
        
        assert results['axial_stress'] == pytest.approx(100.0)
        assert results['axial_force'] == pytest.approx(1000.0)
        # Verify no torsion was induced
        assert results['torque'] == pytest.approx(0.0)

    def test_calculate_stresses_rotated_rod(self):
        """
        @test Stress Recovery with Coordinate Transformation.
        @brief Verifies that global displacements are correctly rotated to local
               before calculating stress.
        @details
        Bar at 90 degrees (Along Y axis).
        Global displacement in Y should cause AXIAL stress.
        Global displacement in X should cause NO stress (shear/bending ignored).
        """
        # Setup: Bar along Y axis, L=1.0
        n1 = MockNode(1, [0., 0., 0.])
        n2 = MockNode(2, [0., 1., 0.])
        prop = MockProp(A=1.0, J=1.0)
        mat = MockMat(E=1000.0, G=1000.0)
        
        solver = RodElementSolverV2(n1, n2, prop, mat)
        mapper = MockDofMapper()
        
        u_global = np.zeros(12)
        
        # Apply Global Displacement in Y at Node 2 (Index 7) -> +0.1
        # In Local System (where X_local aligns with Global Y), this is a local extension.
        u_global[7] = 0.1 
        
        results = solver.calculate_stresses(u_global, mapper)
        
        # Strain = 0.1 / 1.0 = 0.1
        # Stress = 1000 * 0.1 = 100.0
        assert results['axial_stress'] == pytest.approx(100.0)
        
        # --- Case 2: Apply Global X displacement (Index 6) ---
        # This is perpendicular to the bar. Should cause NO axial stress in a Rod.
        u_global_shear = np.zeros(12)
        u_global_shear[6] = 0.1
        
        results_shear = solver.calculate_stresses(u_global_shear, mapper)
        assert results_shear['axial_stress'] == pytest.approx(0.0, abs=1e-9)

    def test_calculate_stresses_rotated_rod_torsion(self):
        """
        @test Stress Recovery with Coordinate Transformation (Rotations).
        @brief Verifies that global rotations are correctly transformed to local twist.
        @details
        Bar at 90 degrees (Along Global Y axis).
        - Global Rotation about Y (Ry) -> Aligns with Local X -> CAUSES TORSION.
        - Global Rotation about X (Rx) -> Perpendicular to Local X -> NO STRESS (Bending ignored).
        """
        # Setup: Bar along Y axis, L=1.0
        n1 = MockNode(1, [0., 0., 0.])
        n2 = MockNode(2, [0., 1., 0.])
        
        # Properties for Torsion
        # J=2.0, C=0.5 (Distance to extreme fiber)
        prop = MockProp(A=10.0, J=2.0, C=0.5)
        mat = MockMat(E=1000.0, G=1000.0) # G=1000
        
        solver = RodElementSolverV2(n1, n2, prop, mat)
        mapper = MockDofMapper()
        
        # --- Case 1: Apply Global Rotation in Y (Ry) ---
        # Node 2 indices: [6, 7, 8, 9(Rx), 10(Ry), 11(Rz)]
        # We apply +0.1 rad to Ry (Index 10)
        u_global_tor = np.zeros(12)
        u_global_tor[10] = 0.1
        
        results = solver.calculate_stresses(u_global_tor, mapper)
        
        # Physics Check:
        # The bar is vertical, so Global Ry == Local Twist (Rx_local)
        # Twist Rate = 0.1 / 1.0 = 0.1
        # Torque = G * J * Rate = 1000 * 2.0 * 0.1 = 200.0
        # Shear Stress = (T * C) / J = (200 * 0.5) / 2.0 = 50.0
        
        assert results['torque'] == pytest.approx(200.0)
        assert results['torsional_stress'] == pytest.approx(50.0)
        
        # --- Case 2: Apply Global Rotation in X (Rx) ---
        # We apply +0.1 rad to Rx (Index 9)
        # In Local system (Bar along Y), Global X is a transversal axis.
        # Ideally, this is bending. Rod element has NO bending stiffness.
        u_global_bend = np.zeros(12)
        u_global_bend[9] = 0.1
        
        results_bend = solver.calculate_stresses(u_global_bend, mapper)
        
        # Should be zero torque
        assert results_bend['torque'] == pytest.approx(0.0, abs=1e-9)
        assert results_bend['torsional_stress'] == pytest.approx(0.0, abs=1e-9)

    def test_stress_recovery_axial_and_strain(self):
        """
        @test  Axial Strain and Stress validation.
        @brief Extend the bar by 0.01m and check epsilon, sigma, and force.
        """
        # Setup: 1m rod along X
        n1, n2 = MockNode(1, [0,0,0]), MockNode(2, [1,0,0])
        prop = MockProp(A=0.1, J=0.0) 
        mat = MockMat(E=200e9, G=0)
        
        solver = RodElementSolverV2(n1, n2, prop, mat)
        mapper = MockDofMapper()

        # Case: Displacement u2 = 0.01 (Traction)
        # Global indices: u1=0..5, u2=6..11. Tx is at 0 and 6.
        u_global = np.zeros(12)
        u_global[6] = 0.01 # Tx in node 2
        
        results = solver.calculate_stresses(u_global, mapper)
        
        expected_strain = 0.01 / 1.0  # deltaL / L
        expected_stress = expected_strain * 200e9 # Hooke Law
        
        assert results['axial_strain'] == pytest.approx(expected_strain)
        assert results['axial_stress'] == pytest.approx(expected_stress)
        assert results['axial_force'] == pytest.approx(expected_stress * 0.1)

    def test_stress_recovery_torsion_and_twist(self):
            """
            @test  Twist Rate and Torque validation.
            @brief Rotate node 2 by 0.1 rad and check the twist_rate.
            """
            # Setup: Vertical rod (to test that it rotates with local axes)
            # L = 2.0m
            n1, n2 = MockNode(1, [0,0,0]), MockNode(2, [0,2.0,0])
            prop = MockProp(A=0.1, J=0.5)
            mat = MockMat(E=0, G=80e9)
            
            solver = RodElementSolverV2(n1, n2, prop, mat)
            mapper = MockDofMapper()
            
            # Case: Global rotation in Y (Ry) on a vertical bar (Y axis).
            # Global Ry corresponds to local Rx (torsion).
            # Index Ry node 2 is 10 (6+4)
            u_global = np.zeros(12)
            u_global[10] = 0.1 # 0.1 radianes
            
            results = solver.calculate_stresses(u_global, mapper)
            
            expected_twist = 0.1 / 2.0 # theta / L = 0.05 rad/m
            expected_torque = 80e9 * 0.5 * expected_twist # G * J * twist
            
            assert results['twist_rate'] == pytest.approx(expected_twist)
            assert results['torque'] == pytest.approx(expected_torque)
            # Check that it does not ‘invent’ axial force.
            assert results['axial_force'] == pytest.approx(0.0)