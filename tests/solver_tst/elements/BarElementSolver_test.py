"""
@file BarElementSolver_test.py
@brief Unit tests for the CBAR Element Solver (Bending, Axial, Torsion). (Physics & Assembly).
@details
Validates the complete physics of the Euler-Bernoulli beam element:
1. Stiffness Matrix Assembly (Axial, Torsion, Bending XY/XZ).
2. Geometric Transformations (Arbitrary 3D rotations).
3. Stress Recovery (Internal Forces, Stresses at points, Strains, Twist Rates).
4. Error Handling (Zero length).
"""

import pytest
import numpy as np
from src.solver.elements.bar_element import BarElementSolver

# --- MOCKS (Data Objects Simulation) ---
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
    """Simulates the global-to-local index mapping."""
    def get_dof_indices(self, nid):
        # Maps Node ID to a block of 6 indices
        # Node 1 -> [0..5], Node 2 -> [6..11] (Assuming sequential numbering for tests)
        start = (nid - 1) * 6
        return list(range(start, start + 6))

# --- TESTS ---

class TestBarElementSolver:
    """
    @class TestBarElementSolver
    @brief Comprehensive test suite for CBAR element.
    """

    # --- BLOCK 1: STIFFNESS MATRIX VALIDATION ---

    def test_stiffness_matrix_bending_xy_analytical(self):
        """
        @test Validate Bending Stiffness (XY Plane).
        @brief Compares FEM matrix terms against analytical Euler-Bernoulli formulas (12EI/L^3).
        """
        L = 2.0
        E = 1000.0
        I1 = 0.001 # Izz (Flexión XY)
        
        n1 = MockNode(1, [0,0,0])
        n2 = MockNode(2, [L,0,0])
        prop = MockProp(A=0.1, I1=I1, I2=0.001, J=0.0)
        mat = MockMat(E=E, G=1000.0)
        # Orientation default Y (v=[0,1,0])
        orient = [0.0, 1.0, 0.0] 
        
        solver = BarElementSolver(n1, n2, prop, mat, orient)
        K = solver.get_stiffness_matrix()
        
        # Analytical formulas (Plane XY - DOFs: Ty, Rz)
        # Ty indices: 1 (N1), 7 (N2)
        # Rz indices: 5 (N1), 11 (N2)
        
        k_force = (12 * E * I1) / L**3
        k_moment = (4 * E * I1) / L
        k_couple = (6 * E * I1) / L**2
        
        # Assert Diagonal Terms
        assert K[1, 1] == pytest.approx(k_force)   # Ty1 - Ty1
        assert K[7, 7] == pytest.approx(k_force)   # Ty2 - Ty2
        assert K[5, 5] == pytest.approx(k_moment)  # Rz1 - Rz1
        assert K[11, 11] == pytest.approx(k_moment)# Rz2 - Rz2
        
        # Assert Off-Diagonal Coupling
        assert K[1, 5] == pytest.approx(k_couple)  # Ty1 - Rz1
        assert K[1, 7] == pytest.approx(-k_force)  # Ty1 - Ty2 (Equilibrium)


    def test_stiffness_matrix_bending_xz_sign_convention(self):
        """
        @test Validate Bending Stiffness (XZ Plane).
        @brief Checks sign convention for Plane XZ (Ry rotation vs Tz displacement).
        """
        L = 1.0
        E = 1.0
        I2 = 1.0 # Iyy
        
        n1, n2 = MockNode(1, [0,0,0]), MockNode(2, [1,0,0])
        prop = MockProp(A=0, I1=0, I2=I2, J=0)
        mat = MockMat(E=E, G=0)
        
        solver = BarElementSolver(n1, n2, prop, mat, [0,1,0])
        K = solver.get_stiffness_matrix()
        
        # DOFs Plane XZ: Tz (2, 8), Ry (4, 10)
        # Diagonal stiffness is always positive
        assert K[4, 4] == pytest.approx(4.0) # 4EI/L
        
        # Coupling term Tz-Ry usually has opposite sign to Ty-Rz 
        # because theta_y = -dw/dx
        # Expected: -6EI/L^2 = -6.0
        assert K[2, 4] == pytest.approx(-6.0)


# --- BLOCK 2: GEOMETRIC TRANSFORMATIONS ---

    def test_geometric_transformation_90_degrees(self):
        """
        @test 90 Degree Rotation.
        @brief A beam aligned with Global Y must map its Axial Stiffness to Global Ty.
        """
        n1 = MockNode(1, [0,0,0])
        n2 = MockNode(2, [0, 10.0, 0]) # Vertical Y Beam
        prop = MockProp(A=1.0, I1=0, I2=0, J=0)
        mat = MockMat(E=100.0, G=0)
        
        # Orientation vector must NOT be parallel to beam axis (Y). Use X.
        orient = [1.0, 0.0, 0.0]
        solver = BarElementSolver(n1, n2, prop, mat, orient)
        K = solver.get_stiffness_matrix()
        
        # K axial = EA/L = 100*1 / 10 = 10.0
        # In Local: Tx. In Global (rotated 90deg): Ty.
        # Global Ty indices: 1, 7
        
        assert K[1, 1] == pytest.approx(10.0)
        assert K[0, 0] == pytest.approx(0.0) # Global Tx (Shear) should be 0

    def test_arbitrary_rotation_3d_kinematics(self):
        """
        @test Arbitrary 3D Rotation & Strain Recovery.
        @brief Places bar on vector V=[1, 2, 3] and verifies that global displacement
               along that vector produces correct Axial Strain and Stress.
        """
        # 1. Setup: Bar from (0,0,0) to (1, 2, 3)
        vec = np.array([1.0, 2.0, 3.0])
        L = np.linalg.norm(vec) # sqrt(14) ~= 3.74
        unit_vec = vec / L
        
        n1 = MockNode(1, [0,0,0])
        n2 = MockNode(2, vec)
        prop = MockProp(A=0.1, I1=1e-4, I2=1e-4, J=1e-4)
        mat = MockMat(E=200e9, G=80e9)
        
        # Orientation: Global Y (solver will project it)
        orient = [0.0, 1.0, 0.0]
        
        solver = BarElementSolver(n1, n2, prop, mat, orient)
        
        # 2. Load Case: Pull Node 2 along the bar axis
        delta = 0.01 # 1 cm extension
        disp_vec = unit_vec * delta 
        
        u_global = np.zeros(12)
        # u_global indices node 2: 6(Tx), 7(Ty), 8(Tz)
        u_global[6:9] = disp_vec
        
        # 3. Solve Stresses
        mapper = MockDofMapper()
        res = solver.calculate_stresses(u_global, mapper)
        
        # 4. Validations
        expected_strain = delta / L
        expected_stress = expected_strain * 200e9
        
        # Check new kinematic outputs
        assert res['axial_strain'] == pytest.approx(expected_strain, rel=1e-5)
        assert res['axial_stress'] == pytest.approx(expected_stress, rel=1e-5)
        
        # Verify NO parasitic Bending or Torsion
        forces_n1 = res['forces_n1']
        assert abs(forces_n1['My']) < 1.0 # Negligible moments
        assert abs(res['torque']) < 1.0
    
    # --- BLOCK 3: STRESS RECOVERY ---

    def test_stress_recovery_points_torsion(self):
        """
        @test Torsion Stress at specific cross-section points.
        @brief Verifies that Shear Stress (tau) scales with radius 'r'.
        """
        n1, n2 = MockNode(1, [0,0,0]), MockNode(2, [2,0,0])
        
        # Define 3 points: Center, Mid-radius, Outer-radius
        # P0(0,0), P1(0, 0.05), P2(0.1, 0)
        points = [(0.0, 0.0), (0.0, 0.05), (0.1, 0.0)]
        prop = MockProp(A=0.1, I1=1, I2=1, J=0.5, stress_points=points)

        mat = MockMat(E=0, G=1000.0)
        
        solver = BarElementSolver(n1, n2, prop, mat, [0,1,0])
        
        # Apply pure Torsion at Node 2 (Rx = 0.1 rad)
        u_global = np.zeros(12)
        u_global[9] = 0.1 
        
        mapper = MockDofMapper()
        results = solver.calculate_stresses(u_global, mapper)
        
        # Expected Twist & Torque
        # twist = delta_theta / L = 0.1 / 2.0 = 0.05
        expected_twist = 0.05
        expected_torque = 1000.0 * 0.5 * expected_twist # G*J*twist = 25.0
        
        assert results['twist_rate'] == pytest.approx(expected_twist)
        assert results['torque'] == pytest.approx(expected_torque)
        
        # Check Stresses at Points
        stress_data = results['stress_points']
        
        # P0 (r=0) -> Tau=0
        assert stress_data['P0']['tau_xy'] == pytest.approx(0.0)
        
        # P1 (r=0.05) -> Tau = T*r/J = 25 * 0.05 / 0.5 = 2.5
        assert stress_data['P1']['tau_xy'] == pytest.approx(2.5)
        
        # P2 (r=0.1) -> Tau = T*r/J = 25 * 0.1 / 0.5 = 5.0
        assert stress_data['P2']['tau_xy'] == pytest.approx(5.0)

    def test_stress_recovery_bending_curvature(self):
        """
        @test Bending Stress distribution (Sigma = -My/I).
        @brief Verifies Tension/Compression symmetry and component breakdown (Axial vs Bending).
        """
        n1, n2 = MockNode(1, [0,0,0]), MockNode(2, [1,0,0])
        
        # Points: Neutral Axis (0,0), Top (0.1, 0), Bottom (-0.1, 0)
        points = [(0.0, 0.0), (0.1, 0.0), (-0.1, 0.0)]
        prop = MockProp(A=0.1, I1=1e-4, I2=1e-4, J=1e-4, stress_points=points)
        
        mat = MockMat(E=200e9, G=0)
        solver = BarElementSolver(n1, n2, prop, mat, [0,1,0])
        
        # Case: Constant Curvature ("Smile" shape / Sagging)
        # To get a "Smile" (Positive Curvature), Ends must rotate UPwards relative to center.
        # Node 1 (Left): Slope down (-theta)
        # Node 2 (Right): Slope up (+theta)
        theta = 0.001
        u_global = np.zeros(12)
        u_global[5] = -theta  # Rz1 (Slope down \ )
        u_global[11] = theta  # Rz2 (Slope up   / )
        
        mapper = MockDofMapper()
        res = solver.calculate_stresses(u_global, mapper)

        # Retrieve result dictionaries for each point
        r_neutral = res['stress_points']['P0'] # y=0.0
        r_top     = res['stress_points']['P1'] # y=0.1
        r_bot     = res['stress_points']['P2'] # y=-0.1
        
        # --- 1. AXIAL COMPONENT VALIDATION ---
        # Since only rotation is applied, Axial Stress must be effectively zero.
        assert r_neutral['sigma_axial'] == pytest.approx(0.0, abs=1e-9)
        assert r_top['sigma_axial']     == pytest.approx(0.0, abs=1e-9)
        
        # --- 2. BENDING COMPONENT VALIDATION ---
        # Neutral axis (y=0) must have zero bending stress.
        assert r_neutral['sigma_bending'] == pytest.approx(0.0, abs=1e-5)
        
        # Top fiber (y > 0) -> Compression (Negative) given positive curvature ("Smile" shape).
        assert r_top['sigma_bending'] < 0
        
        # Bottom fiber (y < 0) -> Tension (Positive).
        assert r_bot['sigma_bending'] > 0
        
        # Symmetry check: Magnitude at Top and Bottom should be equal.
        assert abs(r_top['sigma_bending']) == pytest.approx(abs(r_bot['sigma_bending']))

        # --- 3. TOTAL STRESS VALIDATION (Superposition) ---
        # Verify that Total Stress is exactly the sum of Bending + Axial.
        # sigma_total = sigma_axial + sigma_bending
        val_total_top = r_top['sigma_bending'] + r_top['sigma_axial']
        val_total_bot = r_bot['sigma_bending'] + r_bot['sigma_axial']

        assert r_top['sigma_total'] == pytest.approx(val_total_top)
        assert r_bot['sigma_total'] == pytest.approx(val_total_bot)
        '''
        s_neutral = res['stress_points']['P0']['sigma_total']
        s_top = res['stress_points']['P1']['sigma_total']
        s_bot = res['stress_points']['P2']['sigma_total']

        # 1. Neutral axis stress should be 0 (Pure Bending)
        assert s_neutral == pytest.approx(0.0, abs=1e-5)
        
        # 2. Top fiber (y>0) compressed (Sigma < 0) for positive curvature
        assert s_top < 0
        
        # 3. Bottom fiber (y<0) tensioned (Sigma > 0)
        assert s_bot > 0
        
        # 4. Symmetry check
        assert abs(s_top) == pytest.approx(abs(s_bot))
        '''

    # --- BLOCK 4: ERROR HANDLING ---

    def test_zero_length_error(self):
        """@test Ensure solver raises ValueError for coincident nodes."""
        n1 = MockNode(1, [0,0,0])
        n2 = MockNode(2, [0,0,0]) 
        prop = MockProp(A=1, I1=1, I2=1, J=1)
        mat = MockMat(E=1, G=1)
        
        with pytest.raises(ValueError, match="Zero length"):
            BarElementSolver(n1, n2, prop, mat, [0,1,0])



