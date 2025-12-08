"""
@file test_linear_static_solver.py
@brief Integration tests ensuring the Solver correctly orchestrates the analysis.
@details
Covered Scenarios:
1. End-to-End Solution: Verifying Hooke's Law (PL/AE) on a simple rod model.
2. Enforced Displacements: Verifying non-zero SPCs.
3. Error Handling: Verifying that the solver rejects unsupported entities (Fail Fast).
"""
import pytest
import numpy as np
from src.solver.linear_static_solver import LinearStaticSolver
from src.preprocessing.data_model import (
    FEMModel, Node, IsotropicMaterial, RodProperty, RodElement, 
    SPCConstraint, SinglePointLoad, SinglePointMoment
)

class TestLinearStaticSolverIntegration:
    """
    @class TestLinearStaticSolverIntegration
    @brief Integration test suite for the Linear Static Solver.
    """

    # --- HAPPY PATH TESTS (Physics Validation) ---

    def test_solve_simple_tension_rod(self):
        """
        @test Validate Hooke's Law (Tension).
        @brief Solves a 1-element Rod model and checks displacement against analytical formula.
        @details
        Physical Setup:
            - Clamped at Node 1 (x=0).
            - Tensile Force Fx=1000N at Node 2 (x=2).
            - Properties: L=2.0m, A=1.0m^2, E=1000Pa.
        
        Analytical Solution:
            delta = (P * L) / (A * E) 
            delta = (1000 * 2) / (1 * 1000) = 2.0 meters.
        """
        # 1. Setup Model
        model = FEMModel()
        
        # Nodes
        model.nodes[1] = Node(nid=1, coord_id=0, x1=0., x2=0., x3=0., cd_id=0, ps="", seid=0)
        model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0., 0., 0.
        
        model.nodes[2] = Node(nid=2, coord_id=0, x1=2., x2=0., x3=0., cd_id=0, ps="", seid=0)
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 2., 0., 0.

        # Properties
        model.materials[1] = IsotropicMaterial(mid=1, E=1000.0, nu=0.3)
        model.properties[1] = RodProperty(pid=1, mat_id=1, A=1.0, J=1.0, C=0.0, NSM=0.0)

        # Element
        model.elements[10] = RodElement(eid=10, prop_id=1, nodes=[1, 2])

        # Constraints (Clamp Node 1)
        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))
        # Stability (Roller Node 2)
        model.constraints.append(SPCConstraint(sid=1, node_id=2, components="23456", value=0.0))

        # Load
        model.loads.append(SinglePointLoad(sid=1, appnode_id=2, module=1000.0, vector=[1.0, 0.0, 0.0]))

        # 2. Run Solver
        solver = LinearStaticSolver(model)
        solver.assemble()
        displacements = solver.solve()
        
        # 3. Validate
        dofs_n2 = solver.mapper.get_dof_indices(2)
        u_x_n2 = displacements[dofs_n2[0]] # Index for Tx
        
        expected_u = 2.0 # metres
        assert u_x_n2 == pytest.approx(expected_u, rel=1e-5)

    def test_solve_enforced_displacement(self):
        """
        @test Enforced Displacement (Non-zero SPC).
        @brief Verifies that the solver correctly handles non-zero boundary conditions.
        @details
        Model:
            - Node 1: Fixed (u=0).
            - Node 2: Enforced displacement in X = +0.01 m.
            - No external forces.
        """
        model = FEMModel()
        # REMIND FORCING THE GLOBAL COORDINATE DEFINITION x, y, z of each node. Default = 0., 0., 0.
        model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
        model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0.0, 0.0, 0.0
        model.nodes[2] = Node(2, 0, 1.0, 0., 0., 0, "", 0) # L=1.0m
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 1.0, 0.0, 0.0

        model.materials[1] = IsotropicMaterial(mid=1, E=1000.0, nu=0.3)
        model.properties[1] = RodProperty(pid=1, mat_id=1, A=1.0, J=1.0, C=0.0, NSM=0.0)
        model.elements[1] = RodElement(eid=1, prop_id=1, nodes=[1, 2])

        # SPC Node 1: Fully Fixed
        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))
        
        # SPC Node 2: Enforced Displacement X=0.01
        model.constraints.append(SPCConstraint(sid=1, node_id=2, components="1", value=0.01))
        # Fix other DOFs for stability, leave Dof 4 free to avoid 'Model fully constrain error.'
        # DOFs Ty, Tz cannot be leave free because there is no stiffnes in that direction and therefore solver would try to divide by 0
        model.constraints.append(SPCConstraint(sid=1, node_id=2, components="2356", value=0.0)) # components="2356"

        solver = LinearStaticSolver(model)
        solver.assemble()
        displacements = solver.solve()

        # Check Node 2 X displacement
        dofs_n2 = solver.mapper.get_dof_indices(2)
        u_x_n2 = displacements[dofs_n2[0]]

        assert u_x_n2 == pytest.approx(0.01, abs=1e-9)
    


    # --- ERROR HANDLING TESTS (Fail Fast) ---

    def test_error_unsupported_element(self):
        """@test Fail on Unsupported Element (e.g., QUAD)."""
        from src.preprocessing.data_model import QUADElement
        
        model = FEMModel()
        model.nodes[1] = Node(1, 0, 0.,0.,0., 0, "", 0)
        model.elements[1] = QUADElement(eid=1, prop_id=1, nodes=[1, 2, 3, 4])
        
        solver = LinearStaticSolver(model)
        
        with pytest.raises(NotImplementedError) as excinfo:
            solver.build_stiffness_matrix()
        
        assert "QUADElement" in str(excinfo.value)
    
    def test_error_unsupported_load(self):
        """@test Fail on Unsupported Load Type."""
        from src.preprocessing.data_model import LoadCase
        class ThermalLoad(LoadCase): pass

        model = FEMModel()
        model.nodes[1] = Node(1, 0, 0.,0.,0., 0, "", 0)
        model.loads.append(ThermalLoad(sid=1))
        
        solver = LinearStaticSolver(model)
        
        with pytest.raises(NotImplementedError) as excinfo:
            solver.build_load_vector()
        assert "ThermalLoad" in str(excinfo.value)

    def test_error_invalid_coordinate_system_load(self):
        """@test Fail if Load uses local CID!=0."""
        model = FEMModel()
        model.nodes[1] = Node(1, 0, 0.,0.,0., 0, "", 0)
        # Load with CID = 10 (Not Global)
        model.loads.append(SinglePointLoad(sid=1, appnode_id=1, module=10., vector=[1,0,0], cid=10))
        
        solver = LinearStaticSolver(model)
        
        with pytest.raises(ValueError) as excinfo:
            solver.build_load_vector()
        assert "CID=10" in str(excinfo.value)

    def test_error_unsupported_constraint(self):
        """@test Fail on Unsupported Constraint (e.g., MPC)."""
        from src.preprocessing.data_model import Constraint
        class MPC(Constraint): pass

        model = FEMModel()
        model.nodes[1] = Node(1, 0, 0.,0.,0., 0, "", 0)
        model.constraints.append(MPC(sid=10))
        
        solver = LinearStaticSolver(model)
        # Partitioning happens at the start of solve()
        with pytest.raises(NotImplementedError) as excinfo:
            solver.solve()
        assert "MPC" in str(excinfo.value)
    
    def test_singular_matrix_mechanism(self):
        """
        @test Singular Matrix Detection (Mechanism).
        @brief Verifies that the solver fails when the structure is unstable.
        @details
        Model:
            - Node 1: Pinned (Fixed Tx, Ty, Tz).
            - Node 2: Free.
            - Element: Rod along X.
        
        Why is it singular?
            A CROD element has axial and torsional stiffness, but NO bending stiffness.
            Node 2 is free to rotate around Y and Z axes (rigid body rotation).
            The stiffness matrix row corresponding to Ry and Rz at Node 2 will be zero.
            The solver should raise a singular matrix error.
        """
        model = FEMModel()
        model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0) # (0,0,0)
        model.nodes[2] = Node(2, 0, 1., 0., 0., 0, "", 0) # (1,0,0)
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 1.0, 0.0, 0.0
        
        # Setup standard property
        model.materials[1] = IsotropicMaterial(mid=1, E=1e9, nu=0.3)
        model.properties[1] = RodProperty(pid=1, mat_id=1, A=0.01, J=0.01, C=0.0, NSM=0.0)
        model.elements[1] = RodElement(eid=1, prop_id=1, nodes=[1, 2])

        # Constraint: Pin Node 1 (Tx, Ty, Tz fixed), but Rotations FREE
        # Even if we fix rotations at Node 1, Node 2 rotations are unconstrained
        # and uncoupled from axial/torsion forces in a Rod.
        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))
        
        # Load to try to force a solution
        model.loads.append(SinglePointLoad(sid=1, appnode_id=2, module=100.0, vector=[1,0,0]))

        solver = LinearStaticSolver(model)
        solver.assemble()

        # Scipy's spsolve raises a RuntimeError or MatrixRankWarning for singular matrices.
        # We catch generic Exception to be robust against scipy versions.
        with pytest.raises(Exception) as excinfo:
            solver.solve()
        
        # Verify the error is related to singularity/matrix issues
        error_msg = str(excinfo.value).lower()
        assert "singular" in error_msg or "unstable" in error_msg

    def test_multi_material_series_structure(self):
            """
            @test Multi-Material & Multi-Element Assembly.
            @brief Validates stiffness aggregation with heterogeneous properties.
            @details
            Structure: [Wall] --(Steel)-- [Node 2] --(Aluminum)-- [Node 3] -> F
            
            Data:
                - Elem 1 (Steel): L=1, A=1, E=200e9. k1 = 200e9
                - Elem 2 (Alum):  L=1, A=1, E=70e9.  k2 = 70e9
                - Load Fx = 1e5 N at Node 3.
            
            Theoretical Solution:
                The elements are springs in series.
                Force in both elements is F = 1e5.
                u2 = F / k1 = 1e5 / 200e9 = 5.0e-7 m
                extension_2 = F / k2 = 1e5 / 70e9 = 1.428e-6 m
                u3 = u2 + extension_2
            """
            model = FEMModel()
            
            # Nodes
            model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
            model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0., 0., 0.
            
            model.nodes[2] = Node(2, 0, 1., 0., 0., 0, "", 0)
            model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 1., 0., 0.
            
            model.nodes[3] = Node(3, 0, 2., 0., 0., 0, "", 0)
            model.nodes[3].x, model.nodes[3].y, model.nodes[3].z = 2., 0., 0.

            # Materials
            model.materials[1] = IsotropicMaterial(mid=1, E=200e9, nu=0.3) # Steel
            model.materials[2] = IsotropicMaterial(mid=2, E=70e9, nu=0.3)  # Aluminum

            # Properties (Same Area)
            model.properties[1] = RodProperty(pid=1, mat_id=1, A=1.0, J=0.1, C=0., NSM=0.)
            model.properties[2] = RodProperty(pid=2, mat_id=2, A=1.0, J=0.1, C=0., NSM=0.)

            # Elements
            model.elements[1] = RodElement(eid=1, prop_id=1, nodes=[1, 2]) # Steel Section
            model.elements[2] = RodElement(eid=2, prop_id=2, nodes=[2, 3]) # Aluminum Section

            # BCs: Fix Node 1. 
            model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))
            # Stabilize lateral DOFs for 2 and 3 (Roller)
            model.constraints.append(SPCConstraint(sid=1, node_id=2, components="23456", value=0.0))
            model.constraints.append(SPCConstraint(sid=1, node_id=3, components="23456", value=0.0))

            # Load
            F = 1.0e5
            model.loads.append(SinglePointLoad(sid=1, appnode_id=3, module=F, vector=[1,0,0]))

            # Solve
            solver = LinearStaticSolver(model)
            solver.assemble()
            u = solver.solve()

            # Validation
            idx_u2 = solver.mapper.get_dof_indices(2)[0] # Tx node 2
            idx_u3 = solver.mapper.get_dof_indices(3)[0] # Tx node 3
            
            u2_calc = u[idx_u2]
            u3_calc = u[idx_u3]
            
            k1 = (200e9 * 1.0) / 1.0
            k2 = (70e9 * 1.0) / 1.0
            
            u2_theory = F / k1
            u3_theory = u2_theory + (F / k2)
            
            assert u2_calc == pytest.approx(u2_theory, rel=1e-6)
            assert u3_calc == pytest.approx(u3_theory, rel=1e-6)
    def test_rotated_truss_structure(self):
        """
        @test Rotated Elements (Global Coordinate Assembly).
        @brief Validates that rotated matrices assemble correctly (The "Two-Bar Truss" problem).
        @details
        Geometry: Isosceles Triangle.
             (Fixed) N3
                |
                |  (Elem 2, Vertical)
                |
             (Free) N2  <-- Load Fx
               /
              / (Elem 1, 45 degrees)
             /
        (Fixed) N1
        
        Coordinates:
            N1: (0, 0)
            N2: (1, 1)
            N3: (1, 2)
        
        Load: Fx = 1000 N at N2.
        
        Logic:
            - Elem 2 is vertical: Resists Ty at N2 (Axial stiffness) and Tx via Shear/Bending?
              WAIT! CROD has NO shear/bending stiffness. Elem 2 only provides stiffness in Y.
              It offers NO resistance to movement in X (except negligible torsion/stability).
            - Elem 1 is 45 deg: Provides stiffness in both X and Y.
            
            Therefore, N2 should move significantly in X (restrained by Elem 1 projected stiffness)
            and also in Y (coupling due to 45 deg angle).
        """
        model = FEMModel()
        
        # Nodes
        model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
        model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0., 0., 0.

        model.nodes[2] = Node(2, 0, 1., 1., 0., 0, "", 0) # Junction
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 1., 1., 0.

        model.nodes[3] = Node(3, 0, 1., 2., 0., 0, "", 0)
        model.nodes[3].x, model.nodes[3].y, model.nodes[3].z = 1., 2., 0.

        # Common Property: EA = 1e6
        model.materials[1] = IsotropicMaterial(mid=1, E=1.0e6, nu=0.3)
        model.properties[1] = RodProperty(pid=1, mat_id=1, A=1.0, J=1.0, C=0., NSM=0.)

        # Elements
        model.elements[1] = RodElement(eid=1, prop_id=1, nodes=[1, 2]) # 45 deg
        model.elements[2] = RodElement(eid=2, prop_id=1, nodes=[3, 2]) # Vertical (Node 3 to 2)

        # BCs
        # Fix N1 and N3 fully
        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))
        model.constraints.append(SPCConstraint(sid=1, node_id=3, components="123456", value=0.0))
        
        # N2: Free in X and Y (Planar Truss). Fix Z and Rotations for stability
        model.constraints.append(SPCConstraint(sid=1, node_id=2, components="3456", value=0.0))

        # Load Fx on N2
        F_load = 1000.0
        model.loads.append(SinglePointLoad(sid=1, appnode_id=2, module=F_load, vector=[1,0,0]))

        # Solve
        solver = LinearStaticSolver(model)
        solver.assemble()
        u = solver.solve()

        # --- VALIDATION ---
        # Analytical Check:
        # Since Elem 2 is vertical rod, it has stiffness k2_y = EA/L, but k2_x = 0.
        # Elem 1 is 45 deg. Length L1 = sqrt(2). k1 = EA/sqrt(2).
        # Projects to Global K matrix at Node 2:
        # K1_global = k1 * [[c^2, cs], [cs, s^2]] where c=s=0.707 (cos^2 = 0.5)
        # K1_xx = 0.5 * k1, K1_yy = 0.5 * k1, K1_xy = 0.5 * k1
        #
        # Total Stiffness at Node 2:
        # K_xx_tot = K1_xx + K2_xx(0) = 0.5 * k1
        # K_yy_tot = K1_yy + K2_yy(EA/1) 
        # K_xy_tot = K1_xy
        #
        # System:
        # | Kxx  Kxy | | ux |   | F |
        # | Kxy  Kyy | | uy | = | 0 |
        
        E, A = 1.0e6, 1.0
        k1 = (E*A) / np.sqrt(2) # ~707106
        k2 = (E*A) / 1.0        # 1000000
        
        Kxx = 0.5 * k1
        Kyy = 0.5 * k1 + k2
        Kxy = 0.5 * k1
        
        # Solve 2x2 system analytically
        # Kxx*ux + Kxy*uy = F
        # Kxy*ux + Kyy*uy = 0  => uy = -(Kxy/Kyy)*ux
        #
        # Subs: Kxx*ux - (Kxy^2/Kyy)*ux = F
        # ux * (Kxx - Kxy^2/Kyy) = F
        
        det_sys = Kxx - (Kxy**2 / Kyy)
        expected_ux = F_load / det_sys
        expected_uy = -(Kxy / Kyy) * expected_ux
        
        idx_n2 = solver.mapper.get_dof_indices(2)
        ux_calc = u[idx_n2[0]]
        uy_calc = u[idx_n2[1]]
        
        assert ux_calc == pytest.approx(expected_ux, rel=1e-5)
        assert uy_calc == pytest.approx(expected_uy, rel=1e-5)



