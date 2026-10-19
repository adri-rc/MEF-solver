"""Integration tests for the linear static analysis workflow."""

import numpy as np
import pytest

from src.solver.linear_static_solver import LinearStaticSolver
from src.preprocessing.data_model import (
    FEMModel, Node, IsotropicMaterial, RodProperty, RodElement,
    BarProperty, BarElement,
    SPCConstraint, SinglePointLoad, SinglePointMoment
)


class TestLinearStaticSolverIntegration:

    def test_solve_simple_tension_rod(self):
        """
        A clamped rod with L=2, A=1, E=1000, and P=1000 has extension
        P*L/(A*E) = 2.
        """
        model = FEMModel()

        model.nodes[1] = Node(nid=1, coord_id=0, x1=0., x2=0., x3=0., cd_id=0, ps="", seid=0)
        model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0., 0., 0.
        
        model.nodes[2] = Node(nid=2, coord_id=0, x1=2., x2=0., x3=0., cd_id=0, ps="", seid=0)
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 2., 0., 0.

        model.materials[1] = IsotropicMaterial(mid=1, E=1000.0, nu=0.3)
        model.properties[1] = RodProperty(pid=1, mat_id=1, A=1.0, J=1.0, C=0.0, NSM=0.0)

        model.elements[10] = RodElement(eid=10, prop_id=1, nodes=[1, 2])

        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))
        model.constraints.append(SPCConstraint(sid=1, node_id=2, components="23456", value=0.0))

        model.loads.append(SinglePointLoad(sid=1, appnode_id=2, module=1000.0, vector=[1.0, 0.0, 0.0]))

        solver = LinearStaticSolver(model)
        solver.assemble()
        results = solver.solve()
        displacements = results['displacements']
        
        dofs_n2 = solver.dof_mapper.get_dof_indices(2)
        u_x_n2 = displacements[dofs_n2[0]]

        expected_u = 2.0
        assert u_x_n2 == pytest.approx(expected_u, rel=1e-5)

    def test_solve_enforced_displacement(self):
        """Apply a nonzero axial SPC without external forces."""
        model = FEMModel()
        model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
        model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0.0, 0.0, 0.0
        model.nodes[2] = Node(2, 0, 1.0, 0., 0., 0, "", 0)
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 1.0, 0.0, 0.0

        model.materials[1] = IsotropicMaterial(mid=1, E=1000.0, nu=0.3)
        model.properties[1] = RodProperty(pid=1, mat_id=1, A=1.0, J=1.0, C=0.0, NSM=0.0)
        model.elements[1] = RodElement(eid=1, prop_id=1, nodes=[1, 2])

        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))

        model.constraints.append(SPCConstraint(sid=1, node_id=2, components="1", value=0.01))
        # Stabilize the zero-stiffness transverse and rotational DOFs while
        # retaining one free DOF so the model is not fully constrained.
        model.constraints.append(SPCConstraint(sid=1, node_id=2, components="2356", value=0.0))

        solver = LinearStaticSolver(model)
        solver.assemble()
        results = solver.solve()
        displacements = results['displacements']

        dofs_n2 = solver.dof_mapper.get_dof_indices(2)
        u_x_n2 = displacements[dofs_n2[0]]

        assert u_x_n2 == pytest.approx(0.01, abs=1e-9)
    def test_error_unsupported_load(self):
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
        model = FEMModel()
        model.nodes[1] = Node(1, 0, 0.,0.,0., 0, "", 0)
        model.loads.append(SinglePointLoad(sid=1, appnode_id=1, module=10., vector=[1,0,0], cid=10))
        
        solver = LinearStaticSolver(model)
        
        with pytest.raises(ValueError) as excinfo:
            solver.build_load_vector()
        assert "CID=10" in str(excinfo.value)

    def test_error_unsupported_constraint(self):
        from src.preprocessing.data_model import Constraint
        class MPC(Constraint): pass

        model = FEMModel()
        model.nodes[1] = Node(1, 0, 0.,0.,0., 0, "", 0)
        model.constraints.append(MPC(sid=10))
        
        solver = LinearStaticSolver(model)
        with pytest.raises(NotImplementedError) as excinfo:
            solver.solve()
        assert "MPC" in str(excinfo.value)
    
    def test_singular_matrix_mechanism(self):
        """
        A CROD has no bending stiffness, so free Ry/Rz modes at node 2 must
        make the system singular rather than trigger implicit stabilization.
        """
        model = FEMModel()
        model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
        model.nodes[2] = Node(2, 0, 1., 0., 0., 0, "", 0)
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 1.0, 0.0, 0.0
        
        model.materials[1] = IsotropicMaterial(mid=1, E=1e9, nu=0.3)
        model.properties[1] = RodProperty(pid=1, mat_id=1, A=0.01, J=0.01, C=0.0, NSM=0.0)
        model.elements[1] = RodElement(eid=1, prop_id=1, nodes=[1, 2])

        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))

        model.loads.append(SinglePointLoad(sid=1, appnode_id=2, module=100.0, vector=[1,0,0]))

        solver = LinearStaticSolver(model)
        solver.assemble()

        # SciPy exposes singularity differently across supported versions.
        with pytest.raises(Exception) as excinfo:
            solver.solve()

        error_msg = str(excinfo.value).lower()
        assert "singular" in error_msg or "unstable" in error_msg

    def test_multi_material_series_structure(self):
            """
            Two heterogeneous rods in series carry the same force. Their
            extensions therefore add as F/k1 + F/k2.
            """
            model = FEMModel()
            
            model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
            model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0., 0., 0.
            
            model.nodes[2] = Node(2, 0, 1., 0., 0., 0, "", 0)
            model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 1., 0., 0.
            
            model.nodes[3] = Node(3, 0, 2., 0., 0., 0, "", 0)
            model.nodes[3].x, model.nodes[3].y, model.nodes[3].z = 2., 0., 0.

            model.materials[1] = IsotropicMaterial(mid=1, E=200e9, nu=0.3)
            model.materials[2] = IsotropicMaterial(mid=2, E=70e9, nu=0.3)

            model.properties[1] = RodProperty(pid=1, mat_id=1, A=1.0, J=0.1, C=0., NSM=0.)
            model.properties[2] = RodProperty(pid=2, mat_id=2, A=1.0, J=0.1, C=0., NSM=0.)

            model.elements[1] = RodElement(eid=1, prop_id=1, nodes=[1, 2])
            model.elements[2] = RodElement(eid=2, prop_id=2, nodes=[2, 3])

            model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))
            model.constraints.append(SPCConstraint(sid=1, node_id=2, components="23456", value=0.0))
            model.constraints.append(SPCConstraint(sid=1, node_id=3, components="23456", value=0.0))

            F = 1.0e5
            model.loads.append(SinglePointLoad(sid=1, appnode_id=3, module=F, vector=[1,0,0]))

            solver = LinearStaticSolver(model)
            solver.assemble()
            results = solver.solve()
            u = results['displacements']

            idx_u2 = solver.dof_mapper.get_dof_indices(2)[0]
            idx_u3 = solver.dof_mapper.get_dof_indices(3)[0]
            
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
        Assemble a 45-degree rod and a vertical rod at their common node. The
        vertical CROD contributes only global-Y stiffness, while the inclined
        CROD couples the global X and Y translations.
        """
        model = FEMModel()

        model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
        model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0., 0., 0.

        model.nodes[2] = Node(2, 0, 1., 1., 0., 0, "", 0)
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 1., 1., 0.

        model.nodes[3] = Node(3, 0, 1., 2., 0., 0, "", 0)
        model.nodes[3].x, model.nodes[3].y, model.nodes[3].z = 1., 2., 0.

        model.materials[1] = IsotropicMaterial(mid=1, E=1.0e6, nu=0.3)
        model.properties[1] = RodProperty(pid=1, mat_id=1, A=1.0, J=1.0, C=0., NSM=0.)

        model.elements[1] = RodElement(eid=1, prop_id=1, nodes=[1, 2])
        model.elements[2] = RodElement(eid=2, prop_id=1, nodes=[3, 2])

        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))
        model.constraints.append(SPCConstraint(sid=1, node_id=3, components="123456", value=0.0))

        model.constraints.append(SPCConstraint(sid=1, node_id=2, components="3456", value=0.0))

        F_load = 1000.0
        model.loads.append(SinglePointLoad(sid=1, appnode_id=2, module=F_load, vector=[1,0,0]))

        solver = LinearStaticSolver(model)
        solver.assemble()
        results = solver.solve()
        u = results['displacements']

        # The inclined member contributes k1/2 to Kxx, Kyy, and Kxy; the
        # vertical member adds k2 only to Kyy.
        E, A = 1.0e6, 1.0
        k1 = (E*A) / np.sqrt(2)
        k2 = (E*A) / 1.0

        Kxx = 0.5 * k1
        Kyy = 0.5 * k1 + k2
        Kxy = 0.5 * k1
        
        det_sys = Kxx - (Kxy**2 / Kyy)
        expected_ux = F_load / det_sys
        expected_uy = -(Kxy / Kyy) * expected_ux
        
        idx_n2 = solver.dof_mapper.get_dof_indices(2)
        ux_calc = u[idx_n2[0]]
        uy_calc = u[idx_n2[1]]
        
        assert ux_calc == pytest.approx(expected_ux, rel=1e-5)
        assert uy_calc == pytest.approx(expected_uy, rel=1e-5)

    def test_solve_reactions_and_stresses(self):
        """For a cantilever rod, check R=-F and sigma=F/A."""
        model = FEMModel()
        model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
        model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0., 0., 0.
        model.nodes[2] = Node(2, 0, 2., 0., 0., 0, "", 0)
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 2., 0., 0.

        model.materials[1] = IsotropicMaterial(mid=1, E=1000.0, nu=0.3)
        model.properties[1] = RodProperty(pid=1, mat_id=1, A=10.0, J=1.0, C=0.0, NSM=0.0)
        model.elements[10] = RodElement(eid=10, prop_id=1, nodes=[1, 2])

        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))
        model.constraints.append(SPCConstraint(sid=1, node_id=2, components="23456", value=0.0))

        model.loads.append(SinglePointLoad(sid=1, appnode_id=2, module=500.0, vector=[1.0, 0.0, 0.0]))

        solver = LinearStaticSolver(model)
        solver.assemble()
        solver.solve()

        results = solver.reactions_and_stresses()

        r_node1 = results['reactions'][1]
        assert r_node1[0] == pytest.approx(-500.0, abs=1e-9)

        elem_res = results['stresses'][10]
        
        assert elem_res['axial_force'] == pytest.approx(500.0)
        assert elem_res['axial_stress'] == pytest.approx(50.0)
        assert elem_res['torque'] == pytest.approx(0.0)


    def test_reaction_inclined_rod(self):
        """
        Resolve equilibrium reactions for a 45-degree rod loaded in global X.
        The wall carries (-F, -F), while the Y roller carries (0, +F).
        """
        model = FEMModel()
        model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
        model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0., 0., 0.
        model.nodes[2] = Node(2, 0, 1., 1., 0., 0, "", 0)
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 1., 1., 0.

        model.materials[1] = IsotropicMaterial(mid=1, E=1000.0, nu=0.3)
        model.properties[1] = RodProperty(pid=1, mat_id=1, A=1.0, J=1.0, C=0., NSM=0.)
        model.elements[1] = RodElement(eid=1, prop_id=1, nodes=[1, 2])

        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))

        model.constraints.append(SPCConstraint(sid=1, node_id=2, components="23456", value=0.0))

        F_app = 1000.0
        model.loads.append(SinglePointLoad(sid=1, appnode_id=2, module=F_app, vector=[1, 0, 0]))

        solver = LinearStaticSolver(model)
        solver.assemble()
        solver.solve()
        
        results = solver.reactions_and_stresses()
        
        r1 = results['reactions'][1]
        assert r1[0] == pytest.approx(-F_app, abs=1e-9), "Reaction Rx at Node 1 incorrect"
        assert r1[1] == pytest.approx(-F_app, abs=1e-9), "Reaction Ry at Node 1 incorrect"

        r2 = results['reactions'][2]
        assert r2[1] == pytest.approx(F_app, abs=1e-9), "Reaction Ry at Node 2 incorrect"

    def test_reaction_3d_arbitrary_rod(self):
        """
        A rod along the 3-4-12 direction has length 13. Applying the aligned
        vector (300, 400, 1200) gives an axial force of 1300.
        """
        model = FEMModel()

        model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
        model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0.0, 0.0, 0.0
        
        model.nodes[2] = Node(2, 0, 3., 4., 12., 0, "", 0)
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = 3.0, 4.0, 12.0

        model.materials[1] = IsotropicMaterial(mid=1, E=13000.0, nu=0.3)
        model.properties[1] = RodProperty(pid=1, mat_id=1, A=1.0, J=1.0, C=0., NSM=0.)
        model.elements[10] = RodElement(eid=10, prop_id=1, nodes=[1, 2])

        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))

        # Stabilize zero-stiffness transverse and rotational directions.
        model.constraints.append(SPCConstraint(sid=1, node_id=2, components="23456", value=0.0))

        Fx, Fy, Fz = 300.0, 400.0, 1200.0
        model.loads.append(SinglePointLoad(sid=1, appnode_id=2, module=1.0, vector=[Fx, Fy, Fz]))

        solver = LinearStaticSolver(model)
        solver.assemble()
        results = solver.solve()
        
        r1 = results['reactions'][1]
        
        assert r1[0] == pytest.approx(-Fx, abs=1e-9), "Reaction Rx incorrect"
        assert r1[1] == pytest.approx(-Fy, abs=1e-9), "Reaction Ry incorrect"
        assert r1[2] == pytest.approx(-Fz, abs=1e-9), "Reaction Rz incorrect"

        r2 = results['reactions'][2]
        assert r2[1] == pytest.approx(0.0, abs=1e-9)
        assert r2[2] == pytest.approx(0.0, abs=1e-9)
        
        elem_stress = results['stresses'][10]
        expected_force = 1300.0
        
        assert elem_stress['axial_force'] == pytest.approx(expected_force, abs=1e-9)
        
        # The 1.3 m elongation is the X displacement projected by 3/13.
        idx_n2 = solver.dof_mapper.get_dof_indices(2)
        u = results['displacements']
        
        expected_ux = 16.9 / 3.0
        assert u[idx_n2[0]] == pytest.approx(expected_ux, rel=1e-5)
        assert u[idx_n2[1]] == pytest.approx(0.0, abs=1e-9)
        assert u[idx_n2[2]] == pytest.approx(0.0, abs=1e-9)
        

class TestBarElementIntegration:
    def test_solve_cantilever_beam_tip_load(self):
        """
        Validate the Euler-Bernoulli cantilever tip displacement PL**3/(3EI)
        and the root moment under a transverse tip load.
        """
        model = FEMModel()
        L = 2.0
        E = 210e9
        I1 = 1.0e-5
        
        model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
        model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0., 0., 0.
        
        model.nodes[2] = Node(2, 0, L, 0., 0., 0, "", 0)
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = L, 0., 0.

        model.materials[1] = IsotropicMaterial(mid=1, E=E, nu=0.3)
        model.properties[1] = BarProperty(pid=1, mat_id=1, A=0.1, I1=I1, I2=1e-5, J=1e-5)

        model.elements[10] = BarElement(eid=10, prop_id=1, nodes=[1, 2], orientation_vector=[0., 1., 0.])

        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))

        Fy = -1000.0
        model.loads.append(SinglePointLoad(sid=1, appnode_id=2, module=Fy, vector=[0.0, 1.0, 0.0]))

        solver = LinearStaticSolver(model)
        solver.assemble()
        res = solver.solve()
        
        expected_v = (Fy * L**3) / (3 * E * I1)

        dofs_n2 = solver.dof_mapper.get_dof_indices(2)
        v_calc = res['displacements'][dofs_n2[1]]

        assert v_calc == pytest.approx(expected_v, rel=1e-4)

        stresses = res['stresses'][10]
        assert stresses['type'] == 'CBAR'

        mz_n1 = stresses['forces_n1']['Mz']
        assert mz_n1 == pytest.approx(-Fy * L, rel=1e-4)

    def test_solve_torsion_bar(self):
        """Validate the torsional rotation theta=T*L/(G*J)."""
        model = FEMModel()
        L = 1.0
        G = 80e9
        J = 1.0e-4
        
        model.nodes[1] = Node(1, 0, 0., 0., 0., 0, "", 0)
        model.nodes[1].x, model.nodes[1].y, model.nodes[1].z = 0., 0., 0.
        model.nodes[2] = Node(2, 0, L, 0., 0., 0, "", 0)
        model.nodes[2].x, model.nodes[2].y, model.nodes[2].z = L, 0., 0.

        model.materials[1] = IsotropicMaterial(mid=1, E=200e9, G=G, nu=0.3)
        model.properties[1] = BarProperty(pid=1, mat_id=1, A=0.1, I1=1e-4, I2=1e-4, J=J)
        
        model.elements[20] = BarElement(eid=20, prop_id=1, nodes=[1, 2], orientation_vector=[0., 1., 0.])

        model.constraints.append(SPCConstraint(sid=1, node_id=1, components="123456", value=0.0))

        Mx = 500.0
        model.loads.append(SinglePointMoment(sid=1, appnode_id=2, module=Mx, vector=[1.0, 0.0, 0.0]))
        
        solver = LinearStaticSolver(model)
        solver.assemble()
        res = solver.solve()
        
        expected_theta = (Mx * L) / (G * J)

        dofs_n2 = solver.dof_mapper.get_dof_indices(2)
        theta_calc = res['displacements'][dofs_n2[3]]

        assert theta_calc == pytest.approx(expected_theta, rel=1e-4)

        stresses = res['stresses'][20]
        assert stresses['torque'] == pytest.approx(Mx, rel=1e-4)
