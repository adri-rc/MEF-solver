"""Linear CTRIA3 membrane and Reissner-Mindlin shell formulation.

PSHELL MID1, MID2, and MID3 independently control membrane, bending, and
transverse-shear response; blank fields do not inherit another material.
"""

from typing import Optional

import numpy as np

from src.preprocessing.data_model import IsotropicMaterial, Node, ShellProperty

from ..transformations import CoordinateTransformation


class CTRIA3ElementSolver:
    """Compute the 18-DOF response of a three-node triangular shell."""
    
    def __init__(
        self,
        node1: Node,
        node2: Node,
        node3: Node,
        prop: ShellProperty,
        mat: Optional[IsotropicMaterial],
        mat_bending: Optional[IsotropicMaterial] = None,
        mat_shear: Optional[IsotropicMaterial] = None,
        k6rot: float = 0.0,
        theta: float = 0.0,
    ):
        """Initialize independent PSHELL material and scaling fields."""
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.prop = prop
        self.k6rot = k6rot
        self.theta = float(theta)
        self.mat = mat  # MID1
        self.mat_bending = mat_bending  # MID2
        self.mat_shear = mat_shear  # MID3
        if self.mat is None and self.mat_bending is None:
            raise ValueError("CTRIA3 requires a MID1 or MID2 material.")
        if self.mat_shear is not None and self.mat_bending is None:
            raise ValueError("CTRIA3 cannot define MID3 without MID2.")
        
        self.n1_coords = node1.get_coords()
        self.n2_coords = node2.get_coords()
        self.n3_coords = node3.get_coords()
        
        self._compute_transformation()

        # Local coordinates are N1=(0, 0), N2=(x2, 0), and N3=(x3, y3).
        v12 = self.n2_coords - self.n1_coords
        self.x2 = np.linalg.norm(v12)

        v13 = self.n3_coords - self.n1_coords
        self.x3 = np.dot(v13, self.local_x_axis)
        self.y3 = np.dot(v13, self.local_y_axis)

        self.A = 0.5 * self.x2 * self.y3
        if self.A < 1e-12:
            raise ValueError(f"Zero area triangle: Nodes {node1.nid}, {node2.nid}, {node3.nid}")

    def _compute_transformation(self):
        """Build the right-handed 18-by-18 element transformation."""
        v12 = self.n2_coords - self.n1_coords
        v13 = self.n3_coords - self.n1_coords

        # Local x runs from node 1 to node 2.
        length_12 = np.linalg.norm(v12)
        if length_12 < 1e-12:
            raise ValueError(f"Zero length edge in CTRIA3: Nodes {self.node1.nid}-{self.node2.nid}")
        vx = v12 / length_12

        # The right-hand-rule normal fixes the physical top-face direction.
        vz_cross = np.cross(vx, v13)
        normal_norm = np.linalg.norm(vz_cross)
        if normal_norm < 1e-12:
            raise ValueError(
                f"Zero area triangle: Nodes {self.node1.nid}, {self.node2.nid}, {self.node3.nid}"
            )
        vz = vz_cross / normal_norm

        # Local y completes the right-handed triad.
        vy = np.cross(vz, vx)

        self.local_x_axis = vx
        self.local_y_axis = vy
        self.local_z_axis = vz

        R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)
        self.T = CoordinateTransformation.build_full_transformation_matrix(
            R, nodes_count=3, dofs_per_node=6
        )

    def _get_constitutive_matrices(self):
        """Build independent membrane, bending, and shear matrices."""
        # Membrane response (MID1).
        D_m = np.zeros((3, 3))
        if self.mat is not None:
            E_m, nu_m = self.mat.E, self.mat.nu
            factor_m = E_m / (1.0 - nu_m**2)
            D_m = factor_m * np.array([
                [1.0,  nu_m, 0.0],
                [nu_m, 1.0,  0.0],
                [0.0,  0.0,  (1.0 - nu_m) / 2.0]
            ])

        # Bending response (MID2).
        t = self.prop.t
        D_f = np.zeros((3, 3))
        if self.mat_bending is not None:
            E_f, nu_f = self.mat_bending.E, self.mat_bending.nu
            bending_ratio = float(
                getattr(self.prop, "bending_inertia_ratio", 1.0)
            )
            factor_f = (
                bending_ratio
                * E_f
                * t**3
                / (12.0 * (1.0 - nu_f**2))
            )
            D_f = factor_f * np.array([
                [1.0,  nu_f, 0.0],
                [nu_f, 1.0,  0.0],
                [0.0,  0.0,  (1.0 - nu_f) / 2.0]
            ])

        # Transverse-shear response (MID3).
        D_c = np.zeros((2, 2))
        if self.mat_shear is not None:
            E_s, nu_s = self.mat_shear.E, self.mat_shear.nu
            G_c = getattr(self.mat_shear, 'G', E_s / (2.0 * (1.0 + nu_s)))
            shear_ratio = float(
                getattr(self.prop, "shear_thickness_ratio", 5.0 / 6.0)
            )
            D_c = (shear_ratio * G_c * t) * np.eye(2)

        return D_m, D_f, D_c
    
    def _get_isoparametric_b_matrices(self):
        """Build the strain-displacement matrices at the centroid."""
        x = np.array([0.0, self.x2, self.x3])
        y = np.array([0.0, 0.0, self.y3])
        
        dN_dxi  = np.array([-1.0,  1.0,  0.0])
        dN_deta = np.array([-1.0,  0.0,  1.0])
        
        dx_dxi, dy_dxi   = np.dot(dN_dxi, x), np.dot(dN_dxi, y)
        dx_deta, dy_deta = np.dot(dN_deta, x), np.dot(dN_deta, y)
        
        J = np.array([
            [dx_dxi, dy_dxi],
            [dx_deta, dy_deta]
        ])
        
        detJ = np.abs(np.linalg.det(J))
        J_inv = np.linalg.inv(J)
        
        dN_dnat = np.vstack((dN_dxi, dN_deta))
        dN_dphysical = np.dot(J_inv, dN_dnat)
        dN_dx, dN_dy = dN_dphysical[0, :], dN_dphysical[1, :]
        
        B_m = np.zeros((3, 6))
        for i in range(3):
            B_m[0, 2*i]   = dN_dx[i]
            B_m[1, 2*i+1] = dN_dy[i]
            B_m[2, 2*i]   = dN_dy[i]
            B_m[2, 2*i+1] = dN_dx[i]

        xi, eta = 1.0/3.0, 1.0/3.0
        N = np.array([1.0 - xi - eta, xi, eta])
        
        B_f = np.zeros((3, 9))
        B_c = np.zeros((2, 9))
        
        for i in range(3):
            col = 3 * i
            B_f[0, col+2] = dN_dx[i]
            B_f[1, col+1] = -dN_dy[i]
            B_f[2, col+1] = -dN_dx[i]
            B_f[2, col+2] = dN_dy[i]
            
            B_c[0, col]   = dN_dx[i]
            B_c[0, col+2] = N[i]
            B_c[1, col]   = dN_dy[i]
            B_c[1, col+1] = -N[i]

        return B_m, B_f, B_c, detJ

    def _get_B_bending_Nastran(self):
        return self._get_isoparametric_b_matrices()[1]

    def _get_B_shear_Nastran(self, *_natural_coordinates):
        return self._get_isoparametric_b_matrices()[2]

    def _get_D_bending(self):
        return self._get_constitutive_matrices()[1]

    def _get_D_shear(self):
        return self._get_constitutive_matrices()[2]

    def _build_local_stiffness(self) -> np.ndarray:
        """Build the 18-by-18 local stiffness matrix."""
        K_local = np.zeros((18, 18))
        gauss_weight = 0.5
        
        D_m, D_f, D_c = self._get_constitutive_matrices()
        B_m, B_f, B_c, detJ = self._get_isoparametric_b_matrices()

        # Membrane stiffness.
        k_mem = (self.prop.t * detJ * gauss_weight) * (B_m.T @ D_m @ B_m)
        idx_m = [0, 1, 6, 7, 12, 13]
        for i in range(6):
            for j in range(6):
                K_local[idx_m[i], idx_m[j]] += k_mem[i, j]
                
        # Bending and transverse-shear stiffness.
        has_bending = self.prop.t > 0.0 and self.mat_bending is not None
        if has_bending:
            bending_stiffness = (detJ * gauss_weight) * (B_f.T @ D_f @ B_f)
            shear_stiffness = (detJ * gauss_weight) * (B_c.T @ D_c @ B_c)
            k_plate = bending_stiffness + shear_stiffness
            
            idx_f = [2, 3, 4, 8, 9, 10, 14, 15, 16]
            for i in range(9):
                for j in range(9):
                    K_local[idx_f[i], idx_f[j]] += k_plate[i, j]

        # K6ROT stabilizes the drilling rotation of the membrane block; without
        # MID1 there is no membrane material from which to derive its penalty.
        if self.k6rot > 0.0 and self.mat is not None:
            nu = self.mat.nu
            G = self.mat.E / (2.0 * (1.0 + nu))

            # Nastran penalty: K6ROT * 1e-6 * G * t * detJ.
            k_theta = self.k6rot * 1e-6 * G * self.prop.t * detJ

            # Macro-rotation: omega = 0.5 * (dv/dx - du/dy).
            b_omega = np.array([
                -0.5 * B_m[1, 1],  0.5 * B_m[0, 0],
                -0.5 * B_m[1, 3],  0.5 * B_m[0, 2],
                -0.5 * B_m[1, 5],  0.5 * B_m[0, 4],
            ])

            B_drill = np.zeros((3, 9))
            for i in range(3):
                B_drill[i, 0] = -b_omega[0] # u1
                B_drill[i, 1] = -b_omega[1] # v1
                B_drill[i, 3] = -b_omega[2] # u2
                B_drill[i, 4] = -b_omega[3] # v2
                B_drill[i, 6] = -b_omega[4] # u3
                B_drill[i, 7] = -b_omega[5] # v3
                
                # Add the drilling rotation of node i.
                B_drill[i, 2 + i*3] = 1.0

            K_drill_local = k_theta * (B_drill.T @ B_drill)

            # [u1, v1, theta_z1, u2, v2, theta_z2, u3, v3, theta_z3]
            idx_drill = [0, 1, 5, 6, 7, 11, 12, 13, 17]
            for i in range(9):
                for j in range(9):
                    K_local[idx_drill[i], idx_drill[j]] += K_drill_local[i, j]
        
        return K_local

    def get_stiffness_matrix(self) -> np.ndarray:
        """Return the symmetric stiffness matrix in global coordinates."""
        K_loc = self._build_local_stiffness()
        stiffness = self.T.T @ K_loc @ self.T
        return 0.5 * (stiffness + stiffness.T)

    def calculate_stresses(self, u_global_vector: np.ndarray, dof_mapper) -> dict:
        """Recover shell stresses and force resultants at the centroid."""
        idxs = []
        for n in [self.node1, self.node2, self.node3]:
            idxs.extend(dof_mapper.get_dof_indices(n.nid))
            
        u_elem_global = u_global_vector[idxs]
        u_loc = self.T @ u_elem_global
        
        D_m, D_f, D_c = self._get_constitutive_matrices()
        B_m, B_f, B_c, detJ = self._get_isoparametric_b_matrices()
        
        # Membrane response.
        u_mem = np.array([u_loc[0], u_loc[1], u_loc[6], u_loc[7], u_loc[12], u_loc[13]])
        strains = B_m @ u_mem 
        stresses = D_m @ strains 
        sig_x, sig_y, tau_xy = stresses

        # Nastran reports CTRIA3 force/stress components in the element
        # coordinate system. THETA defines material orientation and must not
        # rotate isotropic MAT1 recovery output.
        
        t = self.prop.t
        f_x, f_y, f_xy = sig_x * t, sig_y * t, tau_xy * t
        
        # Bending and transverse-shear response.
        has_bending = self.prop.t > 0.0 and self.mat_bending is not None
        if has_bending:
            u_bend_RM = np.array([u_loc[2], u_loc[3], u_loc[4], u_loc[8], u_loc[9], u_loc[10], u_loc[14], u_loc[15], u_loc[16]])
            kappa = B_f @ u_bend_RM
            moments = D_f @ kappa
            # Match the Nastran element-result sign convention.
            m_x, m_y, m_xy = -moments[0], -moments[1], -moments[2]
            
            gamma_s = B_c @ u_bend_RM
            shear_forces = D_c @ gamma_s
            q_x, q_y = shear_forces[0], shear_forces[1]
        else:
            m_x, m_y, m_xy, q_x, q_y = 0.0, 0.0, 0.0, 0.0, 0.0
            
        # Public convention: top is local z=+t/2 and bottom is local z=-t/2.
        # Avoid Z1/Z2 aliases because their meaning depends on the external
        # solver's shell-result convention.
        bending_ratio = float(
            getattr(self.prop, "bending_inertia_ratio", 1.0)
        )
        section_modulus = bending_ratio * t**2 / 6.0
        sb_x = m_x / section_modulus
        sb_y = m_y / section_modulus
        sb_xy = m_xy / section_modulus
        sig_x_top, sig_x_bot = sig_x - sb_x, sig_x + sb_x
        sig_y_top, sig_y_bot = sig_y - sb_y, sig_y + sb_y
        tau_xy_top, tau_xy_bot = tau_xy - sb_xy, tau_xy + sb_xy

        # Plane-stress invariants.
        def calc_principal_vm(sx, sy, txy):
            vm = np.sqrt(sx**2 - sx*sy + sy**2 + 3*txy**2)
            s_avg, s_diff = (sx + sy) / 2.0, (sx - sy) / 2.0
            R = np.sqrt(s_diff**2 + txy**2)
            return vm, s_avg + R, s_avg - R

        # Evaluate mid-plane and surface fibers.
        vm_mem, s_maj_mem, s_min_mem = calc_principal_vm(sig_x, sig_y, tau_xy)
        vm_top, s_maj_top, s_min_top = calc_principal_vm(sig_x_top, sig_y_top, tau_xy_top)
        vm_bot, s_maj_bot, s_min_bot = calc_principal_vm(sig_x_bot, sig_y_bot, tau_xy_bot)
        
        return {
            # Membrane / mid-plane.
            'sigma_x': sig_x, 'sigma_y': sig_y, 'tau_xy': tau_xy,
            'von_mises': vm_mem, 'major_principal': s_maj_mem, 'minor_principal': s_min_mem,

            # Section forces and moments.
            'membrane_force_x': f_x, 'membrane_force_y': f_y, 'membrane_force_xy': f_xy,
            'shear_force_x': q_x, 'shear_force_y': q_y,
            'bending_moment_x': m_x, 'bending_moment_y': m_y, 'bending_moment_xy': m_xy,

            # Top Fiber (local z=+t/2)
            'stress_top_fiber_x': sig_x_top, 'stress_top_fiber_y': sig_y_top, 'tau_top_fiber_xy': tau_xy_top,
            'von_mises_top': vm_top, 'major_principal_top': s_maj_top, 'minor_principal_top': s_min_top,

            # Bottom Fiber (local z=-t/2)
            'stress_bot_fiber_x': sig_x_bot, 'stress_bot_fiber_y': sig_y_bot, 'tau_bot_fiber_xy': tau_xy_bot,
            'von_mises_bot': vm_bot, 'major_principal_bot': s_maj_bot, 'minor_principal_bot': s_min_bot
        }
