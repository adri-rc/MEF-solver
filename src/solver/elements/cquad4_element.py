"""
@file cquad4_element.py
@brief Shared mechanics for four-node quadrilateral shell elements.
@details
Contains only geometry, coordinate transformations, constitutive matrices and
stress recovery shared by the independent Reissner-Mindlin and Nastran
implementations.  Stiffness formulations live in their own modules."""

from abc import ABC, abstractmethod
import numpy as np
from src.preprocessing.data_model import Node, ShellProperty, IsotropicMaterial

from ..quadrature import GaussQuadrature
from ..transformations import CoordinateTransformation

class CQUAD4ElementBase(ABC):
    """
    @class CQUAD4ElementBase
    @brief Common logic for a 4-node quadrilateral shell element.
    @details
    Computes the Global Stiffness Matrix (24x24) for a quadrilateral shell.
    
    Supported DOFs per node: 6 (Tx, Ty, Tz, Rx, Ry, Rz).
    Total DOFs: 24.
    
    Local Coordinate System (Nastran standard):
      z: Normal to the plane defined by the two diagonals (N1->N3 x N2->N4).
      x: Bisector of the two element diagonals.
      y: Perpendicular in plane (z cross x).
      Origin: Mean of the 4 nodes.
    
    The concrete stiffness formulation is supplied by a subclass.
    """

    # --- CONSTANTS FOR DOF MAPPING ---
    # Nodes: 1, 2, 3, 4
    # DOFs per node: 0-5
    # Global Local Vector indices: 
    # Node 1: 0-5, Node 2: 6-11, Node 3: 12-17, Node 4: 18-23
    
    def __init__(self, node1: Node, node2: Node, node3: Node, node4: Node, 
                 prop: ShellProperty, mat_membrane: IsotropicMaterial, 
                 mat_bending: IsotropicMaterial = None, mat_shear: IsotropicMaterial = None,
                 k6rot: float = 0.0, theta: float = 0.0):
        """
        @brief Initialize the CQUAD4 Element Solver.
        @param mat_bending Optional material for bending stiffness (if different from membrane).
        """
        self.nodes = [node1, node2, node3, node4]
        self.prop = prop
        self.k6rot = float(k6rot)
        self.theta = float(theta)
        
        # Asignación explícita de materiales
        self.mat_membrane = mat_membrane
        self.mat_bending = mat_bending if mat_bending is not None else self.mat_membrane
        # A blank PSHELL MID3 disables physical transverse-shear flexibility.
        # It must not silently inherit MID2; the Nastran QUAD4 retains only its
        # residual bending flexibility in that case.
        self.mat_shear = mat_shear
        
        self.coords = np.array([n.get_coords() for n in self.nodes])
        
        # 1. Compute Local System and Transformation Matrix
        self._compute_transformation()
        
        # 2. Compute Geometric Properties in Local System
        # Express nodes in local 2D plane (z~0 for slightly warped elements, exact 0 for flat)
        # 2. Compute Geometric Properties in Local System
        self.local_coords_2d = np.zeros((4, 2))
        for i in range(4):
            v = self.coords[i] - self.origin
            self.local_coords_2d[i, 0] = np.dot(v, self.local_x_axis)
            self.local_coords_2d[i, 1] = np.dot(v, self.local_y_axis)
        self.x_loc = self.local_coords_2d[:, 0]
        self.y_loc = self.local_coords_2d[:, 1]
            
    def _compute_transformation(self):
        """
        @brief Computes T matrix (24x24) based on Nastran standard local axes for CQUAD4.
        """
        self.origin = np.mean(self.coords, axis=0)
        
        # Z-axis: Cross product of diagonals
        v13 = self.coords[2] - self.coords[0]
        v24 = self.coords[3] - self.coords[1]
        vz = np.cross(v13, v24)
        vz_norm = np.linalg.norm(vz)
        if vz_norm < 1e-12:
            raise ValueError(f"Degenerate CQUAD4 element. Check node numbering.")
        vz = vz / vz_norm
        
        # X-axis: Nastran's CQUAD4 element x-axis bisects the angle between
        # diagonals G1-G3 and G4-G2.  Normalising each diagonal before adding
        # is essential for an irregular quadrilateral.
        diagonal_13 = self.coords[2] - self.coords[0]
        diagonal_42 = self.coords[1] - self.coords[3]
        diagonal_13 /= np.linalg.norm(diagonal_13)
        diagonal_42 /= np.linalg.norm(diagonal_42)
        vx = diagonal_13 + diagonal_42
        vx_norm = np.linalg.norm(vx)
        if vx_norm < 1e-12:
            raise ValueError("Degenerate CQUAD4 element: local X axis cannot be defined.")
        vx = vx / vx_norm
        
        # Y-axis: Completes the orthogonal system
        vy = np.cross(vz, vx)
        vy = vy / np.linalg.norm(vy)
        
        self.local_x_axis, self.local_y_axis, self.local_z_axis = vx, vy, vz

        projected_z = (self.coords - self.origin) @ vz
        characteristic_length = max(
            np.linalg.norm(self.coords[(i + 1) % 4] - self.coords[i]) for i in range(4)
        )
        if np.ptp(projected_z) > 0.2 * characteristic_length:
            raise ValueError("CQUAD4 element is excessively warped.")
        
        R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)
        self.T = CoordinateTransformation.build_full_transformation_matrix(
            R, nodes_count=4, dofs_per_node=6
        )

    def _get_shape_functions(self, xi: float, eta: float):
        """
        @brief Returns Bilinear Shape Functions and their natural derivatives.
        @return Tuple (N, dN_dxi, dN_deta) arrays of size (4,)
        """
        N = np.array([
            0.25 * (1.0 - xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 - eta),
            0.25 * (1.0 + xi) * (1.0 + eta),
            0.25 * (1.0 - xi) * (1.0 + eta)
        ])
        
        dN_dxi = np.array([
            -0.25 * (1.0 - eta),
             0.25 * (1.0 - eta),
             0.25 * (1.0 + eta),
            -0.25 * (1.0 + eta)
        ])
        
        dN_deta = np.array([
            -0.25 * (1.0 - xi),
            -0.25 * (1.0 + xi),
             0.25 * (1.0 + xi),
             0.25 * (1.0 - xi)
        ])
        
        return N, dN_dxi, dN_deta

    def _get_jacobian(self, dN_dxi: np.ndarray, dN_deta: np.ndarray):
        """
        @brief Computes the 2x2 Jacobian matrix and its determinant at an integration point.
        """
        dx_dxi = np.dot(dN_dxi, self.x_loc)
        dy_dxi = np.dot(dN_dxi, self.y_loc)
        dx_deta = np.dot(dN_deta, self.x_loc)
        dy_deta = np.dot(dN_deta, self.y_loc)
        
        J = np.array([
            [dx_dxi, dy_dxi],
            [dx_deta, dy_deta]
        ])
        
        detJ = np.linalg.det(J)
        if detJ <= 0.0:
            raise ValueError(f"Negative or zero Jacobian in CQUAD4: Nodes {self.nodes[0].nid}-{self.nodes[3].nid}. Element may be severely distorted or nodes not numbered consecutively.")
            
        return J, detJ

    def _get_B_matrices(self, xi: float, eta: float):
        """
        @brief Computes the Strain-Displacement B-matrices at an integration point.
        @return Tuple (B_m, B_b, B_s, detJ)
        """
        N, dN_dxi, dN_deta = self._get_shape_functions(xi, eta)
        J, detJ = self._get_jacobian(dN_dxi, dN_deta)
        invJ = np.linalg.inv(J)
        
        # Derivatives with respect to physical x and y
        dN_dxy = invJ @ np.vstack((dN_dxi, dN_deta))
        dN_dx = dN_dxy[0, :]
        dN_dy = dN_dxy[1, :]
        
        # 1. B_membrane (3 x 8)
        B_m = np.zeros((3, 8))
        for i in range(4):
            B_m[0, 2*i]     = dN_dx[i]
            B_m[1, 2*i + 1] = dN_dy[i]
            B_m[2, 2*i]     = dN_dy[i]
            B_m[2, 2*i + 1] = dN_dx[i]
            
        # 2. B_bending (3 x 12)
        # RM DOFs mapped internally as: [w, th_x, th_y] per node
        B_b = np.zeros((3, 12))
        for i in range(4):
            # kappa_x = d(th_y)/dx
            # kappa_y = -d(th_x)/dy
            # kappa_xy = d(th_y)/dy - d(th_x)/dx
            B_b[0, 3*i + 2] = dN_dx[i]   # kappa_x depende de th_y
            B_b[1, 3*i + 1] = -dN_dy[i]
            B_b[2, 3*i + 1] = -dN_dx[i]
            B_b[2, 3*i + 2] = dN_dy[i]   # kappa_xy depende de th_y
            
        # 3. B_shear (2 x 12)
        B_s = np.zeros((2, 12))
        for i in range(4):
            # gamma_xz = dw/dx + th_y
            # gamma_yz = dw/dy - th_x
            B_s[0, 3*i]     = dN_dx[i]
            B_s[0, 3*i + 2] = N[i]
            B_s[1, 3*i]     = dN_dy[i]
            B_s[1, 3*i + 1] = -N[i]
            
        return B_m, B_b, B_s, detJ

    def _get_D_membrane(self) -> np.ndarray:
        """
        @brief Pure Constitutive Matrix for Plane Stress (Membrane - MID1)
        """
        E, nu = self.mat_membrane.E, self.mat_membrane.nu
        factor = E / (1.0 - nu**2)
        return factor * np.array([
            [1.0,  nu,  0.0],
            [nu,  1.0,  0.0],
            [0.0, 0.0, (1.0 - nu)/2.0]
        ])

    def _get_D_bending(self) -> np.ndarray:
        """
        @brief Pure Constitutive Matrix for Plane Stress (Bending - MID2)
        """
        E, nu = self.mat_bending.E, self.mat_bending.nu
        factor = E / (1.0 - nu**2)
        return factor * np.array([
            [1.0,  nu,  0.0],
            [nu,  1.0,  0.0],
            [0.0, 0.0, (1.0 - nu)/2.0]
        ])

    def _get_D_shear(self) -> np.ndarray:
        """
        @brief Pure Constitutive Matrix for Transverse Shear (MID3)
        """
        if self.mat_shear is None:
            return np.zeros((2, 2))
        E, nu = self.mat_shear.E, self.mat_shear.nu
        G = getattr(self.mat_shear, 'G', E / (2.0 * (1.0 + nu)))
        return G * np.eye(2)

    def _get_bending_inertia_ratio(self) -> float:
        """Return PSHELL 12I/T^3, including the documented default."""
        return float(getattr(self.prop, "bending_inertia_ratio", 1.0))

    def _get_shear_thickness_ratio(self) -> float:
        """Return PSHELL TS/T, including the documented default."""
        return float(getattr(self.prop, "shear_thickness_ratio", 5.0 / 6.0))

    def _add_drilling_stiffness(self, stiffness: np.ndarray) -> None:
        """Add K6ROT with the same nodal penalty used by CTRIA3.

        CTRIA3 penalizes the difference between every nodal drilling rotation
        and the element macro-rotation ``omega = 0.5 * (dv/dx - du/dy)``.
        For CQUAD4 there is one residual equation per corner node. Its
        macro-rotation and Jacobian are evaluated at that corner, which is
        the four-node counterpart of CTRIA3's constant derivatives/Jacobian.
        """
        if self.k6rot <= 0.0:
            return

        shear_modulus = getattr(
            self.mat_membrane,
            "G",
            self.mat_membrane.E / (2.0 * (1.0 + self.mat_membrane.nu)),
        )
        penalty_factor = (
            self.k6rot
            * 1.0e-6
            * shear_modulus
            * self.prop.t
        )
        drill_local = np.zeros((12, 12))
        natural_corners = (
            (-1.0, -1.0),
            (1.0, -1.0),
            (1.0, 1.0),
            (-1.0, 1.0),
        )

        for residual_node, (xi, eta) in enumerate(natural_corners):
            _, dN_dxi, dN_deta = self._get_shape_functions(xi, eta)
            jacobian, det_jacobian = self._get_jacobian(dN_dxi, dN_deta)
            dN_dx, dN_dy = np.linalg.inv(jacobian) @ np.vstack(
                (dN_dxi, dN_deta)
            )

            # [u1, v1, theta_z1, ..., u4, v4, theta_z4]. The row is
            # theta_zi - omega at corner i, as in the CTRIA3 penalty.
            residual = np.zeros(12)
            for node_index in range(4):
                base = 3 * node_index
                residual[base] = 0.5 * dN_dy[node_index]
                residual[base + 1] = -0.5 * dN_dx[node_index]
            residual[3 * residual_node + 2] = 1.0

            # The CQUAD natural domain has area four. Nastran therefore
            # assigns each nodal penalty 4*detJ evaluated at its corner.
            nodal_area = 4.0 * det_jacobian
            drill_local += (
                penalty_factor * nodal_area * np.outer(residual, residual)
            )

        drill_dofs = [0, 1, 5, 6, 7, 11, 12, 13, 17, 18, 19, 23]
        stiffness[np.ix_(drill_dofs, drill_dofs)] += drill_local

    @staticmethod
    def _scatter_block(
        stiffness: np.ndarray, block: np.ndarray, dofs: list[int]
    ) -> None:
        stiffness[np.ix_(dofs, dofs)] += block

    @abstractmethod
    def _build_local_stiffness(self) -> np.ndarray:
        """Build the concrete formulation's 24x24 local stiffness matrix."""

    def get_stiffness_matrix(self) -> np.ndarray:
        """
        @brief Computes K_global = T.T * K_local * T
        """
        K_loc = self._build_local_stiffness()
        stiffness = self.T.T @ K_loc @ self.T
        return 0.5 * (stiffness + stiffness.T)

    def _get_recovery_displacements(
        self, element_displacements_global: np.ndarray
    ) -> np.ndarray:
        """Transform element displacements to the formulation's local frame."""
        return self.T @ element_displacements_global

    def calculate_stresses(self, u_global_vector: np.ndarray, dof_mapper) -> dict:
        """
        @brief Calculates Membrane Stresses, Moments and Shear Forces at element center.
        """
        # 1. Recuperar Desplazamientos Globales
        idxs = []
        for n in self.nodes:
            idxs.extend(dof_mapper.get_dof_indices(n.nid))
            
        u_elem_global = u_global_vector[idxs]
        
        # 2. Desplazamientos Locales
        u_loc = self._get_recovery_displacements(u_elem_global)
        
        u_mem = np.array([
            u_loc[0], u_loc[1], u_loc[6], u_loc[7],
            u_loc[12], u_loc[13], u_loc[18], u_loc[19]
        ])
        
        # Evaluamos las tensiones en el CENTRO del elemento cuadrilátero
        pts_1, _ = GaussQuadrature.get_rule(1)
        xi_c, eta_c = pts_1[0], pts_1[0]
        
        B_m, B_b, B_s, _ = self._get_B_matrices(xi_c, eta_c)
        
        t = self.prop.t
        
        # 3. Tensiones de Membrana
        strains = B_m @ u_mem
        D_mem = self._get_D_membrane()
        stresses = D_mem @ strains 
        
        sig_x, sig_y, tau_xy = stresses
        angle = np.radians(self.theta)
        cosine, sine = np.cos(angle), np.sin(angle)

        def rotate_symmetric(x_value, y_value, xy_value):
            return (
                cosine**2*x_value + sine**2*y_value + 2*cosine*sine*xy_value,
                sine**2*x_value + cosine**2*y_value - 2*cosine*sine*xy_value,
                -cosine*sine*x_value + cosine*sine*y_value
                + (cosine**2-sine**2)*xy_value,
            )

        sig_x, sig_y, tau_xy = rotate_symmetric(sig_x, sig_y, tau_xy)
        vm_stress = np.sqrt(sig_x**2 - sig_x*sig_y + sig_y**2 + 3*tau_xy**2)
        
        # Tensiones Principales
        sig_avg = (sig_x + sig_y) / 2.0
        sig_diff = (sig_x - sig_y) / 2.0
        R = np.sqrt(sig_diff**2 + tau_xy**2)
        s_major = sig_avg + R
        s_minor = sig_avg - R
        
        # Fuerzas en el plano (Fuerza = Tensión * Espesor)
        f_x = sig_x * t
        f_y = sig_y * t
        f_xy = tau_xy * t
        
        # 4. Esfuerzos de Flexión y Cortadura (Reissner-Mindlin)
        has_bending = t > 0.0
        if has_bending:
            u_bend_nastran = np.array([
                u_loc[2], u_loc[3], u_loc[4],
                u_loc[8], u_loc[9], u_loc[10],
                u_loc[14], u_loc[15], u_loc[16],
                u_loc[20], u_loc[21], u_loc[22]
            ])
            
            u_bend_RM = u_bend_nastran
            
            # Momentos Flectores: M = D_b * (t^3 / 12) * Curvatura
            kappa = B_b @ u_bend_RM
            D_b = self._get_D_bending()
            I_bend = (
                (t**3) / 12.0 * self._get_bending_inertia_ratio()
            )
            moments = (D_b * I_bend) @ kappa

            # Nastran's positive shell moment convention is opposite to the
            # curvature convention used by B_b.
            m_x, m_y, m_xy = -moments[0], -moments[1], -moments[2]
            m_x, m_y, m_xy = rotate_symmetric(m_x, m_y, m_xy)
            
            if hasattr(self, "_recover_transverse_shear_forces"):
                shear_forces = self._recover_transverse_shear_forces(
                    u_bend_RM
                )
            else:
                # Classical Reissner-Mindlin recovery.
                gamma_s = B_s @ u_bend_RM
                D_s = self._get_D_shear()
                A_shear = self._get_shear_thickness_ratio() * t
                shear_forces = (D_s * A_shear) @ gamma_s
            
            q_x, q_y = shear_forces[0], shear_forces[1]
            q_x, q_y = cosine*q_x + sine*q_y, -sine*q_x + cosine*q_y
        else:
            m_x, m_y, m_xy, q_x, q_y = 0.0, 0.0, 0.0, 0.0, 0.0
            
        # 5. Combinación de Tensiones Top / Bottom
        # sigma_bend = (M * (t/2)) / I_bend = M * (t/2) / (t^3/12) = 6M / t^2
        sb_x = 6.0 * m_x / t**2 if t > 0 else 0.0
        sb_y = 6.0 * m_y / t**2 if t > 0 else 0.0
        sb_xy = 6.0 * m_xy / t**2 if t > 0 else 0.0
        
        # Z1 is the negative fibre (first F06 row), Z2 the positive fibre.
        sig_x_bot = sig_x + sb_x
        sig_y_bot = sig_y + sb_y
        tau_xy_bot = tau_xy + sb_xy

        sig_x_top = sig_x - sb_x
        sig_y_top = sig_y - sb_y
        tau_xy_top = tau_xy - sb_xy

        def principal_and_vm(sx, sy, txy):
            average = 0.5 * (sx + sy)
            radius = np.hypot(0.5 * (sx - sy), txy)
            von_mises = np.sqrt(sx**2 - sx*sy + sy**2 + 3.0*txy**2)
            return average + radius, average - radius, von_mises

        major_top, minor_top, vm_top = principal_and_vm(
            sig_x_top, sig_y_top, tau_xy_top
        )
        major_bot, minor_bot, vm_bot = principal_and_vm(
            sig_x_bot, sig_y_bot, tau_xy_bot
        )
        
        return {
            'sigma_x': sig_x,
            'sigma_y': sig_y,
            'tau_xy': tau_xy,
            'von_mises': vm_stress,
            'major_principal': s_major,
            'minor_principal': s_minor,
            
            'bending_moment_x': m_x,
            'bending_moment_y': m_y,
            'bending_moment_xy': m_xy,
            'stress_fiber_top_x': sig_x_top,
            'stress_fiber_bot_x': sig_x_bot,
            'stress_fiber_top_y': sig_y_top,
            'stress_fiber_bot_y': sig_y_bot,
            'tau_fiber_top_y': tau_xy_top,
            'tau_fiber_bot_y': tau_xy_bot,
            'stress_top_fiber_x': sig_x_top,
            'stress_top_fiber_y': sig_y_top,
            'tau_top_fiber_xy': tau_xy_top,
            'stress_z1_fiber_x': sig_x_bot,
            'stress_z1_fiber_y': sig_y_bot,
            'tau_z1_fiber_xy': tau_xy_bot,
            'stress_bot_fiber_x': sig_x_bot,
            'stress_bot_fiber_y': sig_y_bot,
            'tau_bot_fiber_xy': tau_xy_bot,
            'stress_z2_fiber_x': sig_x_top,
            'stress_z2_fiber_y': sig_y_top,
            'tau_z2_fiber_xy': tau_xy_top,
            'major_principal_top': major_top,
            'minor_principal_top': minor_top,
            'von_mises_top': vm_top,
            'major_principal_bot': major_bot,
            'minor_principal_bot': minor_bot,
            'von_mises_bot': vm_bot,
            
            'membrane_force_x': f_x,
            'membrane_force_y': f_y,
            'membrane_force_xy': f_xy,
            'shear_force_x': q_x,
            'shear_force_y': q_y
        }

