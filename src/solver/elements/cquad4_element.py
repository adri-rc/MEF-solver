"""Shared mechanics for four-node quadrilateral shell elements.

Contains geometry, coordinate transformations, constitutive matrices, and
stress recovery shared by the independent Reissner-Mindlin and Nastran
implementations.  Stiffness formulations live in their own modules."""

from abc import ABC, abstractmethod

import numpy as np

from src.preprocessing.data_model import IsotropicMaterial, Node, ShellProperty

from ..quadrature import GaussQuadrature
from ..transformations import CoordinateTransformation


class CQUAD4ElementBase(ABC):
    """Provide geometry and recovery shared by 24-DOF CQUAD4 formulations."""
    
    def __init__(self, node1: Node, node2: Node, node3: Node, node4: Node, 
                 prop: ShellProperty, mat_membrane: IsotropicMaterial, 
                 mat_bending: IsotropicMaterial = None, mat_shear: IsotropicMaterial = None,
                 k6rot: float = 0.0, theta: float = 0.0):
        """Initialize geometry and independent PSHELL material fields."""
        self.nodes = [node1, node2, node3, node4]
        self.prop = prop
        self.k6rot = float(k6rot)
        self.theta = float(theta)

        if mat_membrane is None:
            raise ValueError("CQUAD4 requires a PSHELL MID1 membrane material.")

        # PSHELL material fields are independent.
        self.mat_membrane = mat_membrane
        # A blank PSHELL MID2 is a membrane-only property.  It must not
        # silently inherit MID1, otherwise a membrane CQUAD4 receives a
        # fictitious [T3,R1,R2] stiffness that an active-DOF-only DMIG check
        # can easily hide.
        self.mat_bending = mat_bending
        # A blank PSHELL MID3 disables physical transverse-shear flexibility.
        # It must not silently inherit MID2; the Nastran QUAD4 retains only its
        # residual bending flexibility in that case.
        self.mat_shear = mat_shear
        
        self.coords = np.array([n.get_coords() for n in self.nodes])
        
        self._compute_transformation()
        
        # Express the four coplanar nodes in the local two-dimensional plane.
        self.local_coords_2d = np.zeros((4, 2))
        for i in range(4):
            v = self.coords[i] - self.origin
            self.local_coords_2d[i, 0] = np.dot(v, self.local_x_axis)
            self.local_coords_2d[i, 1] = np.dot(v, self.local_y_axis)
        self.x_loc = self.local_coords_2d[:, 0]
        self.y_loc = self.local_coords_2d[:, 1]

        # Validate the complete bilinear map independently of the active
        # property blocks.  Otherwise a membrane-only element with K6ROT=0
        # might evaluate only the centre and fail to expose an inverted
        # corner until a different property activates bending or drilling.
        for xi, eta in (
            (-1.0, -1.0),
            (1.0, -1.0),
            (1.0, 1.0),
            (-1.0, 1.0),
        ):
            _, dN_dxi, dN_deta = self._get_shape_functions(xi, eta)
            self._get_jacobian(dN_dxi, dN_deta)
            
    def _compute_transformation(self):
        """Build the 24-by-24 transformation from Nastran CQUAD4 axes."""
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
        coplanarity_tolerance = max(1.0e-10, 1.0e-8 * characteristic_length)
        if np.ptp(projected_z) > coplanarity_tolerance:
            raise ValueError(
                "Non-coplanar CQUAD4 element: the current formulation only "
                "supports planar quadrilaterals."
            )
        
        R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)
        self.T = CoordinateTransformation.build_full_transformation_matrix(
            R, nodes_count=4, dofs_per_node=6
        )

    def _get_shape_functions(self, xi: float, eta: float):
        """Return bilinear shape functions and natural derivatives."""
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
        """Return the two-dimensional Jacobian and determinant."""
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
        """Return membrane, bending, shear, and Jacobian terms."""
        N, dN_dxi, dN_deta = self._get_shape_functions(xi, eta)
        J, detJ = self._get_jacobian(dN_dxi, dN_deta)
        invJ = np.linalg.inv(J)
        
        # Derivatives with respect to physical x and y
        dN_dxy = invJ @ np.vstack((dN_dxi, dN_deta))
        dN_dx = dN_dxy[0, :]
        dN_dy = dN_dxy[1, :]
        
        # Membrane block (3 x 8).
        B_m = np.zeros((3, 8))
        for i in range(4):
            B_m[0, 2*i]     = dN_dx[i]
            B_m[1, 2*i + 1] = dN_dy[i]
            B_m[2, 2*i]     = dN_dy[i]
            B_m[2, 2*i + 1] = dN_dx[i]
            
        # Bending block (3 x 12).
        # RM DOFs mapped internally as: [w, th_x, th_y] per node
        B_b = np.zeros((3, 12))
        for i in range(4):
            # kappa_x = d(th_y)/dx
            # kappa_y = -d(th_x)/dy
            # kappa_xy = d(th_y)/dy - d(th_x)/dx
            B_b[0, 3*i + 2] = dN_dx[i]   # kappa_x depends on theta_y
            B_b[1, 3*i + 1] = -dN_dy[i]
            B_b[2, 3*i + 1] = -dN_dx[i]
            B_b[2, 3*i + 2] = dN_dy[i]   # kappa_xy depends on theta_y
            
        # Transverse-shear block (2 x 12).
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
        """Return the MID1 plane-stress constitutive matrix."""
        E, nu = self.mat_membrane.E, self.mat_membrane.nu
        factor = E / (1.0 - nu**2)
        return factor * np.array([
            [1.0,  nu,  0.0],
            [nu,  1.0,  0.0],
            [0.0, 0.0, (1.0 - nu)/2.0]
        ])

    def _get_D_bending(self) -> np.ndarray:
        """Return the MID2 plane-stress constitutive matrix."""
        if self.mat_bending is None:
            return np.zeros((3, 3))
        E, nu = self.mat_bending.E, self.mat_bending.nu
        factor = E / (1.0 - nu**2)
        return factor * np.array([
            [1.0,  nu,  0.0],
            [nu,  1.0,  0.0],
            [0.0, 0.0, (1.0 - nu)/2.0]
        ])

    def _get_D_shear(self) -> np.ndarray:
        """Return the MID3 transverse-shear constitutive matrix."""
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
        """Return the symmetric stiffness matrix in global coordinates."""
        K_loc = self._build_local_stiffness()
        stiffness = self.T.T @ K_loc @ self.T
        return 0.5 * (stiffness + stiffness.T)

    def _get_recovery_displacements(
        self, element_displacements_global: np.ndarray
    ) -> np.ndarray:
        """Transform element displacements to the formulation's local frame."""
        return self.T @ element_displacements_global

    def calculate_stresses(self, u_global_vector: np.ndarray, dof_mapper) -> dict:
        """Recover stresses and force resultants at the element center."""
        # Gather the element displacement vector in basic coordinates.
        idxs = []
        for n in self.nodes:
            idxs.extend(dof_mapper.get_dof_indices(n.nid))
            
        u_elem_global = u_global_vector[idxs]
        
        # Transform displacements to element coordinates.
        u_loc = self._get_recovery_displacements(u_elem_global)
        
        u_mem = np.array([
            u_loc[0], u_loc[1], u_loc[6], u_loc[7],
            u_loc[12], u_loc[13], u_loc[18], u_loc[19]
        ])
        
        # Recover the response at the element center.
        pts_1, _ = GaussQuadrature.get_rule(1)
        xi_c, eta_c = pts_1[0], pts_1[0]
        
        B_m, B_b, B_s, _ = self._get_B_matrices(xi_c, eta_c)
        
        t = self.prop.t
        
        # Membrane stresses and resultants.
        strains = B_m @ u_mem
        D_mem = self._get_D_membrane()
        stresses = D_mem @ strains 
        
        sig_x, sig_y, tau_xy = stresses
        # MSC Nastran reports CQUAD4 forces and stresses in the element
        # coordinate system.  THETA or MCID orients the material axes; it
        # does not rotate isotropic MAT1 output away from the element axes.
        # Keeping recovery in (x_element, y_element) is also essential when
        # adjacent elements use different THETA values.
        vm_stress = np.sqrt(sig_x**2 - sig_x*sig_y + sig_y**2 + 3*tau_xy**2)
        
        # Principal membrane stresses.
        sig_avg = (sig_x + sig_y) / 2.0
        sig_diff = (sig_x - sig_y) / 2.0
        radius = np.sqrt(sig_diff**2 + tau_xy**2)
        s_major = sig_avg + radius
        s_minor = sig_avg - radius
        
        # In-plane force resultants (stress times thickness).
        f_x = sig_x * t
        f_y = sig_y * t
        f_xy = tau_xy * t
        
        # Bending and transverse-shear response.
        has_bending = t > 0.0 and self.mat_bending is not None
        if has_bending:
            u_bend_nastran = np.array([
                u_loc[2], u_loc[3], u_loc[4],
                u_loc[8], u_loc[9], u_loc[10],
                u_loc[14], u_loc[15], u_loc[16],
                u_loc[20], u_loc[21], u_loc[22]
            ])
            
            u_bend_RM = u_bend_nastran
            
            # Bending moments: M = D_b * (t^3 / 12) * curvature.
            kappa = B_b @ u_bend_RM
            D_b = self._get_D_bending()
            I_bend = (
                (t**3) / 12.0 * self._get_bending_inertia_ratio()
            )
            moments = (D_b * I_bend) @ kappa

            # Nastran's positive shell moment convention is opposite to the
            # curvature convention used by B_b.
            m_x, m_y, m_xy = -moments[0], -moments[1], -moments[2]
            
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
        else:
            m_x, m_y, m_xy, q_x, q_y = 0.0, 0.0, 0.0, 0.0, 0.0
            
        # Combine membrane and bending stresses at the two surfaces.
        # sigma_bend = (M * (t/2)) / I_bend = M * (t/2) / (t^3/12) = 6M / t^2
        sb_x = 6.0 * m_x / t**2 if t > 0 else 0.0
        sb_y = 6.0 * m_y / t**2 if t > 0 else 0.0
        sb_xy = 6.0 * m_xy / t**2 if t > 0 else 0.0
        
        # Public convention: bottom is local z=-t/2 (Nastran Z1) and top is
        # local z=+t/2 (Nastran Z2).  The result keys deliberately use the
        # geometric top/bottom convention instead of the solver-specific Z1/Z2
        # labels.
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
            'stress_top_fiber_x': sig_x_top,
            'stress_top_fiber_y': sig_y_top,
            'tau_top_fiber_xy': tau_xy_top,
            'stress_bot_fiber_x': sig_x_bot,
            'stress_bot_fiber_y': sig_y_bot,
            'tau_bot_fiber_xy': tau_xy_bot,
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

