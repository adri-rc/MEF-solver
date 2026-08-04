"""
@file cquad4_element.py
@brief Finite Element Solver for the CQUAD4 (Quadrilateral Shell) element.
@details
Implements the stiffness matrix calculation for a 4-node Quadrilateral Shell.
        - Membrane Stiffness: Bilinear Isoparametric Formulation (2x2 Gauss Quadrature, or SRI for shear components).
        - Bending Stiffness: Reissner-Mindlin Formulation (2x2 Gauss Quadrature).
        - Transverse Shear: Selective Reduced Integration (1x1 Gauss Quadrature) to avoid Shear Locking."""

import numpy as np
from src.preprocessing.data_model import Node, ShellProperty, IsotropicMaterial

from ..quadrature import GaussQuadrature
from ..transformations import CoordinateTransformation

class CQUAD4ElementSolver:
    """
    @class CQUAD4ElementSolver
    @brief Solver logic for a 4-node Quadrilateral Shell element.
    @details
    Computes the Global Stiffness Matrix (24x24) for a quadrilateral shell.
    
    Supported DOFs per node: 6 (Tx, Ty, Tz, Rx, Ry, Rz).
    Total DOFs: 24.
    
    Local Coordinate System (Nastran standard):
      z: Normal to the plane defined by the two diagonals (N1->N3 x N2->N4).
      x: Vector from the midpoint of edge 1-4 to the midpoint of edge 2-3.
      y: Perpendicular in plane (z cross x).
      Origin: Mean of the 4 nodes.
    
    Physics:
      - Membrane (In-plane): u (tx), v (ty). Bilinear isoparametric formulation (2x2 integration).
      - Bending (Out-of-plane): w, theta_x, theta_y. Reissner-Mindlin (2x2 integration).
      - Transverse Shear: Selective Reduced Integration (1x1 integration) to avoid shear locking.
      - Drilling: theta_z. Uses artificial stiffness.
    """

    # --- CONSTANTS FOR DOF MAPPING ---
    # Nodes: 1, 2, 3, 4
    # DOFs per node: 0-5
    # Global Local Vector indices: 
    # Node 1: 0-5, Node 2: 6-11, Node 3: 12-17, Node 4: 18-23
    
    def __init__(self, node1: Node, node2: Node, node3: Node, node4: Node, 
                 prop: ShellProperty, mat_membrane: IsotropicMaterial, 
                 mat_bending: IsotropicMaterial = None, mat_shear: IsotropicMaterial = None,
                 use_membrane_SRI: bool = False):
        """
        @brief Initialize the CQUAD4 Element Solver.
        @param mat_bending Optional material for bending stiffness (if different from membrane).
        @param use_membrane_SRI Toggle Selective Reduced Integration for In-Plane Shear.
        """
        self.nodes = [node1, node2, node3, node4]
        self.prop = prop
        self.use_membrane_SRI = use_membrane_SRI
        
        # Asignación explícita de materiales
        self.mat_membrane = mat_membrane
        self.mat_bending = mat_bending if mat_bending is not None else self.mat_membrane
        self.mat_shear = mat_shear if mat_shear is not None else self.mat_bending
        
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
        
        # X-axis: Connects midpoint of side 1-4 to midpoint of side 2-3
        mid_14 = 0.5 * (self.coords[0] + self.coords[3])
        mid_23 = 0.5 * (self.coords[1] + self.coords[2])
        vx = mid_23 - mid_14
        vx = vx / np.linalg.norm(vx)
        
        # Y-axis: Completes the orthogonal system
        vy = np.cross(vz, vx)
        vy = vy / np.linalg.norm(vy)
        
        self.local_x_axis, self.local_y_axis, self.local_z_axis = vx, vy, vz
        
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
            B_b[1, 3*i + 1] = dN_dy[i]  #antes -dN_dy[i]  # kappa_y depende de -th_x
            B_b[2, 3*i + 1] = dN_dx[i]  #antes -dN_dx[i]  # kappa_xy depende de -th_x
            B_b[2, 3*i + 2] = dN_dy[i]   # kappa_xy depende de th_y
            
        # 3. B_shear (2 x 12)
        B_s = np.zeros((2, 12))
        for i in range(4):
            # gamma_xz = dw/dx + th_y
            # gamma_yz = dw/dy - th_x
            B_s[0, 3*i]     = dN_dx[i]
            B_s[0, 3*i + 2] = N[i]
            B_s[1, 3*i]     = dN_dy[i]
            B_s[1, 3*i + 1] = N[i] #Antes: -N[i]
            
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
        E, nu = self.mat_shear.E, self.mat_shear.nu
        G = E / (2.0 * (1.0 + nu))
        return G * np.eye(2)

    '''
    def _build_local_stiffness_old(self) -> np.ndarray:
        """
        @brief Computes Element Stiffness in Local Coordinates (24x24).
        """
        K_local = np.zeros((24, 24))
        
        # --- GAUSS QUADRATURE INTEGRATION RULES ---
        pts_2, weights_2 = GaussQuadrature.get_rule(2) # 2x2 para Membrana y Flexión
        pts_1, weights_1 = GaussQuadrature.get_rule(1) # 1x1 para Cortante
        
        K_mem = np.zeros((8, 8))
        K_bend = np.zeros((12, 12))
        K_shear = np.zeros((12, 12))
        
        D_m = self._get_D_membrane()
        D_b = self._get_D_bending()
        D_s = self._get_D_shear()
        
        t = self.prop.t
        has_bending = t > 0.0
        
        # Propiedades seccionales separadas de la ley constitutiva
        I_bend = (t**3) / 12.0
        k_shear = 5.0 / 6.0 
        A_shear = k_shear * t
        
        # --- 1. MEMBRANE & BENDING (2x2 Full Integration) ---
        for i in range(2):
            for j in range(2):
                xi = pts_2[i]
                eta = pts_2[j]
                w = weights_2[i] * weights_2[j]
                
                B_m, B_b, _, detJ = self._get_B_matrices(xi, eta)
                
                # Elemento Diferencial de Área
                dA = detJ * w
                
                # Ensamblaje Membrana: K_m = (B^T * D_m * B) * (Area * t)
                K_mem += (B_m.T @ D_m @ B_m) * (dA * t)
                
                # Ensamblaje Flexión: K_b = (B^T * D_b * B) * (Area * I_bend)
                if has_bending:
                    K_bend += (B_b.T @ D_b @ B_b) * (dA * I_bend)
                    
        # --- 2. TRANSVERSE SHEAR (1x1 Selective Reduced Integration) ---
        if has_bending:
            xi_s, eta_s = pts_1[0], pts_1[0]
            w_s = weights_1[0] * weights_1[0] 
            
            _, _, B_s, detJ_s = self._get_B_matrices(xi_s, eta_s)
            dA_s = detJ_s * w_s
            
            # Ensamblaje Cortante: K_s = (B^T * D_s * B) * (Area * A_shear)
            K_shear += (B_s.T @ D_s @ B_s) * (dA_s * A_shear)

        # --- 3. SCATTER INTO LOCAL 24x24 MATRIX ---
        # Membrane: DOFs 0 (Tx), 1 (Ty) per node
        mem_idx = [0, 1, 6, 7, 12, 13, 18, 19]
        for i, r in enumerate(mem_idx):
            for j, c in enumerate(mem_idx):
                K_local[r, c] += K_mem[i, j]
                
        # Bending & Shear
        if has_bending:
            K_plate_RM = K_bend + K_shear
            
            # Map RM DOFs [w, th_x, th_y] to Nastran DOFs [Tz, Rx, Ry]
            T_node = np.array([
                [1.0,  0.0,  0.0],
                [0.0,  0.0,  1.0],
                [0.0, -1.0,  0.0]
            ])
            T_map = np.kron(np.eye(4), T_node)
            K_plate_nastran = T_map.T @ K_plate_RM @ T_map
            
            # Indices for Bending DOFs in 24x24 matrix
            bend_idx = [2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22]
            for i, r in enumerate(bend_idx):
                for j, c in enumerate(bend_idx):
                    K_local[r, c] += K_plate_nastran[i, j]
                    
        # --- 4. STABILIZATION (Drilling DOFs) ---
        k_diag_small = 0.3845 * self.mat_membrane.E * (t**3)
        drill_idx = [5, 11, 17, 23]
        for idx in drill_idx:
            K_local[idx, idx] += k_diag_small
            
        return K_local
        '''
    
    def _build_local_stiffness(self) -> np.ndarray:
        K_local = np.zeros((24, 24))
        
        pts_2, weights_2 = GaussQuadrature.get_rule(2)
        pts_1, weights_1 = GaussQuadrature.get_rule(1)
        
        K_mem = np.zeros((8, 8))
        K_bend = np.zeros((12, 12))
        K_shear = np.zeros((12, 12))
        
        D_m = self._get_D_membrane()
        D_b = self._get_D_bending()
        D_s = self._get_D_shear()
        
        t = self.prop.t
        has_bending = t > 0.0
        
        I_bend = (t**3) / 12.0
        k_shear = 5.0 / 6.0 
        A_shear = k_shear * t
        
        # Pre-allocate para Integración Reducida Selectiva de Membrana
        K_mem_axial = np.zeros((8, 8))
        if self.use_membrane_SRI:
            # Separamos las matrices D en Axial (2x2) y Shear (1x1)
            D_m_axial = D_m[0:2, 0:2]
            D_m_shear = np.array([[D_m[2, 2]]]) # Aseguramos forma 2D
        
        # --- 1. FULL INTEGRATION LOOP (2x2) ---
        for i in range(2):
            for j in range(2):
                xi, eta = pts_2[i], pts_2[j]
                w = weights_2[i] * weights_2[j]
                
                B_m, B_b, _, detJ = self._get_B_matrices(xi, eta)
                dA = detJ * w
                
                # Ensamblaje Membrana
                if self.use_membrane_SRI:
                    # SRI: Integramos con 2x2 SOLO la parte axial (u, v tracciones)
                    B_m_axial = B_m[0:2, :]
                    K_mem_axial += (B_m_axial.T @ D_m_axial @ B_m_axial) * (dA * t)
                else:
                    # Integración Completa clásica
                    K_mem += (B_m.T @ D_m @ B_m) * (dA * t)
                
                # Ensamblaje Flexión (Reissner-Mindlin bending)
                if has_bending:
                    K_bend += (B_b.T @ D_b @ B_b) * (dA * I_bend)
                    
        # --- 2. SELECTIVE REDUCED INTEGRATION (1x1) ---
        # A) Membrana: Evaluamos el cortante parasito en el centro
        if self.use_membrane_SRI:
            xi_c, eta_c = pts_1[0], pts_1[0]
            w_c = weights_1[0] * weights_1[0] # 2.0 * 2.0 = 4.0
            
            B_m_c, _, _, detJ_c = self._get_B_matrices(xi_c, eta_c)
            B_m_shear = B_m_c[2:3, :] # Extraemos la 3ª fila (cortante en el plano)
            dA_c = detJ_c * w_c
            
            K_mem_shear = (B_m_shear.T @ D_m_shear @ B_m_shear) * (dA_c * t)
            K_mem = K_mem_axial + K_mem_shear

        # B) Flexión: Evaluamos el cortante transversal (Evita out-of-plane shear locking)
        if has_bending:
            xi_s, eta_s = pts_1[0], pts_1[0]
            w_s = weights_1[0] * weights_1[0] 
            
            _, _, B_s, detJ_s = self._get_B_matrices(xi_s, eta_s)
            dA_s = detJ_s * w_s
            K_shear += (B_s.T @ D_s @ B_s) * (dA_s * A_shear)

        # --- 3. SCATTER INTO LOCAL 24x24 MATRIX ---
        mem_idx = [0, 1, 6, 7, 12, 13, 18, 19]
        for i, r in enumerate(mem_idx):
            for j, c in enumerate(mem_idx):
                K_local[r, c] += K_mem[i, j]
                
        if has_bending:
            K_plate_RM = K_bend + K_shear
            T_node = np.array([
                [1.0,  0.0,  0.0],
                [0.0,  0.0,  -1.0], # Rx = ty  -> th_y = Rx -> No, Nastran Ry = -th_x
                [0.0,  1.0,  0.0]   # Nastran Rx = th_y
            ])
            T_map = np.kron(np.eye(4), T_node)
            K_plate_nastran = T_map.T @ K_plate_RM @ T_map
            
            bend_idx = [2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22]
            for i, r in enumerate(bend_idx):
                for j, c in enumerate(bend_idx):
                    K_local[r, c] += K_plate_nastran[i, j]
                    
        '''
        # --- 4. STABILIZATION (Drilling DOFs) ---
        k_diag_small = 0.3845 * self.mat_membrane.E * (t**3)
        drill_idx = [5, 11, 17, 23]
        for idx in drill_idx:
            K_local[idx, idx] += k_diag_small
        '''
            
        return K_local

    def get_stiffness_matrix(self) -> np.ndarray:
        """
        @brief Computes K_global = T.T * K_local * T
        """
        K_loc = self._build_local_stiffness()
        return self.T.T @ K_loc @ self.T

    '''
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
        u_loc = self.T @ u_elem_global
        
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
            
            T_node = np.array([
                [1.0,  0.0,  0.0],
                [0.0,  0.0,  -1.0],
                [0.0, 1.0,  0.0]
            ])
            T_map = np.kron(np.eye(4), T_node)
            u_bend_RM = T_map @ u_bend_nastran
            
            # Momentos Flectores: M = D_b * (t^3 / 12) * Curvatura
            kappa = B_b @ u_bend_RM
            D_b = self._get_D_bending()
            I_bend = (t**3) / 12.0
            moments = (D_b * I_bend) @ kappa
            
            m_x, m_y, m_xy = moments[0], moments[1], moments[2]
            
            # Fuerzas Cortantes Transversales: Q = D_s * (k * t) * Deformacion_Angular
            gamma_s = B_s @ u_bend_RM
            D_s = self._get_D_shear()
            k_shear = 5.0 / 6.0
            A_shear = k_shear * t
            shear_forces = (D_s * A_shear) @ gamma_s
            
            q_x, q_y = shear_forces[0], shear_forces[1]
        else:
            m_x, m_y, m_xy, q_x, q_y = 0.0, 0.0, 0.0, 0.0, 0.0
            
        # 5. Combinación de Tensiones Top / Bottom
        # sigma_bend = (M * (t/2)) / I_bend = M * (t/2) / (t^3/12) = 6M / t^2
        sb_x = 6.0 * m_x / t**2 if t > 0 else 0.0
        sb_y = 6.0 * m_y / t**2 if t > 0 else 0.0
        sb_xy = 6.0 * m_xy / t**2 if t > 0 else 0.0
        
        sig_x_top = sig_x - sb_x
        sig_y_top = sig_y - sb_y
        tau_xy_top = tau_xy - sb_xy
        
        sig_x_bot = sig_x + sb_x
        sig_y_bot = sig_y + sb_y
        tau_xy_bot = tau_xy + sb_xy
        
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
            
            'membrane_force_x': f_x,
            'membrane_force_y': f_y,
            'membrane_force_xy': f_xy,
            'shear_force_x': q_x,
            'shear_force_y': q_y
        }
    '''

