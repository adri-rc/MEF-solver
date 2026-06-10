"""
@file ctria3_element.py
@brief Finite Element Solver for the CTRIA3 (Triangular Shell) element.
@details
Implements the stiffness matrix calculation for a 3-node Triangular Shell.
        - Membrane Stiffness: Constant Strain Triangle (CST).
        - Bending Stiffness: Reissner-Mindlin Formulation (MIN3).
"""

import numpy as np
from src.preprocessing.data_model import Node, ShellProperty, IsotropicMaterial
from ..transformations import CoordinateTransformation

class CTRIA3ElementSolver:
    """
    @class CTRIA3ElementSolver
    @brief Solver logic for a 3-node Triangle element.
    @details
    Computes the Global Stiffness Matrix (18x18) for a triangular shell.
    
    Supported DOFs per node: 6 (Tx, Ty, Tz, Rx, Ry, Rz).
    Total DOFs: 18.
    
    Local Coordinate System (defined by transformations.compute_local_axes_triangle):
      x: Node 1 -> Node 2
      z: Normal to plane
      y: Perpendicular in plane
    
    Physics:
      - Membrane (In-plane): u (tx), v (ty). Formulation: CST (Linear shape functions).
      - Bending (Out-of-plane) & Shear: w, theta_x, theta_y. Formulation: Reissner-Mindlin.
      - Drilling: theta_z. (Often needs artificial stiffness).
    """

    # --- CONSTANTS FOR DOF MAPPING ---
    # Nodes: 1, 2, 3
    # DOFs per node: 0-5
    # Global Local Vector indices: 
    # Node 1: 0-5, Node 2: 6-11, Node 3: 12-17
    
    def __init__(self, node1: Node, node2: Node, node3: Node, prop: ShellProperty, mat: IsotropicMaterial, mat_bending: IsotropicMaterial = None, mat_shear: IsotropicMaterial = None): #, k6rot: float = 0.0
        """
        @brief Initialize the CTRIA3 Element Solver.
        @param mat_bending Optional material for bending stiffness (if different from membrane).
        """
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.prop = prop
        # self.k6rot = k6rot
        self.mat = mat # MID1
        
        # Use bending material if provided, otherwise use the same as membrane
        self.mat_bending = mat_bending if mat_bending is not None else mat # MID2
        self.mat_shear = mat_shear if mat_shear is not None else self.mat_bending # MID3
        
        self.n1_coords = node1.get_coords()
        self.n2_coords = node2.get_coords()
        self.n3_coords = node3.get_coords()
        
        # 1. Compute Local System and Transformation Matrix
        self._compute_transformation()
        
        # 2. Compute Geometric Properties in Local System
        # Express nodes in local coords
        # N1 is at (0,0,0) by definition of local system
        # N2 is at (x2, 0, 0)
        # N3 is at (x3, y3, 0)
        
        v12 = self.n2_coords - self.n1_coords
        self.x2 = np.linalg.norm(v12) # Distance N1-N2
        
        # Project N3 onto local axes to get (x3, y3)
        v13 = self.n3_coords - self.n1_coords
        self.x3 = np.dot(v13, self.local_x_axis)
        self.y3 = np.dot(v13, self.local_y_axis)

        if self.y3 < 0:
            self.y3 = -self.y3
            self.local_y_axis = -self.local_y_axis
            self.local_z_axis = -self.local_z_axis
            # Recompute Transformation matrix with corrected axes
            R = CoordinateTransformation.get_rotation_matrix_3x3(self.local_x_axis, self.local_y_axis, self.local_z_axis)
            self.T = CoordinateTransformation.build_full_transformation_matrix(R, 3, 6)
        
        # Area of triangle
        self.A = 0.5 * self.x2 * self.y3
        if self.A < 1e-12:
            raise ValueError(f"Zero area triangle: Nodes {node1.nid}, {node2.nid}, {node3.nid}")

    def _compute_transformation(self):
        """@brief Computes T matrix (18x18) based on element geometry."""
        vx, vy, vz = CoordinateTransformation.compute_local_axes_triangle(
            self.n1_coords, self.n2_coords, self.n3_coords
        )
        self.local_x_axis = vx
        self.local_y_axis = vy
        self.local_z_axis = vz
        
        R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)
        self.T = CoordinateTransformation.build_full_transformation_matrix(
            R, nodes_count=3, dofs_per_node=6
        )

    def _get_matrices_constitutivas(self):
        """
        @brief Calcula las matrices D utilizando materiales independientes para cada comportamiento.
        """
        # --- COMPORTAMIENTO A MEMBRANA (MID1) ---
        E_m, nu_m = self.mat.E, self.mat.nu
        factor_m = E_m / (1.0 - nu_m**2)
        D_m = factor_m * np.array([
            [1.0,  nu_m, 0.0],
            [nu_m, 1.0,  0.0],
            [0.0,  0.0,  (1.0 - nu_m) / 2.0]
        ])

        # --- COMPORTAMIENTO A FLEXIÓN (MID2) ---
        t = self.prop.t
        E_f, nu_f = self.mat_bending.E, self.mat_bending.nu
        factor_f = (E_f * t**3) / (12.0 * (1.0 - nu_f**2))
        D_f = factor_f * np.array([
            [1.0,  nu_f, 0.0],
            [nu_f, 1.0,  0.0],
            [0.0,  0.0,  (1.0 - nu_f) / 2.0]
        ])

        # --- COMPORTAMIENTO A CORTANTE TRANSVERSAL (MID3) ---
        E_s, nu_s = self.mat_shear.E, self.mat_shear.nu
        G_c = E_s / (2.0 * (1.0 + nu_s)) # Deducción del módulo G para material isótropo
        k_shear = 5.0 / 6.0              # Factor de corrección de Mindlin
        factor_c = k_shear * G_c * t
        
        D_c = factor_c * np.array([
            [1.0, 0.0],
            [0.0, 1.0]
        ])

        return D_m, D_f, D_c
    
    def _get_B_matrices_isoparametricas(self):
        """
        @brief Construye las matrices B usando la formulación Isoparamétrica.
        """
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
        dN_dfis = np.dot(J_inv, dN_dnat)
        dN_dx, dN_dy = dN_dfis[0, :], dN_dfis[1, :]
        
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

    def _build_local_stiffness(self) -> np.ndarray:
        """
        @brief Computes Element Stiffness in Local Coordinates (18x18).
        """
        K_local = np.zeros((18, 18))
        peso_gauss = 0.5 
        
        D_m, D_f, D_c = self._get_matrices_constitutivas()
        B_m, B_f, B_c, detJ = self._get_B_matrices_isoparametricas()

        # 1. Rigidez de Membrana
        k_mem = (self.prop.t * detJ * peso_gauss) * (B_m.T @ D_m @ B_m)
        idx_m = [0, 1, 6, 7, 12, 13]
        for i in range(6):
            for j in range(6):
                K_local[idx_m[i], idx_m[j]] += k_mem[i, j]
                
        # 2. Rigidez de Flexión y Cortante Transversal
        has_bending = self.prop.t > 0.0
        if has_bending:
            k_flex = (detJ * peso_gauss) * (B_f.T @ D_f @ B_f)
            k_cort = (detJ * peso_gauss) * (B_c.T @ D_c @ B_c)
            k_plate = k_flex + k_cort
            
            idx_f = [2, 3, 4, 8, 9, 10, 14, 15, 16]
            for i in range(9):
                for j in range(9):
                    K_local[idx_f[i], idx_f[j]] += k_plate[i, j]

        
        # ---------------------------------------------------------------------
        # 3. Estabilización Drilling Dinámica (Sintonizada con MSC Nastran)
        # ---------------------------------------------------------------------
        ''' NATES DE INTENTAR METER EL K6ROT COMO PARÁMETRO EXTERNO.
        if self.k6rot > 0.0:
            # Escalamos tu fórmula original usando el factor exacto deducido de Nastran
            # El factor 0.34417 ajusta tu constante al valor real por unidad de K6ROT
            k_drilling = self.k6rot * 0.34417 * self.mat.E * (self.prop.t**3)
        else:
            # Si K6ROT es 0.0, dejamos la micro-estabilización residual para evitar singularidades
            k_drilling = 1e-6 * self.mat.E * (self.prop.t**3)
        '''
        k_drilling = 1e-6 * self.mat.E * (self.prop.t**3)    
        K_local[5, 5] += k_drilling
        K_local[11, 11] += k_drilling
        K_local[17, 17] += k_drilling
        

        
        # ---------------------------------------------------------------------
        # 4. NUEVO: Estabilización Cortante Transversal (Pseudo-AUTOSPC)
        # Evita singularidades (mecanismos) si el usuario define MID3 nulo
        # en elementos inclinados respecto al sistema global.
        # ---------------------------------------------------------------------
        if self.mat_shear.E < 1e-5:
            # Añadimos una rigidez residual minúscula (1 millón de veces menor que la membrana).
            # Suficiente para que SciPy pueda invertir la matriz sin alterar los resultados de tensión.
            k_shear_small = 1e-6 * self.mat.E * self.prop.t 
            
            # Grados de libertad de traslación Z local (w1, w2, w3) son los índices 2, 8, 14
            K_local[2, 2] += k_shear_small
            K_local[8, 8] += k_shear_small
            K_local[14, 14] += k_shear_small
        
        
        return K_local

    def get_stiffness_matrix(self) -> np.ndarray:
        """
        @brief Computes K_global = T.T * K_local * T
        """
        K_loc = self._build_local_stiffness()
        return self.T.T @ K_loc @ self.T

    def calculate_stresses(self, u_global_vector: np.ndarray, dof_mapper) -> dict:
        """
        @brief Calculates Element Stresses and Internal Forces using exact constitutive matrices.
        """
        idxs = []
        for n in [self.node1, self.node2, self.node3]:
            idxs.extend(dof_mapper.get_dof_indices(n.nid))
            
        u_elem_global = u_global_vector[idxs]
        u_loc = self.T @ u_elem_global
        
        D_m, D_f, D_c = self._get_matrices_constitutivas()
        B_m, B_f, B_c, detJ = self._get_B_matrices_isoparametricas()
        
        # Esfuerzos Membrana
        u_mem = np.array([u_loc[0], u_loc[1], u_loc[6], u_loc[7], u_loc[12], u_loc[13]])
        strains = B_m @ u_mem 
        stresses = D_m @ strains 
        sig_x, sig_y, tau_xy = stresses
        
        vm_stress = np.sqrt(sig_x**2 - sig_x*sig_y + sig_y**2 + 3*tau_xy**2)
        sig_avg, sig_diff = (sig_x + sig_y) / 2.0, (sig_x - sig_y) / 2.0
        R = np.sqrt(sig_diff**2 + tau_xy**2)
        s_major, s_minor = sig_avg + R, sig_avg - R

        t = self.prop.t
        f_x, f_y, f_xy = sig_x * t, sig_y * t, tau_xy * t
        
        # Esfuerzos Flexión y Cortante
        has_bending = self.prop.t > 0.0
        if has_bending:
            u_bend_RM = np.array([u_loc[2], u_loc[3], u_loc[4], u_loc[8], u_loc[9], u_loc[10], u_loc[14], u_loc[15], u_loc[16]])
            kappa = B_f @ u_bend_RM
            moments = D_f @ kappa
            m_x, m_y, m_xy = moments[0], moments[1], moments[2]
            
            gamma_s = B_c @ u_bend_RM
            shear_forces = D_c @ gamma_s
            q_x, q_y = shear_forces[0], shear_forces[1]
        else:
            m_x, m_y, m_xy, q_x, q_y = 0.0, 0.0, 0.0, 0.0, 0.0
            
        sb_x, sb_y, sb_xy = 6.0 * m_x / t**2, 6.0 * m_y / t**2, 6.0 * m_xy / t**2
        sig_x_top, sig_x_bot = sig_x - sb_x, sig_x + sb_x
        sig_y_top, sig_y_bot = sig_y - sb_y, sig_y + sb_y
        tau_xy_top, tau_xy_bot = tau_xy - sb_xy, tau_xy + sb_xy
        
        return {
            'sigma_x': sig_x, 
            'sigma_y': sig_y, 
            'tau_xy': tau_xy,
            'von_mises': vm_stress, 
            'major_principal': s_major, 
            'minor_principal': s_minor,

            # Additional Bending Info
            'bending_moment_x': m_x, 
            'bending_moment_y': m_y, 
            'bending_moment_xy': m_xy,
            'stress_fiber_top_x': sig_x_top, 
            'stress_fiber_bot_x': sig_x_bot,
            'stress_fiber_top_y': sig_y_top, 
            'stress_fiber_bot_y': sig_y_bot,
            'tau_fiber_top_y': tau_xy_top, 
            'tau_fiber_bot_y': tau_xy_bot,

            # Forces
            'membrane_force_x': f_x, 
            'membrane_force_y': f_y, 
            'membrane_force_xy': f_xy,
            'shear_force_x': q_x, 
            'shear_force_y': q_y
        }
