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
    
    def __init__(self, node1: Node, node2: Node, node3: Node, prop: ShellProperty, mat: IsotropicMaterial, mat_bending: IsotropicMaterial = None, mat_shear: IsotropicMaterial = None, k6rot: float = 0.0):
        """
        @brief Initialize the CTRIA3 Element Solver.
        @param mat_bending Optional material for bending stiffness (if different from membrane).
        """
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.prop = prop
        self.k6rot = k6rot
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

        
        '''
        INTENTO SIMPLE:
        # ---------------------------------------------------------------------
        # 3. Estabilización Drilling Dinámica (Sintonizada con MSC Nastran)
        # ---------------------------------------------------------------------
        
        if self.k6rot > 0.0:
            # Escalamos tu fórmula original usando el factor exacto deducido de Nastran
            # El factor 0.34417 ajusta tu constante al valor real por unidad de K6ROT
            k_drilling = self.k6rot * 0.34417 * self.mat.E * (self.prop.t**3)
        else:
            # Si K6ROT es 0.0, dejamos la micro-estabilización residual para evitar singularidades
            k_drilling = 1e-6 * self.mat.E * (self.prop.t**3)
            
        K_local[5, 5] += k_drilling
        K_local[11, 11] += k_drilling
        K_local[17, 17] += k_drilling
        '''
        # ---------------------------------------------------------------------
        # INTENTO COMPLICADO
        # '''''''SOLUCIÓN ANALÍTICA EXACTA (FORMULACIÓN DE NASTRAN)''''''
        # 3. MATRIZ DE PENALIZACIÓN DE DRILLING CONSISTENTE (PARAM, K6ROT)
        # ---------------------------------------------------------------------
        # Grados de libertad implicados localmente en el plano (9 DOFs):
        # [u1, v1, tz1, u2, v2, tz2, u3, v3, tz3]
        
        K_drill_9x9 = np.zeros((9, 9))
        A = self.A
        
        if self.k6rot > 0.0:
            # Constante de penalización exacta calibrada de MSC Nastran
            G = self.mat.G
            t = self.prop.t
            gamma = 1e-4 * self.k6rot * G * t
            
            # Constantes analíticas de la geometría del triángulo local
            b1 = -self.y3
            b2 = self.y3
            b3 = 0.0
            
            c1 = self.x3 - self.x2
            c2 = -self.x3
            c3 = self.x2
            
            # Componentes del vector de vorticidad mecánica de membrana w_z
            # w_z = (1/4A) * sum(-c_j * u_j + b_j * v_j)
            w_u1 = -c1 / (4.0 * A)
            w_v1 =  b1 / (4.0 * A)
            w_u2 = -c2 / (4.0 * A)
            w_v2 =  b2 / (4.0 * A)
            w_u3 = -c3 / (4.0 * A)
            w_v3 =  b3 / (4.0 * A)
            
            # Nastran define tres vectores de restricción independientes (uno por nodo)
            # para evitar acoplamientos espurios entre theta_zi y theta_zj
            h1 = np.array([-w_u1, -w_v1, 1.0, -w_u2, -w_v2, 0.0, -w_u3, -w_v3, 0.0])
            h2 = np.array([-w_u1, -w_v1, 0.0, -w_u2, -w_v2, 1.0, -w_u3, -w_v3, 0.0])
            h3 = np.array([-w_u1, -w_v1, 0.0, -w_u2, -w_v2, 0.0, -w_u3, -w_v3, 1.0])
            
            # La rigidez es la suma del producto exterior de las restricciones nodales
            K_drill_9x9 = gamma * (np.outer(h1, h1) + np.outer(h2, h2) + np.outer(h3, h3))
        else:
            # Si K6ROT = 0.0, la rigidez de drilling es nula según el estándar analítico
            k_residual = 0.0
            for idx in [2, 5, 8]:
                K_drill_9x9[idx, idx] = k_residual

        # Ensamble de la matriz de drilling de 9x9 dentro de la K_local de 18x18
        map_18 = [0, 1, 5, 6, 7, 11, 12, 13, 17]
        for i in range(9):
            for j in range(9):
                K_local[map_18[i], map_18[j]] += K_drill_9x9[i, j]
        

        
        '''MIRAR SI ESTO TIENE SENTIDO O ES METER MÁS RUIDO
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
        '''
        
        
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
