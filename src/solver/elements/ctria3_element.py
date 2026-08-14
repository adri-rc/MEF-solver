"""
@file ctria3_element.py
@brief Finite Element Solver for the CTRIA3 (Triangular Shell) element.
@details
Implements the stiffness matrix calculation for a 3-node Triangular Shell.
        - Membrane Stiffness: Constant Strain Triangle (CST).
        - Bending Stiffness: Reissner-Mindlin Formulation.
IMPORTANT: 
Pure shear stiffnes does not mach with Nastrnas vales.
Therefore to match Nastran results MID3 with E=0.0 value HAVE TO BE DEFINED in PSHELL card.
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
        Pure shear stiffnes does not mach with Nastrnas vales.
        Therefore to match Nastran results MID3 with E=0.0 value HAVE TO BE DEFINED in PSHELL card.
      - Drilling: theta_z. (Defined with Nastran formula, K6ROT).
    """

    # --- CONSTANTS FOR DOF MAPPING ---
    # Nodes: 1, 2, 3
    # DOFs per node: 0-5
    # Global Local Vector indices: 
    # Node 1: 0-5, Node 2: 6-11, Node 3: 12-17
    
    def __init__(self, node1: Node, node2: Node, node3: Node, prop: ShellProperty, mat: IsotropicMaterial, mat_bending: IsotropicMaterial = None, mat_shear: IsotropicMaterial = None, k6rot: float = 0.0, theta: float = 0.0):
        """
        @brief Initialize the CTRIA3 Element Solver.
        @param mat_bending Optional material for bending stiffness (if different from membrane).
        """
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.prop = prop
        self.k6rot = k6rot
        self.theta = float(theta)
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

        '''
        ESTE IF PUEDE HACER QUE SE INVIERTA LA CARA TOP DEL ELEMENTO POR LO QEU LOS RESULTADOS DE STRESS NO CUADREN CON NASTRAN:
        if self.y3 < 0:
            self.y3 = -self.y3
            self.local_y_axis = -self.local_y_axis
            self.local_z_axis = -self.local_z_axis
            # Recompute Transformation matrix with corrected axes
            R = CoordinateTransformation.get_rotation_matrix_3x3(self.local_x_axis, self.local_y_axis, self.local_z_axis)
            self.T = CoordinateTransformation.build_full_transformation_matrix(R, 3, 6)
        '''
        # Area of triangle
        self.A = 0.5 * self.x2 * self.y3
        if self.A < 1e-12:
            raise ValueError(f"Zero area triangle: Nodes {node1.nid}, {node2.nid}, {node3.nid}")

    def _compute_transformation(self):
        """
        @brief Computes T matrix (18x18) based on element geometry.
        @details 
            Ensures strict compliance with MSC Nastran Right-Hand Rule to 
            define Top/Bottom faces correctly.
        """
        v12 = self.n2_coords - self.n1_coords
        v13 = self.n3_coords - self.n1_coords

        # 1. Eje X Local: Directamente del Nodo 1 al Nodo 2
        length_12 = np.linalg.norm(v12)
        if length_12 < 1e-12:
            raise ValueError(f"Zero length edge in CTRIA3: Nodes {self.node1.nid}-{self.node2.nid}")
        vx = v12 / length_12

        # 2. Eje Z Local (Normal): Regla de la mano derecha (X cruz v13)
        # Esto define físicamente hacia dónde mira la cara superior (Top Face / Z1)
        vz_cross = np.cross(vx, v13)
        normal_norm = np.linalg.norm(vz_cross)
        if normal_norm < 1e-12:
            raise ValueError(
                f"Zero area triangle: Nodes {self.node1.nid}, {self.node2.nid}, {self.node3.nid}"
            )
        vz = vz_cross / normal_norm

        # 3. Eje Y Local: Producto cruzado para cerrar el triedro (Z cruz X)
        vy = np.cross(vz, vx)

        # Guardamos los ejes en la instancia
        self.local_x_axis = vx
        self.local_y_axis = vy
        self.local_z_axis = vz
        
        # Construimos la matriz de transformación
        R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)
        self.T = CoordinateTransformation.build_full_transformation_matrix(
            R, nodes_count=3, dofs_per_node=6
        )

        '''
        VERSIÓN ANTIGUA:
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
        '''

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
        G_c = getattr(self.mat_shear, 'G', E_s / (2.0 * (1.0 + nu_s)))
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

    def _get_B_bending_Nastran(self):
        return self._get_B_matrices_isoparametricas()[1]

    def _get_B_shear_Nastran(self, *_natural_coordinates):
        return self._get_B_matrices_isoparametricas()[2]

    def _get_D_bending(self):
        return self._get_matrices_constitutivas()[1]

    def _get_D_shear(self):
        return self._get_matrices_constitutivas()[2]

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
        
        if self.k6rot > 0.0:
            # 3.1 Calcular el módulo de cortante G
            nu = self.mat.nu  # Asegúrate de tener acceso al coeficiente de Poisson
            G = self.mat.E / (2.0 * (1.0 + nu))
            
            # 3.2 Constante de penalización exacta de Nastran
            # Nastran usa: k_theta = K6ROT * 1e-6 * G * t * detJ
            k_theta = self.k6rot * 1e-6 * G * self.prop.t * detJ
            
            # 3.3 Extraer coeficientes de la macro-rotación omega
            # omega = 0.5 * (dv/dx - du/dy)
            # En la matriz B_m estándar: B_m[0, 2i] es dNi/dx, B_m[1, 2i+1] es dNi/dy
            b_omega = np.array([
                -0.5 * B_m[1, 1],  0.5 * B_m[0, 0],  # coeficientes para u1, v1
                -0.5 * B_m[1, 3],  0.5 * B_m[0, 2],  # coeficientes para u2, v2
                -0.5 * B_m[1, 5],  0.5 * B_m[0, 4]   # coeficientes para u3, v3
            ])
            
            # 3.4 Construir matriz de acoplamiento B_drill (3 nodos x 9 GdL: u, v, theta_z)
            B_drill = np.zeros((3, 9))
            for i in range(3):
                # Restamos la macro-rotación al nodo i
                B_drill[i, 0] = -b_omega[0] # u1
                B_drill[i, 1] = -b_omega[1] # v1
                B_drill[i, 3] = -b_omega[2] # u2
                B_drill[i, 4] = -b_omega[3] # v2
                B_drill[i, 6] = -b_omega[4] # u3
                B_drill[i, 7] = -b_omega[5] # v3
                
                # Y sumamos la rotación nodal theta_zi (1.0 en su respectiva columna)
                B_drill[i, 2 + i*3] = 1.0   
                
            # 3.5 Matriz de rigidez de drilling local (9x9)
            K_drill_local = k_theta * (B_drill.T @ B_drill)
            
            # 3.6 Ensamblar en la matriz local de 18x18
            # Índices correspondientes a [u1, v1, tz1, u2, v2, tz2, u3, v3, tz3]
            idx_drill = [0, 1, 5, 6, 7, 11, 12, 13, 17]
            for i in range(9):
                for j in range(9):
                    K_local[idx_drill[i], idx_drill[j]] += K_drill_local[i, j]
        
        return K_local

    def get_stiffness_matrix(self) -> np.ndarray:
        """
        @brief Computes K_global = T.T * K_local * T
        """
        K_loc = self._build_local_stiffness()
        stiffness = self.T.T @ K_loc @ self.T
        return 0.5 * (stiffness + stiffness.T)

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
        
        # --- Esfuerzos Membrana ---
        u_mem = np.array([u_loc[0], u_loc[1], u_loc[6], u_loc[7], u_loc[12], u_loc[13]])
        strains = B_m @ u_mem 
        stresses = D_m @ strains 
        sig_x, sig_y, tau_xy = stresses

        # Nastran reports CTRIA3 force/stress components in the element
        # coordinate system. THETA defines material orientation and must not
        # rotate isotropic MAT1 recovery output.
        
        '''        
        vm_stress = np.sqrt(sig_x**2 - sig_x*sig_y + sig_y**2 + 3*tau_xy**2)
        sig_avg, sig_diff = (sig_x + sig_y) / 2.0, (sig_x - sig_y) / 2.0
        R = np.sqrt(sig_diff**2 + tau_xy**2)
        s_major, s_minor = sig_avg + R, sig_avg - R
        '''

        t = self.prop.t
        f_x, f_y, f_xy = sig_x * t, sig_y * t, tau_xy * t
        
        # --- Esfuerzos Flexión y Cortante ---
        has_bending = self.prop.t > 0.0
        if has_bending:
            u_bend_RM = np.array([u_loc[2], u_loc[3], u_loc[4], u_loc[8], u_loc[9], u_loc[10], u_loc[14], u_loc[15], u_loc[16]])
            kappa = B_f @ u_bend_RM
            moments = D_f @ kappa
            # INVERSIÓN DE SIGNOS PARA COINCIDIR CON NASTRAN
            m_x, m_y, m_xy = -moments[0], -moments[1], -moments[2]
            
            gamma_s = B_c @ u_bend_RM
            shear_forces = D_c @ gamma_s
            q_x, q_y = shear_forces[0], shear_forces[1]
        else:
            m_x, m_y, m_xy, q_x, q_y = 0.0, 0.0, 0.0, 0.0, 0.0
            
        # Public convention: top is local z=+t/2 and bottom is local z=-t/2.
        # Avoid Z1/Z2 aliases because their meaning depends on the external
        # solver's shell-result convention.
        sb_x, sb_y, sb_xy = 6.0 * m_x / t**2, 6.0 * m_y / t**2, 6.0 * m_xy / t**2
        sig_x_top, sig_x_bot = sig_x - sb_x, sig_x + sb_x
        sig_y_top, sig_y_bot = sig_y - sb_y, sig_y + sb_y
        tau_xy_top, tau_xy_bot = tau_xy - sb_xy, tau_xy + sb_xy

        # Función auxiliar para cálculo de Von Mises y Principales
        def calc_principal_vm(sx, sy, txy):
            vm = np.sqrt(sx**2 - sx*sy + sy**2 + 3*txy**2)
            s_avg, s_diff = (sx + sy) / 2.0, (sx - sy) / 2.0
            R = np.sqrt(s_diff**2 + txy**2)
            return vm, s_avg + R, s_avg - R

        # Evaluar fibras
        vm_mem, s_maj_mem, s_min_mem = calc_principal_vm(sig_x, sig_y, tau_xy)
        vm_top, s_maj_top, s_min_top = calc_principal_vm(sig_x_top, sig_y_top, tau_xy_top)
        vm_bot, s_maj_bot, s_min_bot = calc_principal_vm(sig_x_bot, sig_y_bot, tau_xy_bot)
        
        return {
            # Membrana / Midplane
            'sigma_x': sig_x, 'sigma_y': sig_y, 'tau_xy': tau_xy,
            'von_mises': vm_mem, 'major_principal': s_maj_mem, 'minor_principal': s_min_mem,

            # Fuerzas y Momentos
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
