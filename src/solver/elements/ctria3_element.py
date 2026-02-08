
"""
@file ctria3_element.py
@brief Finite Element Solver for the CTRIA3 (Triangular Shell) element.
@details
Implements the stiffness matrix calculation for a 3-node Triangular Shell.
Currently implements:
- Membrane Stiffness: Constant Strain Triangle (CST).
- Bending Stiffness: Placeholder/Basic implementation (or zero with warnings).
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
      - Bending (Out-of-plane): w (tz), theta_x, theta_y. 
      - Drilling: theta_z. (Often needs artificial stiffness).
    """

    # --- CONSTANTS FOR DOF MAPPING ---
    # Nodes: 1, 2, 3
    # DOFs per node: 0-5
    # Global Local Vector indices: 
    # Node 1: 0-5, Node 2: 6-11, Node 3: 12-17
    
    def __init__(self, node1: Node, node2: Node, node3: Node, prop: ShellProperty, mat: IsotropicMaterial):
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.prop = prop
        self.mat = mat
        
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

    def _get_B_membrane(self) -> np.ndarray:
        """
        @brief Strain-Displacement matrix for Membrane (CST).
        @details
        Strains: [eps_x, eps_y, gam_xy]
        Displacements: [u1, v1, u2, v2, u3, v3]
        
        B = (1/2A) * [ ... ]
        """
        # Beta_i = y_j - y_k
        # Gamma_i = x_k - x_j
        # Cyclic permutation: 1-2-3
        
        # Local coords:
        # 1: (0,0)
        # 2: (x2, 0)
        # 3: (x3, y3)
        x1, y1 = 0.0, 0.0
        x2, y2 = self.x2, 0.0
        x3, y3 = self.x3, self.y3
        
        # i=1, j=2, k=3
        b1 = y2 - y3
        g1 = x3 - x2
        
        # i=2, j=3, k=1
        b2 = y3 - y1
        g2 = x1 - x3
        
        # i=3, j=1, k=2
        b3 = y1 - y2
        g3 = x2 - x1
        
        # B matrix (3 rows x 6 cols)
        # [ u1 v1 u2 v2 u3 v3 ]
        # Row 1 (eps_x): b_i, 0
        # Row 2 (eps_y): 0, g_i
        # Row 3 (gam_xy): g_i, b_i
        
        B = np.zeros((3, 6))
        
        # Node 1
        B[0, 0] = b1
        B[1, 1] = g1
        B[2, 0] = g1
        B[2, 1] = b1
        
        # Node 2
        B[0, 2] = b2
        B[1, 3] = g2
        B[2, 2] = g2
        B[2, 3] = b2
        
        # Node 3
        B[0, 4] = b3
        B[1, 5] = g3
        B[2, 4] = g3
        B[2, 5] = b3
        
        return B * (1.0 / (2.0 * self.A))

    def _get_D_membrane(self) -> np.ndarray:
        """
        @brief Constitutive Matrix for Plane Stress (Isotropic).
        """
        E = self.mat.E
        nu = self.mat.nu
        factor = E / (1 - nu**2)
        
        D = np.array([
            [1,  nu, 0],
            [nu, 1,  0],
            [0,  0,  (1-nu)/2]
        ])
        return D * factor

    def _build_local_stiffness(self) -> np.ndarray:
        """
        @brief Computes Element Stiffness in Local Coordinates (18x18).
        """
        K_local = np.zeros((18, 18))
        
        # --- 1. MEMBRANE STIFFNESS (CST) ---
        # DOFs: u, v at each node (Tx, Ty)
        # Indices in 18x18 matrix:
        # N1: 0 (Tx), 1 (Ty)
        # N2: 6 (Tx), 7 (Ty)
        # N3: 12 (Tx), 13 (Ty)
        # Map 0->0, 1->1, 2->6, 3->7, 4->12, 5->13
        mem_indices_map = [0, 1, 6, 7, 12, 13]
        
        B_mem = self._get_B_membrane()
        D_mem = self._get_D_membrane()
        t = self.prop.t
        
        # K_mem = t * A * B.T * D * B
        k_mem = (self.prop.t * self.A) * (B_mem.T @ D_mem @ B_mem)
        
        # Scatter into K_local
        for i, row_idx in enumerate(mem_indices_map):
            for j, col_idx in enumerate(mem_indices_map):
                K_local[row_idx, col_idx] += k_mem[i, j]
                
        # --- 2. BENDING STIFFNESS ---
        # TODO: Implement DKT or Mindlin Plate. 
        # For now, we leave zero stiffness for bending to test Membrane first.
        # WARNING: This will cause singularity if out-of-plane forces exist or 
        # if the structure is 3D.
        # To avoid pure singularity in 3D, adding very small stiffness to diagonals of unhandled DOFs
        # eps = 1e-6 * np.trace(K_local) / 18.0
        # if eps == 0: eps = 1.0
        
        # Unhandled DOFs in this simplified versions: 
        # Tz(2, 8, 14), Rx(3, 9, 15), Ry(4, 10, 16), Rz(5, 11, 17)
        # We will NOT add artificial stiffness blindly, trusting the user uses this for plane stress
        # or accepting the singularity for now until Bending is implemented.
        # *Self-Correction*: User asked for a solver. If it's singular, it fails.
        # Adding a small number to diagonals of drilling rotation (Rz around local Z) is standard.
        # Drilling DOFs: 5, 11, 17.
        
        k_diag_small = 1e-8 * self.mat.E * self.A # Crude approx
        
        # Drilling Stabilization
        K_local[5, 5] += k_diag_small
        K_local[11, 11] += k_diag_small
        K_local[17, 17] += k_diag_small
        
        return K_local

    def get_stiffness_matrix(self) -> np.ndarray:
        """
        @brief Computes K_global = T.T * K_local * T
        """
        K_loc = self._build_local_stiffness()
        return self.T.T @ K_loc @ self.T

    def calculate_stresses(self, u_global_vector: np.ndarray, dof_mapper) -> dict:
        """
        @brief Calculates Membrane Stresses (Von Mises).
        """
        # 1. Get Global Displacements
        idxs = []
        for n in [self.node1, self.node2, self.node3]:
            idxs.extend(dof_mapper.get_dof_indices(n.nid))
            
        u_elem_global = u_global_vector[idxs]
        
        # 2. Local Displacements
        u_loc = self.T @ u_elem_global
        
        # 3. Membrane Strains
        # Extract u,v for CST
        # u_loc indices: 0,1, 6,7, 12,13
        u_mem = np.array([
            u_loc[0], u_loc[1],
            u_loc[6], u_loc[7],
            u_loc[12], u_loc[13]
        ])
        
        B_mem = self._get_B_membrane()
        strains = B_mem @ u_mem # [eps_x, eps_y, gam_xy]
        
        # 4. Membrane Stresses
        D_mem = self._get_D_membrane()
        stresses = D_mem @ strains # [sig_x, sig_y, tau_xy]
        
        sig_x, sig_y, tau_xy = stresses
        
        # Von Mises (Plane Stress)
        # sig_vm = sqrt(sig_x^2 - sig_x*sig_y + sig_y^2 + 3*tau^2)
        vm_stress = np.sqrt(sig_x**2 - sig_x*sig_y + sig_y**2 + 3*tau_xy**2)
        
        return {
            'sigma_x': sig_x,
            'sigma_y': sig_y,
            'tau_xy': tau_xy,
            'von_mises': vm_stress
        }
