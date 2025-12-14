"""
@file rod_element_v2.py
@brief Finite Element Solver for the CROD element.
@details 
Implements the stiffness matrix calculation for a 2-node Rod element
considering Axial and Torsional degrees of freedom.
"""
import numpy as np
from src.preprocessing.data_model import Node, RodProperty, IsotropicMaterial
from ..quadrature import GaussQuadrature
from ..transformations import CoordinateTransformation 


class RodElementSolverV2:
    """
    @class RodElementSolverV2
    @brief GeSolver logic for a 2-node Rod element (Linear interpolation).
    @details
    Computes the Global Stiffness Matrix (12x12) by:
    1. Calculating Local Stiffness (Integration of B.T * D * B).
    2. Applying Coordinate Transformation (T.T * K_loc * T).
    
    Supported DOFs per node: 6 (Tx, Ty, Tz, Rx, Ry, Rz).
    Active DOFs in Local System: Tx (Axial), Rx (Torsion).
    """
    # Class constants for DOF mapping (Local indices mapped to 12x12 matrix)
    # Node 1: Tx=0, Ty=1, Tz=2, Rx=3...
    # Node 2: Tx=6, Ty=7, ...
    IDX_AXIAL = [0, 6] 
    IDX_TORSION = [3, 9]

    def __init__(self, node1: Node, node2: Node, prop: RodProperty, mat: IsotropicMaterial, n_gaussp: int = 1):
        
        self.node1 = node1
        self.node2 = node2
        self.prop = prop
        self.mat = mat
        self.n_gaussp = n_gaussp

        self.n1_coords = node1.get_coords()
        self.n2_coords = node2.get_coords()
        diff = self.n2_coords - self.n1_coords
        self.L = np.linalg.norm(diff)
        if self.L < 1e-12:
            raise ValueError(f"[ERROR] Zero length element with id nodes {node1.nid}-{node2.nid}")

        self._compute_transformation() # Pre-compute Transformation Matrix (Expensive operation done once)

    def _compute_transformation(self):
        """Helper to build the 12x12 Transformation Matrix T."""
        vx, vy, vz = CoordinateTransformation.compute_local_axes(
            self.n1_coords, self.n2_coords, orientation_vector=None
        )
        R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)
        self.T = CoordinateTransformation.build_full_transformation_matrix(
            R, nodes_count=2, dofs_per_node=6
        )

    def _get_B_matrix(self, xi: float, J: float) -> np.ndarray:
        """
        @brief Calculates Strain-Displacement matrix B at natural coordinate xi.
        @details
        B = [dN1/dx, dN2/dx] = [dN1/dxi, dN2/dxi] * (1/J)
        """
        # Derivatives of shape linear functions linear with respect to xi
        # For linear 2-node rod, shape functions:
        # N1 = 0.5*(1-xi) -> dN1 = -0.5
        # N2 = 0.5*(1+xi) -> dN2 =  0.5
        dNe_dxi = np.array([-0.5, 0.5]) #linear 2-node rod
        
        return dNe_dxi * (1.0 / J)
    
    def get_stiffness_matrix(self) -> np.ndarray:
        """
        @brief Computes the element stiffness matrix in Global Coordinates.
        @return Numpy array (12, 12)
        """
        points, weights = GaussQuadrature.get_rule(self.n_gaussp) # Default n_gaussp =1
        
        # Initialize Local Stiffness Sub-matrices (2x2)
        k_axial_local = np.zeros((2, 2))
        k_tor_local = np.zeros((2, 2))

        # Numerical Integration Loop
        # Jacobian is constant for a straight linear rod: L/2
        J = self.L / 2.0
        for xi, w in zip(points, weights):
                      
            B = self._get_B_matrix(xi, J) # Strain-Displacement Matrix (1x2)
            
            # --- AXIAL ---
            # K = Integral( B.T * E * A * B * detJ * dxi )
            # Outer product handles the matrix multiplication (2x1) * (1x2) = (2x2)
            integrand_axial = np.outer(B, B) * self.mat.E * self.prop.A * J
            k_axial_local += integrand_axial * w
            
            # --- TORSIONAL ---
            # K = Integral( B.T * G * J_polar * B * detJ * dxi )
            # same B because the interpolation of the rotation angle is linear, like displacements (u).
            integrand_tor = np.outer(B, B) * self.mat.G * self.prop.J * J
            k_tor_local += integrand_tor * w
        
        # Assemble into 12x12 Local Matrix
        K_local = self._assemble_local_matrix(k_axial_local, k_tor_local)
        
        # Transform to Global Coordinates: K_global = T.T * K_local * T
        return self.T.T @ K_local @ self.T
    
    def _assemble_local_matrix(self, k_axial, k_tor):
        """
        @brief Assembles the 2x2 physics matrices into the sparse 12x12 element matrix.
        """
        K = np.zeros((12, 12))

        # k_axial is 2x2. relations (ij_loc_reduced->ij_loc_full): (0,0)->(0,0), (0,1)->(0,6), (1.0)->(6,0), (1,1)->(6,6)
        # K reduced is the 2x2 matrix, meanwhile the full K_loc_full is withs all DOFs of the element.
        for i_local, i_global in enumerate(self.IDX_AXIAL):
            for j_local, j_global in enumerate(self.IDX_AXIAL):
                K[i_global, j_global] = k_axial[i_local, j_local]
        # k_tor is 2x2. relations (ij_loc_reduced->ij_loc_full): (0,0)->(3,3), (0,1)->(3,9), (1,0)->(9,3), (1,1)->(9,9)
        for i_local, i_global in enumerate(self.IDX_TORSION):
            for j_local, j_global in enumerate(self.IDX_TORSION):
                K[i_global, j_global] = k_tor[i_local, j_local]
                
        return K
    
    def calculate_stresses(self, u_global_vector: np.ndarray, dof_mapper) -> dict:
        """
        @brief Calculates element internal forces and stresses from global displacements.
        @details
        Performs the following steps:
          1. Extracts global displacements for the element nodes.
          2. Transforms global displacements to local coordinates (u_local = T * u_global).
          3. Computes strains, stresses, and forces based on constitutive equations.
        
        Calculated outputs:
          - Axial Force (N) and Axial Stress (sigma)
          - Torque (T) and Torsional Shear Stress (tau)

        @param u_global_vector The full solution vector of the system (size N_dofs).
        @param dof_mapper Object responsible for mapping Node IDs to global indices.
        @return Dictionary with keys: 'axial_force', 'axial_stress', 'torque', 'torsional_stress'.
        """
        # 1. Extract Global Displacements for this element
        idxs_n1 = dof_mapper.get_dof_indices(self.node1.nid)
        idxs_n2 = dof_mapper.get_dof_indices(self.node2.nid)
        
        # Concatenate to get the 12x1 vector [u_n1, u_n2]
        u_elem_global = np.concatenate([
            u_global_vector[idxs_n1], 
            u_global_vector[idxs_n2]
        ])
        # 2. Transform to Local Coordinates: 
        # u_local = T * u_global
        u_local = self.T @ u_elem_global
        
        # 3. Physics Calculation
        # --- AXIAL (Tx) ---
        # Indices: 0 (Node 1 Tx), 6 (Node 2 Tx) in local system
        u_ax_1, u_ax_2 = u_local[0], u_local[6]
        axial_strain = (u_ax_2 - u_ax_1) / self.L
        axial_stress = self.mat.E * axial_strain # Hooke's Law: sigma = E * epsilon
        axial_force = axial_stress * self.prop.A

        # --- TORSION (Rx) ---
        # Indices: 3 (Node 1 Rx), 9 (Node 2 Rx) in local system
        th_x1, th_x2 = u_local[3], u_local[9]
        twist_rate = (th_x2 - th_x1) / self.L
        torque = self.mat.G * self.prop.J * twist_rate
        
        # Torsional Shear Stress (Tau) = (T * C) / J
        c_val = getattr(self.prop, 'C', 0.0) # Coeficiente de fibra más alejada
        
        torsional_stress = 0.0
        if self.prop.J > 1e-12:
            torsional_stress = (torque * c_val) / self.prop.J

        return {
            'axial_force': axial_force,
            'axial_stress': axial_stress,
            'torque': torque,
            'torsional_stress': torsional_stress
        }



