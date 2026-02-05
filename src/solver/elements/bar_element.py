
"""
@file bar_element.py
@brief Finite Element Solver for the CBAR (Beam) element.
@details
Implements the stiffness matrix calculation for a 2-node Euler-Bernoulli Beam.
Includes Axial, Torsional, and Bending (2 planes) stiffness.
"""

import numpy as np
from src.preprocessing.data_model import Node, BarProperty, IsotropicMaterial
from ..quadrature import GaussQuadrature
from ..transformations import CoordinateTransformation

class BarElementSolver:
    """
    @class BarElementSolver
    @brief Solver logic for a 2-node Bar element (Euler-Bernoulli Beam).
    @details
    Computes the Global Stiffness Matrix (12x12) by integrating coupled physics:
    - Axial Stiffness (EA) (Linear Shape Functions)
    - Torsional Stiffness (GJ) (Linear Shape Functions)
    - Bending XY & XZ (Cubic Hermite Shape Functions)
        - XY Plane Bending Stiffness (EIzz) -> DOFs: v, theta_z
        - XZ Plane Bending Stiffness (EIyy) -> DOFs: w, theta_y

    Supported DOFs per node: 6 (Tx, Ty, Tz, Rx, Ry, Rz).
    
    Active DOFs in Local System: 
      - Tx (Axial)
      - Rx (Torsion)
      - Ty, Rz (Bending Plane XY)
      - Tz, Ry (Bending Plane XZ)
    """

    # --- CONSTANTS FOR DOF MAPPING (Local indices mapped to 12x12 matrix) ---
    # Global Local Vector indices: [Tx1, Ty1, Tz1, Rx1, Ry1, Rz1, Tx2, Ty2, Tz2, Rx2, Ry2, Rz2]
    # Indices:                      0    1    2    3    4    5    6    7    8    9    10   11
    
    
    IDX_AXIAL = [0, 6] # Axial depends on Tx at Node 1 (idx 0) and Node 2 (idx 6) -> Matrix 2x2
    IDX_TORSION = [3, 9] # Torsion depends on Rx at Node 1 (idx 3) and Node 2 (idx 9) -> Matrix 2x2
    IDX_BEND_XY = [1, 5, 7, 11] # Bending XY depends on Ty and Rz at both nodes -> Matrix 4x4, Order: [Ty1, Rz1, Ty2, Rz2]
    IDX_BEND_XZ = [2, 4, 8, 10]     # Bending XZ depends on Tz and Ry at both nodes -> Matrix 4x4, Order: [Tz1, Ry1, Tz2, Ry2]
    
    def __init__(self, node1: Node, node2: Node, prop: BarProperty, mat: IsotropicMaterial, orientation: list, n_gaussp: int = 2):
        self.node1 = node1
        self.node2 = node2
        self.prop = prop
        self.mat = mat
        self.n_gaussp = max(n_gaussp, 2) # Requireed a minimum of 2 gauss points to integrate exact cubic polynomials

        self.n1_coords = node1.get_coords()
        self.n2_coords = node2.get_coords()
        self.orientation_vector = np.array(orientation)
        
        diff = self.n2_coords - self.n1_coords
        self.L = np.linalg.norm(diff)
        
        if self.L < 1e-12:
            raise ValueError(f"[ERROR] Zero length element with id nodes {node1.nid}-{node2.nid}, in CBAR {prop.pid}")

        self._compute_transformation()
        
        # Pre-compute internal geometric props for speed
        self.J = self.L / 2.0 # Jacobian


    def _compute_transformation(self):
        """@brief Computes T matrix (12x12) based on orientation vector."""
        vx, vy, vz = CoordinateTransformation.compute_local_axes(
            self.n1_coords, self.n2_coords, orientation_vector=self.orientation_vector
        )
        R = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)
        self.T = CoordinateTransformation.build_full_transformation_matrix(
            R, nodes_count=2, dofs_per_node=6
        )
    
    # --- SHAPE FUNCTION DERIVATIVES ---

    def _get_B_axial_torsion(self) -> np.ndarray:
        """
        @brief Calculates Strain-Displacement matrix B at natural coordinate xi for linear problems (Axial/Torsion).
        @details 
        B = [dN1/dx, dN2/dx] = [dN1/dxi, dN2/dxi] * (1/J)
        @param J: Jacobian (L/2).
        @return Vector B = [-1/L, 1/L] in physical coords.
        """
        # dN/dxi = [-0.5, 0.5] -> dN/dx = dN/dxi * (1/J)
        dNe_dxi = np.array([-0.5, 0.5])
        return dNe_dxi * (1.0 / self.J)

    def _get_B_bending(self, xi: float) -> np.ndarray:
        """
        @brief Curvature-displacement matrix for bending (Hermite Cubic).
        @details
        Returns B vector of size 4 corresponding to [v1, th1, v2, th2].
        B contains d^2H/dx^2.
        @param xi: Natural coordinate [-1, 1].
        @param J: Jacobian (L/2).
        """
        # 1. Second derivatives w.r.t natural coord xi (d^2/dxi^2)
        # H1'' = 1.5*xi
        # H2'' = -0.5 + 1.5*xi
        # H3'' = -1.5*xi
        # H4'' = 0.5 + 1.5*xi
        
        ddH = np.array([
            1.5 * xi,           # v1
            -0.5 + 1.5 * xi,    # theta1 (associated with shape function * J)
            -1.5 * xi,          # v2
            0.5 + 1.5 * xi      # theta2 (associated with shape function * J)
        ])
        
        
        # 2. Chain Rule for Curvature: d^2/dx^2 = (1/J^2) * d^2/dxi^2
        inv_J2 = 1.0 / (self.J**2) # factor = (1/J)^2 for second derivative chain rule

        
        # 3. Scaling for Rotation DOFs:
        # Physical rotation theta = dv/dx.
        # The Shape functions for theta nodes are scaled by J implies N_theta = J * H_hat.
        # Thus, their derivatives must be multiplied by J.
        # Total scale for Rot terms: (1/J^2) * J = 1/J
        
        B = np.array([
            ddH[0] * inv_J2,        # v1 -> scaling 1/J^2
            ddH[1] * inv_J2 * self.J,    # theta1 -> scaling 1/J
            ddH[2] * inv_J2,        # v2 -> scaling 1/J^2
            ddH[3] * inv_J2 * self.J     # theta2 -> scaling 1/J
        ])
        return B
    
    def _build_local_stiffness(self) -> np.ndarray:
        """ 
        @brief Computes the local element stiffness matrix in Local Coordinates (12x12),
        needed for both Global K and Stress results recovery.
        @return Local element stiffness matrix: Numpy array (12, 12)
        """
        
        points, weights = GaussQuadrature.get_rule(self.n_gaussp)
        
        # Initialize sub-matrices
        k_axial = np.zeros((2, 2))
        k_torsion = np.zeros((2, 2))
        k_bend_xy = np.zeros((4, 4))
        k_bend_xz = np.zeros((4, 4))
        
        # Numerical Integration Loop
        for xi, w in zip(points, weights):
            detJ_w = self.J * w 
            
            # 1. Axial & Torsion
            # --- 1. AXIAL & TORSION ---
            B_lin = self._get_B_axial_torsion()
            # B_lin: due to in axial and torsion B is linear and the same value

            k_axial += np.outer(B_lin, B_lin) * self.mat.E * self.prop.A * detJ_w
            k_torsion += np.outer(B_lin, B_lin) * self.mat.G * self.prop.J * detJ_w
            
            # --- 2. BENDING ---
            B_bend = self._get_B_bending(xi)
            
            # Plane XY (I1)
            k_bend_xy += np.outer(B_bend, B_bend) * self.mat.E * self.prop.I1 * detJ_w
            
            # Plane XZ (I2) - Sign flip for RotY
            # Note: For XZ plane, positive rotation Ry causes negative displacement in Z 
            # based on right-hand rule, effectively changing signs in B-matrix for rot terms.
            B_bend_xz = B_bend.copy()
            B_bend_xz[1] *= -1.0 
            B_bend_xz[3] *= -1.0 
            k_bend_xz += np.outer(B_bend_xz, B_bend_xz) * self.mat.E * self.prop.I2 * detJ_w
            
        # --- ASSEMBLY INTO 12x12 LOCAL MATRIX ---
        K_local = np.zeros((12, 12))
        
        # Use the constants to map sub-matrices to the full matrix
        self._add_to_block(K_local, k_axial, self.IDX_AXIAL)
        self._add_to_block(K_local, k_torsion, self.IDX_TORSION)
        self._add_to_block(K_local, k_bend_xy, self.IDX_BEND_XY)
        self._add_to_block(K_local, k_bend_xz, self.IDX_BEND_XZ)
        
        return K_local
   
    def get_stiffness_matrix(self) -> np.ndarray:
        """
        @brief Transforms the local element stifness matrix into global coordinate system.
        @details Computes the element stiffness matrix in Global Coordinates (12x12) as:
        K_global = T.T * K_local * T
        @return Global Element Stiffness Matrix in global coordinates: Numpy array (12, 12)
        """
        K_loc = self._build_local_stiffness()
        
        return self.T.T @ K_loc @ self.T

    def _add_to_block(self, K_full: np.ndarray, K_sub: np.ndarray, indices: list[int]):
        """
        @brief Helper to scatter small matrix values into the 12x12 matrix.
        @param K_full The full 12x12 Local Stiffness Matrix being built.
        @param K_sub The smaller sub-matrix (2x2 for axial/torsion, 4x4 for bending) containing specific physics.
        @param indices: List of indices in the 12x12 matrix corresponding to the submatrix rows/cols.
                        Example: [0, 6] for axial maps K_sub[0,0] to K_full[0,0] and K_sub[0,1] to K_full[0,6].
        """
        for i, row_idx in enumerate(indices):
            for j, col_idx in enumerate(indices):
                K_full[row_idx, col_idx] += K_sub[i, j]
    

    def calculate_stresses(self, u_global_vector: np.ndarray, dof_mapper) -> dict:
        """
        @brief Recovers internal forces and stresses at specific points.
        @return Dict with keys 'forces' (dict) and 'stress_points' (dict).
        """
        # 1. Get Global Displacements
        idxs_n1 = dof_mapper.get_dof_indices(self.node1.nid)
        idxs_n2 = dof_mapper.get_dof_indices(self.node2.nid)
        u_elem_global = np.concatenate([u_global_vector[idxs_n1], u_global_vector[idxs_n2]])
        
        # 2. Transform global displacements to local element system
        u_loc = self.T @ u_elem_global
        
        
        # 2. Calculate Internal Forces via Stiffness: F_loc = K_loc * u_loc
        # This captures all coupled effects (if any) and correct signs.
        K_local = self._build_local_stiffness()
        f_local = K_local @ u_loc

        # Extract forces at Node 1 (End A) and Node 2 (End B)
        # Naming: Fx (Axial), Fy (Shear Y), Fz (Shear Z), Mx (Torque), My (Bend XZ), Mz (Bend XY)
        # Indices Node 1: 0,1,2,3,4,5
        forces_n1 = {
            'Fx': f_local[0], 'Fy': f_local[1], 'Fz': f_local[2],
            'Mx': f_local[3], 'My': f_local[4], 'Mz': f_local[5]
        }
        # Indices Node 2: 6,7,8,9,10,11 (Note: These are reaction forces on the element)
        forces_n2 = {
            'Fx': f_local[6], 'Fy': f_local[7], 'Fz': f_local[8],
            'Mx': f_local[9], 'My': f_local[10], 'Mz': f_local[11]
        }
        
        # 3. KINEMATIC RECOVERY (Better for Constant Strain/Twist)
        # --------------------------------------------------------
        
        # --- AXIAL ---
        # u_loc[0] = u_x1, u_loc[6] = u_x2
        axial_strain = (u_loc[6] - u_loc[0]) / self.L
        sigma_ax = self.mat.E * axial_strain
        axial_force_val = sigma_ax * self.prop.A
        # Note: The kinematically calculated axial force (sigma * A) should match
        # f_loc[6] (reaction at node 2). We use kinematics for stress reporting.
        
        # --- TORSION ---
        # u_loc[3] = theta_x1, u_loc[9] = theta_x2
        twist_rate = (u_loc[9] - u_loc[3]) / self.L
        torque_val = self.mat.G * self.prop.J * twist_rate
       
       
       
       
        # 4. STRESS RECOVERY AT POINTS (Navier's Formula)
        # --------------------------------------------------------
            # Sigma = Fx/A - (Mz * y)/I1 + (My * z)/I2
            # Tau   = Mx * r / J  (Simplified torsion)
        # --------------------------------------------------------

        stress_results = {}
        # If no points defined, use default point (0,0)
        points_to_calc = self.prop.stress_points if self.prop.stress_points else [(0.0, 0.0)]

        # Bending moment calculation at node 1
        # We retrieve the local nodal forces previously calculated (f_loc = K_loc * u_loc).
        # These forces represent the "External Reactions" required at the nodes to maintain equilibrium.
        
        # Sign Logic (Case Mz - XY Plane):
        # 1. FEM Convention: A positive Mz at Node 1 is a Counter-Clockwise rotation.
        # 2. Physical Effect: A Counter-Clockwise rotation at the left end lifts the beam,
        #    creating a positive curvature ("Smile" / Sagging shape).
        # 3. Navier Formula: sigma = - (Mz * y) / I
        #    - If we input a positive Mz (+) and evaluate at the top fiber where y is positive (+):
        #    - Result: sigma = - (+) * (+) = (-) -> COMPRESSION.
        #
        # Conclusion: The raw sign from the FEM matrix is compatible with the standard formula.
        # The sign MUST NOT be inverted.
        
        M_y = f_local[4] # Moment Y at Node 1 (reaction)
        M_z = f_local[5] # Moment Z at Node 1 (reaction)
        
        for i, (y, z) in enumerate(points_to_calc):
            # Normal Stress = (Axial + Bending)
            # Normal Stress = Axial + Bending
            # Navier Formula: Sigma = -Mz*y/I1 + My*z/I2

            sigma_bend_z = -(M_z * y) / self.prop.I1 if self.prop.I1 else 0
            sigma_bend_y = (M_y * z) / self.prop.I2 if self.prop.I2 else 0
            
            sigma_tot = sigma_ax + sigma_bend_z + sigma_bend_y
            
            # Shear Stress (Torsion)
            r = np.sqrt(y**2 + z**2)
            tau_tor = (torque_val * r) / self.prop.J if self.prop.J else 0
            
            stress_results[f"P{i}"] = {
                'sigma_x': sigma_tot,
                'tau_xy': tau_tor, # simplified
                'y': y, 'z': z
            }

        return {
            # --- CAMPOS UNIFICADOS (Compatibilidad con CROD) ---
            'axial_force': axial_force_val,
            'axial_stress': sigma_ax,
            'torque': torque_val,
            'torsional_stress': 0.0, # (Opcional: max shear stress)
            
            # --- CAMPOS ESPECÍFICOS DE VIGA ---
            'axial_strain': axial_strain,
            'twist_rate': twist_rate,
            'forces_n1': forces_n1,
            'forces_n2': forces_n2,
            'stress_points': stress_results
        }
