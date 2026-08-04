
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
    
    def __init__(
        self, node1: Node, node2: Node, prop: BarProperty,
        mat: IsotropicMaterial, orientation: list, n_gaussp: int = 2,
        wa=None, wb=None, pa: int = 0, pb: int = 0,
    ):
        self.node1 = node1
        self.node2 = node2
        self.prop = prop
        self.mat = mat
        self.n_gaussp = max(n_gaussp, 2) # Requireed a minimum of 2 gauss points to integrate exact cubic polynomials

        self.n1_coords = node1.get_coords()
        self.n2_coords = node2.get_coords()
        self.orientation_vector = np.asarray(orientation, dtype=float)
        if self.orientation_vector.shape != (3,) or not np.all(np.isfinite(self.orientation_vector)):
            raise ValueError("CBAR orientation must contain three finite components.")
        self.wa = np.zeros(3) if wa is None else np.asarray(wa, dtype=float)
        self.wb = np.zeros(3) if wb is None else np.asarray(wb, dtype=float)
        if self.wa.shape != (3,) or self.wb.shape != (3,):
            raise ValueError("CBAR offsets WA/WB must contain three components.")
        self.pa = int(pa or 0)
        self.pb = int(pb or 0)
        
        diff = self.n2_coords - self.n1_coords
        self.L = np.linalg.norm(diff)
        
        if self.L < 1e-12:
            raise ValueError(f"[ERROR] Zero length element with id nodes {node1.nid}-{node2.nid}, in CBAR {prop.pid}")

        self._compute_transformation()
        
        # Pre-compute internal geometric props for speed
        self.J = self.L / 2.0 # Jacobian

    @staticmethod
    def _skew(vector: np.ndarray) -> np.ndarray:
        x, y, z = vector
        return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])

    def _offset_transform_global(self) -> np.ndarray:
        """Map GRID DOFs to the two offset CBAR end DOFs."""
        transform = np.eye(12)
        transform[0:3, 3:6] = -self._skew(self.wa)
        transform[6:9, 9:12] = -self._skew(self.wb)
        return transform

    @staticmethod
    def _timoshenko_bending(EI: float, kGA: float, length: float, flip_rotation=False):
        phi = 12.0 * EI / (kGA * length**2)
        scale = EI / (length**3 * (1.0 + phi))
        matrix = scale * np.array([
            [12.0, 6.0*length, -12.0, 6.0*length],
            [6.0*length, (4.0+phi)*length**2, -6.0*length, (2.0-phi)*length**2],
            [-12.0, -6.0*length, 12.0, -6.0*length],
            [6.0*length, (2.0-phi)*length**2, -6.0*length, (4.0+phi)*length**2],
        ])
        if flip_rotation:
            signs = np.diag([1.0, -1.0, 1.0, -1.0])
            matrix = signs @ matrix @ signs
        return matrix

    def _apply_releases(self, stiffness: np.ndarray) -> np.ndarray:
        released = []
        for flag, offset in ((self.pa, 0), (self.pb, 6)):
            for component in str(flag):
                if component in "123456":
                    released.append(offset + int(component) - 1)
        released = sorted(set(released))
        if not released:
            return stiffness
        retained = [index for index in range(12) if index not in released]
        k_rr = stiffness[np.ix_(retained, retained)]
        k_rq = stiffness[np.ix_(retained, released)]
        k_qq = stiffness[np.ix_(released, released)]
        condensed = k_rr - k_rq @ np.linalg.pinv(k_qq, rcond=1e-12) @ k_rq.T
        result = np.zeros_like(stiffness)
        result[np.ix_(retained, retained)] = condensed
        return 0.5 * (result + result.T)


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
        k_bend_coupling = np.zeros((8, 8))
        
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

            i12 = getattr(self.prop, 'I12', 0.0)
            if i12 != 0.0:
                B_pair = np.zeros((2, 8))
                B_pair[0, :4] = B_bend
                B_pair[1, 4:] = B_bend_xz
                inertia_coupling = self.mat.E * np.array([
                    [0.0, i12], [i12, 0.0]
                ])
                k_bend_coupling += B_pair.T @ inertia_coupling @ B_pair * detJ_w

        if getattr(self.prop, 'K1', 0.0) > 0.0:
            k_bend_xy = self._timoshenko_bending(
                self.mat.E * self.prop.I1,
                self.prop.K1 * self.prop.A * self.mat.G,
                self.L,
            )
        if getattr(self.prop, 'K2', 0.0) > 0.0:
            k_bend_xz = self._timoshenko_bending(
                self.mat.E * self.prop.I2,
                self.prop.K2 * self.prop.A * self.mat.G,
                self.L,
                flip_rotation=True,
            )
            
        # --- ASSEMBLY INTO 12x12 LOCAL MATRIX ---
        K_local = np.zeros((12, 12))
        
        # Use the constants to map sub-matrices to the full matrix
        self._add_to_block(K_local, k_axial, self.IDX_AXIAL)
        self._add_to_block(K_local, k_torsion, self.IDX_TORSION)
        self._add_to_block(K_local, k_bend_xy, self.IDX_BEND_XY)
        self._add_to_block(K_local, k_bend_xz, self.IDX_BEND_XZ)
        self._add_to_block(
            K_local, k_bend_coupling, self.IDX_BEND_XY + self.IDX_BEND_XZ
        )
        
        return self._apply_releases(K_local)
   
    def get_stiffness_matrix(self) -> np.ndarray:
        """
        @brief Transforms the local element stifness matrix into global coordinate system.
        @details Computes the element stiffness matrix in Global Coordinates (12x12) as:
        K_global = T.T * K_local * T
        @return Global Element Stiffness Matrix in global coordinates: Numpy array (12, 12)
        """
        K_loc = self._build_local_stiffness()
        
        k_offset_ends = self.T.T @ K_loc @ self.T
        offset_transform = self._offset_transform_global()
        return offset_transform.T @ k_offset_ends @ offset_transform

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
        @return Dict with keys 'forces', 'stress_points', 'max_total_stress', 'min_total_stress'.
        """
        # 1. Get Global Displacements
        idxs_n1 = dof_mapper.get_dof_indices(self.node1.nid)
        idxs_n2 = dof_mapper.get_dof_indices(self.node2.nid)
        u_elem_global = np.concatenate([u_global_vector[idxs_n1], u_global_vector[idxs_n2]])
        
        # 2. Transform global displacements to local element system
        u_offset_global = self._offset_transform_global() @ u_elem_global
        u_loc = self.T @ u_offset_global
        
        
        # 3. Calculate Internal Forces via Stiffness: F_loc = K_loc * u_loc
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
        
        # 4. KINEMATIC RECOVERY (Better for Constant Strain/Twist)
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
       
       
       
       
        # 5. STRESS RECOVERY AT POINTS (Navier's Formula)
        # --------------------------------------------------------
            # Sigma = Fx/A - (Mz * y)/I1 + (My * z)/I2
            # Tau   = Mx * r / J  (Simplified torsion)
            # Mz positive compress +y fiber.
            # My positive compress +z fiber.
        # --------------------------------------------------------
        # Convention Conversion:
        # f_local contains NODAL forces.
        # Internal Section Forces at End A = - (Nodal Forces at End A)
        # Internal Section Forces at End B = + (Nodal Forces at End B)

        # We calculate stresses usually at End A (or we could interpolate).
        # Here we use End A Internal Forces:
        M_y_int = -f_local[4]  # Internal Bending Moment Plane 2 (about local y)
        M_z_int = -f_local[5]  # Internal Bending Moment Plane 1 (about local z)
        
        # Note regarding Axial contribution:
        # We use sigma_ax calculated kinematically above for consistency,
        # but technically P/A = (-f_local[0]) / A.
        
        stress_results = {}
        stress_results_end_b = {}
        # If no points defined, use default point (0,0)
        points_to_calc = self.prop.stress_points if self.prop.stress_points else [(0.0, 0.0)]

        # List to collect all total stresses for Max/Min calculation
        computed_total_stresses = []
        i12 = getattr(self.prop, 'I12', 0.0)
        inertia = np.array([[self.prop.I1, i12], [i12, self.prop.I2]])
        end_moments = (
            (-f_local[5], -f_local[4]),
            ( f_local[11], f_local[10]),
        )
        stress_sets = (stress_results, stress_results_end_b)
        
        for result_set, (moment_z, moment_y) in zip(stress_sets, end_moments):
            # Coefficients of sigma_b = a*y + b*z for a general I12.
            coefficients = np.linalg.solve(inertia, np.array([-moment_z, moment_y]))
            for i, (y, z) in enumerate(points_to_calc):
                sigma_bending_total = coefficients[0] * y + coefficients[1] * z
                sigma_tot = sigma_ax + sigma_bending_total
                computed_total_stresses.append(sigma_tot)
                radius = np.hypot(y, z)
                tau_tor = (torque_val * radius) / self.prop.J if self.prop.J else 0.0
                result_set[f"P{i}"] = {
                    'sigma_total': sigma_tot,
                    'sigma_bending': sigma_bending_total,
                    'sigma_axial': sigma_ax,
                    'tau_xy': tau_tor,
                    'y': y, 'z': z
                }
        
        # Calculate Max and Min Total Stress (Equivalent to SA-MAX / SA-MIN from .f06  CBAR results)
        end_a_totals = [point['sigma_total'] for point in stress_results.values()]
        max_total_stress = max(end_a_totals) if end_a_totals else 0.0
        min_total_stress = min(end_a_totals) if end_a_totals else 0.0
        max_element_stress = max(computed_total_stresses) if computed_total_stresses else 0.0
        min_element_stress = min(computed_total_stresses) if computed_total_stresses else 0.0

        return {
            # --- CAMPOS UNIFICADOS (CROD compatibility) ---
            'axial_force': axial_force_val,
            'axial_stress': sigma_ax,
            'torque': torque_val,
            'torsional_stress': max(
                (abs(point['tau_xy']) for points in stress_sets for point in points.values()),
                default=0.0,
            ),
            
            # --- CAMPOS ESPECÍFICOS DE VIGA ---
            'axial_strain': axial_strain,
            'twist_rate': twist_rate,
            'forces_n1': forces_n1,
            'forces_n2': forces_n2,
            'stress_points': stress_results,
            'stress_points_end_a': stress_results,
            'stress_points_end_b': stress_results_end_b,

            # --- EXTREMES (Comparison with SA-MAX / SA-MIN from .f06  CBAR results) ---
            'max_total_stress': max_total_stress,
            'min_total_stress': min_total_stress,
            'max_element_stress': max_element_stress,
            'min_element_stress': min_element_stress,
        }
