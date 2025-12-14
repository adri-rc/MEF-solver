import numpy as np
from ...preprocessing.data_model import Node, RodProperty, IsotropicMaterial
from ..quadrature import GaussQuadrature


class RodElementSolver:
    """
    @class RodElementSolver
    @brief Gets the element stiffness matrix of a Rod element (CROD) (2nodes; 1 integration point)
    @details
    Rod supports axial load (Ttension/compression) and torsional.
    Calculate the stiffness matrix by integrating B.T * D * B.
    """
    def __init__(self, node1: Node, node2: Node, prop: RodProperty, mat: IsotropicMaterial, n_gaussp: int = 1):
        """
        @param n_gaussp Number of Gauss integration points (default 1 for constant CROD).
        """
        self.n1 = node1
        self.n2 = node2
        self.prop = prop
        self.mat = mat
        self.n_gaussp = n_gaussp

        self.n1_coords = node1.get_coords()
        self.n2_coords = node2.get_coords()
        
        diff = self.n1_coords - self.n2_coords
        self.L = np.linalg.norm(diff)
        if self.L == 0:
            raise ValueError(f"Element with zero length between nodes {node1.nid} and {node2.nid}")
        
        self.n_vec = diff / self.L # Unitary Director Vector (local x-axis)
       
    def _get_B_matrix(self, xi: float, J: float) -> np.ndarray:
        """
        Calculate matrix Be (Strain-Displacement) at point xi.
        Be = [dN1/dx, dN2/dx] = [dN1/dxi, dN2/dxi] * (1/J)
        """
        # Derivatives of shape linear functions linear with respect to xi
        # N1 = 0.5*(1-xi) -> dN1 = -0.5
        # N2 = 0.5*(1+xi) -> dN2 =  0.5
        dNe_dxi = np.array([-0.5, 0.5])
        
        # Transformación al espacio físico (dx)????
        return dNe_dxi * (1.0 / J)
    
    def get_stiffness_matrix(self) -> np.ndarray:
        """
        @brief Calculate and assemble the 12x12 global stiffness matrix (6 DOFs per node).
        """
        
        points, weights = GaussQuadrature.get_rule(self.n_gaussp) # Default n_gaussp =1
        
        # 1. Initialization of local scalar rigidities. (2x2)
        # k_axial_local --> u_local_1, u_local_2 (displacements)
        # k_tor_local   --> theta_local_1, theta_local_2 (rotations)
        k_axial_local = np.zeros((2, 2))
        k_tor_local = np.zeros((2, 2))

        # 2. NUMERICAL INTEGRATION LOOP
        for xi, w in zip(points, weights):
            
            J = self.L / 2.0 # Jacobian for linear element with 2 nodes (L/2)
            Be = self._get_B_matrix(xi, J) # Be matrix at the integration point (Row vector of size 2)
            
            # --- INTEGRACIÓN AXIAL ---
            # K = Integral( Be.T * E * A * Be * detJ * dxi )
            # np.outer(B, B) genera la matriz 2x2 (B.T * B)
            integrand_axial = np.outer(Be, Be) * self.mat.E * self.prop.A * J
            k_axial_local += integrand_axial * w
            
            # --- INTEGRACIÓN TORSIONAL ---
            # K = Integral( B.T * G * J_polar * B * detJ * dxi )
            # Nota: Usamos la misma B porque la interpolación del ángulo de giro es lineal igual que u.
            integrand_tor = np.outer(Be, Be) * self.mat.G * self.prop.J * J
            k_tor_local += integrand_tor * w
        
        # 3. TRANSFORMACIÓN A GLOBAL Y ENSAMBLAJE (12x12)
        # Ahora tenemos las matrices 2x2 puras. Necesitamos proyectarlas al espacio 3D.
        
        # Extraemos el coeficiente escalar de rigidez (k11 = EA/L numérico)
        k_ax_val = k_axial_local[0, 0]
        k_tor_val = k_tor_local[0, 0]
        
        # Matriz de Proyección P = n * n.T (3x3)
        # Proyecta un vector 3D sobre el eje de la barra
        P = np.outer(self.n_vec, self.n_vec)
        
        # Submatrices Globales 3x3
        K_sub_ax = k_ax_val * P   # Rigidez traslacional
        K_sub_tor = k_tor_val * P # Rigidez rotacional
        
        # Inicializamos matriz global 12x12
        # Estructura de DOFs: [Tx1, Ty1, Tz1, Rx1, Ry1, Rz1, Tx2, Ty2, Tz2, Rx2, Ry2, Rz2]
        K_global = np.zeros((12, 12))
        
        # Bloques de índices para los nodos 1 y 2
        # T (Traslación): 0:3, R (Rotación): 3:6
        idx_T1, idx_R1 = slice(0, 3), slice(3, 6)
        idx_T2, idx_R2 = slice(6, 9), slice(9, 12)

        # --- ENSAMBLAJE AXIAL (Afecta a Traslaciones T) ---
        # Fuerza axial opuesta en nodos: [ K  -K ]
        #                                [ -K  K ]
        K_global[idx_T1, idx_T1] += K_sub_ax
        K_global[idx_T1, idx_T2] -= K_sub_ax
        K_global[idx_T2, idx_T1] -= K_sub_ax
        K_global[idx_T2, idx_T2] += K_sub_ax
        
        # --- ENSAMBLAJE TORSIONAL (Afecta a Rotaciones R) ---
        # Momento torsor opuesto en nodos
        K_global[idx_R1, idx_R1] += K_sub_tor
        K_global[idx_R1, idx_R2] -= K_sub_tor
        K_global[idx_R2, idx_R1] -= K_sub_tor
        K_global[idx_R2, idx_R2] += K_sub_tor
        
        return K_global

















