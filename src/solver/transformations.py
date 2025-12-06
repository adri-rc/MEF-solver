import numpy as np

class CoordinateTransformation:
    """
    @class CoordinateTransformation
    @brief Calculation of local coordinate systems and transformations to global ones.
    """

    @staticmethod
    def compute_local_axes(node1_coords: np.ndarray, node2_coords: np.ndarray, orientation_vector: np.ndarray = None):
        """
        @brief Generate the 3 unit vectors (vx, vy, vz) of the local system.
        
        @param orientation_vector: Nastran vector 'v' for CBAR. 
                                   If it is None (CROD case), one is generated automatically.
        @return (vx, vy, vz) Unitary vectors numpy arrays (3,)
        """
        # 1. Local X-axis (Longitudinal)
        diff = node2_coords - node1_coords
        L = np.linalg.norm(diff)
        if L == 0:
            raise ValueError(f"Element with zero length between node1 ({node1_coords})  and node2 ({node2_coords})")
        vx = diff / L

        #2. Orientation Vector (Auxiliary)
        v_aux = orientation_vector # if provided
        
        if v_aux is None:
            # logic for CROD (Standard convention), try to use the Global Axis Y [0, 1, 0]
            v_aux = np.array([0.0, 1.0, 0.0])
            
            # Check if the bar is parallel to Y, the cross product would be 0. 
            # and Global X axis is used in that case. 
            if np.abs(np.dot(vx, v_aux)) > 0.99: # ~parallel
                v_aux = np.array([1.0, 0.0, 0.0])
        else:
            # For CBAR: Use the given vector, normalized
            v_aux = v_aux / np.linalg.norm(v_aux)

        # 3. Local Z axis (Normal to the X-Aux plane)
        # v_z = x CROSS v_aux
        vz = np.cross(vx, v_aux)
        norm_z = np.linalg.norm(vz)
        
        if norm_z < 1e-6:
             raise ValueError("The orientation vector is parallel to the axis of the element.")
        
        vz = vz / norm_z

        # 4. Local Y-axis (Completes the trihedron)
        # v_y = z CROSS x
        vy = np.cross(vz, vx)
        # vy ya es unitario
        
        return vx, vy, vz

    @staticmethod
    def get_rotation_matrix_3x3(vx, vy, vz):
        """
        @brief Construct the rotation matrix R (3x3) where v_local = R * v_global
        @details The rows are the unit vectors expressed in globals."""
        return np.array([vx, vy, vz])

    @staticmethod
    def build_full_transformation_matrix(R_3x3: np.ndarray, nodes_count: int = 2, dofs_per_node: int = 6):
        """
        Expand the 3x3 rotation matrix to the complete matrix of the element (T).
        For a beam element with 2 nodes and 6 DOFs/node, it generates a 12x12 diagonal by blocks.
        
        T = [ R  0  0  0 ]
            [ 0  R  0  0 ]
            [ 0  0  R  0 ]
            [ 0  0  0  R ]
        """
        # size: (nodos * dofs) x (nodos * dofs)
        size = nodes_count * dofs_per_node
        T = np.zeros((size, size))
        
        # Filling diagonal blocks
        # assume that the DOFs are ordered as: [Tx, Ty, Tz, Rx, Ry, Rz]
        # Therefore, R applies to [Tx, Ty, Tz] and also to [Rx, Ry, Rz].
        
        # num_blocks: Number of 3x3 blocks fitting
        # Para 2 nodos y 6 dofs: 12 dofs total -> 4 bloques de 3x3 (Trans1, Rot1, Trans2, Rot2)
        num_blocks = size // 3
        
        for i in range(num_blocks):
            idx = i * 3
            T[idx : idx+3, idx : idx+3] = R_3x3
            
        return T