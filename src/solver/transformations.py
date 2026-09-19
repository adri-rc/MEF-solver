"""Coordinate-system utilities for element matrices and vectors."""

import numpy as np


class CoordinateTransformation:
    """Build orthonormal element frames and block transformations."""
    
    @staticmethod
    def compute_local_axes_triangle(
        node1_coords: np.ndarray,
        node2_coords: np.ndarray,
        node3_coords: np.ndarray,
    ):
        """Return a right-handed local basis for a triangular element.

        The origin is node 1, local x follows node 1 to node 2, and local z is
        normal to the element plane.
        """
        # Local x follows node 1 to node 2.
        v12 = node2_coords - node1_coords
        L = np.linalg.norm(v12)
        if L < 1e-12:
             raise ValueError("Zero length edge N1-N2")
        x_axis = v12 / L

        # Node 1 to node 3 defines the element plane.
        v13 = node3_coords - node1_coords
        
        # Local z is normal to the element plane.
        z_axis = np.cross(x_axis, v13)
        norm_z = np.linalg.norm(z_axis)
        if norm_z < 1e-12:
            raise ValueError("Nodes are collinear, cannot define a triangle plane.")
        z_axis = z_axis / norm_z

        # Local y completes the in-plane right-handed triad.
        y_axis = np.cross(z_axis, x_axis) 

        # Keep node 3 on the positive side of local y.
        if np.dot(v13, y_axis) < 0:
             y_axis = -y_axis
             # Flip z as well to preserve a right-handed basis.
             z_axis = -z_axis

        
        return x_axis, y_axis, z_axis
    

    @staticmethod
    def compute_local_axes(
        node1_coords: np.ndarray,
        node2_coords: np.ndarray,
        orientation_vector: np.ndarray = None,
    ):
        """Return the local basis for a two-node element.

        ``orientation_vector`` is the Nastran orientation vector for CBAR. A
        stable auxiliary axis is generated automatically for CROD.
        """
        # Local x is the longitudinal direction.
        diff = node2_coords - node1_coords
        L = np.linalg.norm(diff)
        if L == 0:
            raise ValueError(f"Element with zero length between node1 ({node1_coords})  and node2 ({node2_coords})")
        vx = diff / L

        # Auxiliary orientation vector.
        v_aux = orientation_vector
        
        if v_aux is None:
            # Prefer basic y, except when the element is almost parallel to it.
            v_aux = np.array([0.0, 1.0, 0.0])
            
            if np.abs(np.dot(vx, v_aux)) > 0.99:
                v_aux = np.array([1.0, 0.0, 0.0])
        else:
            # For CBAR: Use the given vector, normalized
            v_aux = np.asarray(v_aux, dtype=float)
            norm_aux = np.linalg.norm(v_aux)
            if not np.isfinite(norm_aux) or norm_aux < 1e-12:
                raise ValueError("The orientation vector must be finite and non-zero.")
            v_aux = v_aux / norm_aux

        # Local z is normal to the x/auxiliary plane.
        # v_z = x CROSS v_aux
        vz = np.cross(vx, v_aux)
        norm_z = np.linalg.norm(vz)
        
        if norm_z < 1e-6:
             raise ValueError("The orientation vector is parallel to the axis of the element.")
        
        vz = vz / norm_z

        # Local y completes the right-handed triad.
        # v_y = z CROSS x
        vy = np.cross(vz, vx)
        # vy is already normalized because vx and vz are orthonormal.
        
        return vx, vy, vz

    @staticmethod
    def get_rotation_matrix_3x3(vx, vy, vz):
        """Return ``R`` such that ``v_local = R @ v_global``."""
        return np.array([vx, vy, vz])

    @staticmethod
    def build_full_transformation_matrix(R_3x3: np.ndarray, nodes_count: int = 2, dofs_per_node: int = 6):
        """Expand a 3x3 rotation into the element block-diagonal matrix."""
        # Matrix size: (nodes * DOFs) x (nodes * DOFs).
        size = nodes_count * dofs_per_node
        T = np.zeros((size, size))
        
        # Filling diagonal blocks
        # assume that the DOFs are ordered as: [Tx, Ty, Tz, Rx, Ry, Rz]
        # Therefore, R applies to [Tx, Ty, Tz] and also to [Rx, Ry, Rz].
        
        # num_blocks: Number of 3x3 blocks fitting
        # Two six-DOF nodes require four 3x3 blocks: translation and rotation
        # for each node.
        num_blocks = size // 3
        
        for i in range(num_blocks):
            idx = i * 3
            T[idx : idx+3, idx : idx+3] = R_3x3
            
        return T
