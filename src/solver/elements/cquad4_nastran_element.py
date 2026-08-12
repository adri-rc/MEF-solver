"""MSC/Nastran-style QUAD4 stiffness formulation.

The in-plane and curvature blocks use a uniform-strain field plus the two
orthogonal bilinear modes used by QUAD4.  The geometric hourglass metric is
built from the two element mid-lines, which is essential for skewed and
distorted quadrilaterals.  Non-coplanar nodes receive an equilibrium-based
force and moment correction after the mean-plane stiffness is formed.
"""

import numpy as np

from ..quadrature import GaussQuadrature
from .cquad4_element import CQUAD4ElementBase


class CQUAD4NastranElementSolver(CQUAD4ElementBase):
    """Four-node shell reconstructing the Nastran QUAD4 stiffness.

    Flat elements use uniform membrane/curvature strains plus the two
    bilinear stabilization modes.  Warped elements additionally move the
    mean-plane forces and moments to the four physical nodes before the
    usual element-to-global coordinate transformation.
    """

    formulation = "nastran"

    def _build_nastran_plane_stiffness(
        self, constitutive: np.ndarray, young_modulus: float, section_scale: float
    ) -> np.ndarray:
        """Build one Nastran-style in-plane block on the element mean plane."""
        B_center, _, _, det_jacobian = self._get_B_matrices(0.0, 0.0)
        area = 4.0 * det_jacobian
        stiffness = (
            B_center.T @ constitutive @ B_center
        ) * (area * section_scale)

        _, dN_dxi, dN_deta = self._get_shape_functions(0.0, 0.0)
        jacobian, _ = self._get_jacobian(dN_dxi, dN_deta)
        center_gradients = np.linalg.inv(jacobian) @ np.vstack(
            (dN_dxi, dN_deta)
        )
        base_mode = np.array([1.0, -1.0, 1.0, -1.0])
        hourglass_vector = base_mode - center_gradients.T @ (
            self.local_coords_2d.T @ base_mode
        )
        hourglass_operator = np.zeros((8, 2))
        hourglass_operator[0::2, 0] = hourglass_vector
        hourglass_operator[1::2, 1] = hourglass_vector

        # Vectors joining opposite-side midpoints.  Unlike a metric based on
        # the inverse Jacobian at Gauss points, this expression retains the
        # skew coupling present in Nastran's distorted QUAD4 stiffness.
        midline_x = 0.5 * (
            self.local_coords_2d[1]
            + self.local_coords_2d[2]
            - self.local_coords_2d[0]
            - self.local_coords_2d[3]
        )
        midline_y = 0.5 * (
            self.local_coords_2d[2]
            + self.local_coords_2d[3]
            - self.local_coords_2d[0]
            - self.local_coords_2d[1]
        )
        metric = area * (
            np.outer(midline_x, midline_x) / np.dot(midline_x, midline_x) ** 2
            + np.outer(midline_y, midline_y) / np.dot(midline_y, midline_y) ** 2
        )
        hourglass_constitutive = young_modulus * metric / 12.0
        stiffness += (
            hourglass_operator
            @ hourglass_constitutive
            @ hourglass_operator.T
            * section_scale
        )
        return 0.5 * (stiffness + stiffness.T)

    def _build_local_stiffness(self) -> np.ndarray:
        """Assemble membrane, bending, shear and drilling mean-plane blocks."""
        stiffness = np.zeros((24, 24))
        thickness = self.prop.t
        bending_inertia = thickness**3 / 12.0
        shear_area = (5.0 / 6.0) * thickness

        membrane = self._build_nastran_plane_stiffness(
            self._get_D_membrane(), self.mat_membrane.E, thickness
        )
        bending_plane = self._build_nastran_plane_stiffness(
            self._get_D_bending(), self.mat_bending.E, bending_inertia
        )
        rotation_to_plane = np.zeros((8, 12))
        for node_index in range(4):
            rotation_to_plane[2 * node_index, 3 * node_index + 2] = 1.0
            rotation_to_plane[2 * node_index + 1, 3 * node_index + 1] = -1.0
        bending = rotation_to_plane.T @ bending_plane @ rotation_to_plane

        points, weights = GaussQuadrature.get_rule(2)
        shear_full = np.zeros((12, 12))
        for i, xi in enumerate(points):
            for j, eta in enumerate(points):
                _, _, B_shear, det_jacobian = self._get_B_matrices(xi, eta)
                differential_area = det_jacobian * weights[i] * weights[j]
                shear_full += (
                    B_shear.T @ self._get_D_shear() @ B_shear
                ) * (differential_area * shear_area)

        _, center_weights = GaussQuadrature.get_rule(1)
        _, _, B_shear_center, det_jacobian = self._get_B_matrices(0.0, 0.0)
        shear_reduced = (
            B_shear_center.T @ self._get_D_shear() @ B_shear_center
        ) * (det_jacobian * center_weights[0] ** 2 * shear_area)
        shear = shear_reduced + 1.0e-3 * (shear_full - shear_reduced)

        membrane_dofs = [0, 1, 6, 7, 12, 13, 18, 19]
        plate_dofs = [2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22]
        self._scatter_block(stiffness, membrane, membrane_dofs)
        if thickness > 0.0:
            self._scatter_block(stiffness, bending + shear, plate_dofs)
        self._add_drilling_stiffness(stiffness)
        return 0.5 * (stiffness + stiffness.T)

    def _build_warping_force_transformation(self) -> np.ndarray:
        """Map mean-plane resultants to the actual warped corner nodes.

        Nastran treats a non-coplanar quadrilateral as a flat element on the
        plane through the four edge midpoints.  Moving its resultants to the
        real nodes requires two equilibrium corrections:

        * in-plane forces are resolved along both incident edges and produce
          self-equilibrated transverse ``kick-off`` forces;
        * bending moments are made tangent to each nodal surface and their
          added normal components are balanced by an in-plane force system.

        The returned matrix acts on nodal *forces*.  Consequently the warped
        stiffness is obtained by the congruence ``W @ K_mean @ W.T``.
        """
        warping = np.eye(24)
        heights = (self.coords - self.origin) @ self.local_z_axis
        characteristic_length = max(
            np.linalg.norm(
                self.local_coords_2d[(node + 1) % 4]
                - self.local_coords_2d[node]
            )
            for node in range(4)
        )
        if np.max(np.abs(heights)) <= 1.0e-12 * characteristic_length:
            return warping

        # At each corner resolve (Fx, Fy) into the directions of its outgoing
        # and incoming edges.  Each edge component is then transferred from
        # the mean-plane endpoints to their true elevations while preserving
        # the edge's force and moment equilibrium.
        edge_resolution = []
        for node in range(4):
            outgoing = (
                self.local_coords_2d[(node + 1) % 4]
                - self.local_coords_2d[node]
            )
            incoming = (
                self.local_coords_2d[(node - 1) % 4]
                - self.local_coords_2d[node]
            )
            outgoing /= np.linalg.norm(outgoing)
            incoming /= np.linalg.norm(incoming)
            edge_resolution.append(
                np.linalg.inv(np.column_stack((outgoing, incoming)))
            )

        for node in range(4):
            next_node = (node + 1) % 4
            edge_length = np.linalg.norm(
                self.local_coords_2d[next_node]
                - self.local_coords_2d[node]
            )
            force_at_node = edge_resolution[node][0]
            force_at_next = edge_resolution[next_node][1]
            for component in range(2):
                from_node = (
                    -heights[node] * force_at_node[component] / edge_length
                )
                from_next = (
                    heights[next_node]
                    * force_at_next[component]
                    / edge_length
                )
                row_node = 6 * node + 2
                row_next = 6 * next_node + 2
                column_node = 6 * node + component
                column_next = 6 * next_node + component
                warping[row_node, column_node] += from_node
                warping[row_node, column_next] += from_next
                warping[row_next, column_node] -= from_node
                warping[row_next, column_next] -= from_next

        # A moment initially tangent to the mean plane receives a local Mz
        # component so that it remains tangent to the two physical edges at
        # its node.  Store the coefficients because their total Mz must be
        # balanced by the self-equilibrating in-plane pattern below.
        rotation = self.T[:3, :3]
        normal_moment_coefficients = []
        for node in range(4):
            outgoing = self.coords[(node + 1) % 4] - self.coords[node]
            incoming = self.coords[(node - 1) % 4] - self.coords[node]
            nodal_normal = np.cross(outgoing, incoming)
            nodal_normal /= np.linalg.norm(nodal_normal)
            if np.dot(nodal_normal, self.local_z_axis) < 0.0:
                nodal_normal *= -1.0
            local_normal = rotation @ nodal_normal
            coefficients = -local_normal[:2] / local_normal[2]
            normal_moment_coefficients.append(coefficients)
            warping[
                6 * node + 5, 6 * node + 3 : 6 * node + 5
            ] += coefficients

        # Equal x/y forces arranged tangentially around the four corners have
        # zero resultant and supply only a moment about the mean-plane normal.
        force_pattern = np.array(
            [[-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0], [1.0, 1.0]]
        )
        pattern_moment = sum(
            x * force[1] - y * force[0]
            for (x, y), force in zip(self.local_coords_2d, force_pattern)
        )
        for force_node in range(4):
            for moment_node in range(4):
                warping[
                    6 * force_node : 6 * force_node + 2,
                    6 * moment_node + 3 : 6 * moment_node + 5,
                ] += np.outer(
                    force_pattern[force_node],
                    -normal_moment_coefficients[moment_node] / pattern_moment,
                )
        return warping

    def get_stiffness_matrix(self) -> np.ndarray:
        """Return the warped and globally oriented 24x24 stiffness matrix."""
        mean_plane_stiffness = self._build_local_stiffness()
        warping = self._build_warping_force_transformation()
        warped_stiffness = warping @ mean_plane_stiffness @ warping.T
        stiffness = self.T.T @ warped_stiffness @ self.T
        return 0.5 * (stiffness + stiffness.T)
