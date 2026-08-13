"""MSC/Nastran-style QUAD4 stiffness formulation.

The in-plane and curvature blocks use a uniform-strain field plus the two
orthogonal bilinear modes used by QUAD4.  The geometric hourglass metric is
built from the two element mid-lines, which is essential for skewed and
distorted quadrilaterals.  Non-coplanar nodes receive an equilibrium-based
force and moment correction after the mean-plane stiffness is formed.
"""

import numpy as np

from .cquad4_element import CQUAD4ElementBase


class CQUAD4NastranElementSolver(CQUAD4ElementBase):
    """Four-node shell reconstructing the Nastran QUAD4 stiffness.

    Flat elements use uniform membrane/curvature strains plus the two
    bilinear stabilization modes.  Warped elements additionally move the
    mean-plane forces and moments to the four physical nodes before the
    usual element-to-global coordinate transformation.
    """

    formulation = "nastran"

    @staticmethod
    def _selective_shear_sampling_points():
        gauss = 1.0 / np.sqrt(3.0)
        return (
            (0.0, -gauss, 0),
            (0.0, gauss, 0),
            (-gauss, 0.0, 1),
            (gauss, 0.0, 1),
        )

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
        """Assemble membrane, bending-flexibility and drilling blocks."""
        stiffness = np.zeros((24, 24))
        thickness = self.prop.t
        bending_inertia = (
            thickness**3 / 12.0 * self._get_bending_inertia_ratio()
        )

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
        shear = self._build_augmented_shear_stiffness(bending_inertia)

        membrane_dofs = [0, 1, 6, 7, 12, 13, 18, 19]
        plate_dofs = [2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22]
        self._scatter_block(stiffness, membrane, membrane_dofs)
        if thickness > 0.0:
            self._scatter_block(stiffness, bending + shear, plate_dofs)
        self._add_drilling_stiffness(stiffness)
        return 0.5 * (stiffness + stiffness.T)

    def _build_augmented_shear_stiffness(
        self, bending_inertia: float
    ) -> np.ndarray:
        """Build MacNeal QUAD4 shear stiffness from augmented flexibility.

        Four transverse strains are sampled at the selective-integration
        points: gamma_x at (0,+/-1/sqrt(3)) and gamma_y at
        (+/-1/sqrt(3),0).  Their flexibility is the sum of MacNeal's
        residual bending flexibility and, when MID3 is present, the physical
        transverse-shear flexibility.  This is intentionally a flexibility
        formulation: adding MID3 makes the element softer rather than adding
        a Reissner-Mindlin stiffness in parallel.
        """
        sampling_points = self._selective_shear_sampling_points()
        shear_operator = np.zeros((4, 12))
        jacobian_determinants = np.zeros(4)
        for row, (xi, eta, component) in enumerate(sampling_points):
            _, _, B_shear, determinant = self._get_B_matrices(xi, eta)
            shear_operator[row] = B_shear[component]
            jacobian_determinants[row] = determinant

        coordinates = self.local_coords_2d
        delta_x = 0.5 * abs(
            coordinates[1, 0]
            + coordinates[2, 0]
            - coordinates[0, 0]
            - coordinates[3, 0]
        )
        delta_y = 0.5 * abs(
            coordinates[3, 1]
            + coordinates[2, 1]
            - coordinates[0, 1]
            - coordinates[1, 1]
        )
        area = 0.5 * abs(
            np.dot(coordinates[:, 0], np.roll(coordinates[:, 1], -1))
            - np.dot(coordinates[:, 1], np.roll(coordinates[:, 0], -1))
        )

        bending_constitutive = self._get_D_bending() * bending_inertia
        if bending_constitutive[0, 0] <= 0.0:
            return np.zeros((12, 12))
        # Released MSC QUAD4 uses epsilon=0.025 in MacNeal's aspect-ratio
        # correction.  Cases 65/69 and 67/72 recover this value directly.
        epsilon = 0.025
        alpha = epsilon / (
            epsilon + (1.0 - epsilon) * (delta_x / delta_y) ** 2
        )
        beta = epsilon / (
            epsilon + (1.0 - epsilon) * (delta_y / delta_x) ** 2
        )
        flexibility_x = (
            delta_x**2 / (12.0 * area * bending_constitutive[0, 0])
        )
        flexibility_y = (
            delta_y**2 / (12.0 * area * bending_constitutive[1, 1])
        )
        residual_flexibility = np.zeros((4, 4))
        residual_flexibility[:2, :2] = flexibility_x * np.array(
            [[1.0 + alpha, 1.0 - alpha],
             [1.0 - alpha, 1.0 + alpha]]
        )
        residual_flexibility[2:, 2:] = flexibility_y * np.array(
            [[1.0 + beta, 1.0 - beta],
             [1.0 - beta, 1.0 + beta]]
        )

        total_flexibility = residual_flexibility
        if self.mat_shear is not None:
            shear_modulus = getattr(
                self.mat_shear,
                "G",
                self.mat_shear.E / (2.0 * (1.0 + self.mat_shear.nu)),
            )
            shear_thickness = (
                self.prop.t * self._get_shear_thickness_ratio()
            )
            physical_flexibility = 1.0 / (
                2.0
                * jacobian_determinants
                * shear_thickness
                * shear_modulus
            )
            total_flexibility = residual_flexibility + np.diag(
                physical_flexibility
            )

        # Solve Z X = D instead of explicitly forming Z^-1.  This preserves
        # the MacNeal expression D.T @ Z^-1 @ D with better conditioning.
        flexibility_solution = np.linalg.solve(
            total_flexibility, shear_operator
        )
        return shear_operator.T @ flexibility_solution

    def _recover_transverse_shear_forces(
        self, bending_displacements: np.ndarray
    ) -> np.ndarray:
        """Recover centroidal Qx/Qy from the augmented-flexibility modes."""
        bending_inertia = (
            self.prop.t**3 / 12.0 * self._get_bending_inertia_ratio()
        )
        shear_stiffness = self._build_augmented_shear_stiffness(
            bending_inertia
        )
        nodal_shear_forces = shear_stiffness @ bending_displacements

        shear_operator = np.zeros((4, 12))
        jacobian_determinants = np.zeros(4)
        for row, (xi, eta, component) in enumerate(
            self._selective_shear_sampling_points()
        ):
            _, _, B_shear, determinant = self._get_B_matrices(xi, eta)
            shear_operator[row] = B_shear[component]
            jacobian_determinants[row] = determinant

        generalized_forces = np.linalg.lstsq(
            shear_operator.T, nodal_shear_forces, rcond=None
        )[0]
        # Each modal force contains the integration measure 2*det(J) of its
        # selective sampling line. Remove that measure before averaging the
        # two samples of each centroidal force resultant.
        sampled_resultants = generalized_forces / (
            2.0 * jacobian_determinants
        )
        return np.array(
            [
                0.5 * (sampled_resultants[0] + sampled_resultants[1]),
                0.5 * (sampled_resultants[2] + sampled_resultants[3]),
            ]
        )

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

    def _get_recovery_displacements(
        self, element_displacements_global: np.ndarray
    ) -> np.ndarray:
        """Return work-conjugate displacements on the CQUAD mean plane."""
        warping = self._build_warping_force_transformation()
        return warping.T @ self.T @ element_displacements_global
