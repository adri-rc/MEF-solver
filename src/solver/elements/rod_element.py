"""Linear CROD element formulation."""

import numpy as np

from src.preprocessing.data_model import IsotropicMaterial, Node, RodProperty

from ..quadrature import GaussQuadrature
from ..transformations import CoordinateTransformation


class RodElementSolver:
    """Compute linear stiffness and stress recovery for a two-node CROD.

    Each node exposes six degrees of freedom, but only local axial translation
    and local torsional rotation contribute to this element.
    """

    AXIAL_DOFS = (0, 6)
    TORSIONAL_DOFS = (3, 9)

    def __init__(
        self,
        node1: Node,
        node2: Node,
        prop: RodProperty,
        mat: IsotropicMaterial,
        n_gaussp: int = 1,
    ):
        self.node1 = node1
        self.node2 = node2
        self.prop = prop
        self.mat = mat
        self.n_gaussp = n_gaussp

        self.n1_coords = node1.get_coords()
        self.n2_coords = node2.get_coords()
        self.L = np.linalg.norm(self.n2_coords - self.n1_coords)
        if self.L < 1.0e-12:
            raise ValueError(
                f"Zero length CROD between nodes {node1.nid} and {node2.nid}."
            )

        self._compute_transformation()

    def _compute_transformation(self) -> None:
        """Build the 12x12 element-to-global transformation matrix."""
        vx, vy, vz = CoordinateTransformation.compute_local_axes(
            self.n1_coords,
            self.n2_coords,
            orientation_vector=None,
        )
        rotation = CoordinateTransformation.get_rotation_matrix_3x3(vx, vy, vz)
        self.T = CoordinateTransformation.build_full_transformation_matrix(
            rotation,
            nodes_count=2,
            dofs_per_node=6,
        )

    @staticmethod
    def _get_b_matrix(jacobian: float) -> np.ndarray:
        """Return the constant strain-displacement operator."""
        shape_derivatives = np.array([-0.5, 0.5])
        return shape_derivatives / jacobian

    def get_stiffness_matrix(self) -> np.ndarray:
        """Return the 12x12 element stiffness matrix in basic coordinates."""
        points, weights = GaussQuadrature.get_rule(self.n_gaussp)
        axial_stiffness = np.zeros((2, 2))
        torsional_stiffness = np.zeros((2, 2))

        jacobian = self.L / 2.0
        for _, weight in zip(points, weights):
            b_matrix = self._get_b_matrix(jacobian)
            outer_b = np.outer(b_matrix, b_matrix) * jacobian * weight
            axial_stiffness += outer_b * self.mat.E * self.prop.A
            torsional_stiffness += outer_b * self.mat.G * self.prop.J

        local_stiffness = self._assemble_local_matrix(
            axial_stiffness,
            torsional_stiffness,
        )
        return self.T.T @ local_stiffness @ self.T

    def _assemble_local_matrix(
        self,
        axial_stiffness: np.ndarray,
        torsional_stiffness: np.ndarray,
    ) -> np.ndarray:
        """Embed the axial and torsional 2x2 blocks in the 12x12 matrix."""
        stiffness = np.zeros((12, 12))
        stiffness[np.ix_(self.AXIAL_DOFS, self.AXIAL_DOFS)] = axial_stiffness
        stiffness[np.ix_(self.TORSIONAL_DOFS, self.TORSIONAL_DOFS)] = (
            torsional_stiffness
        )
        return stiffness

    def calculate_stresses(
        self,
        u_global_vector: np.ndarray,
        dof_mapper,
    ) -> dict:
        """Recover axial and torsional response from basic displacements."""
        node1_dofs = dof_mapper.get_dof_indices(self.node1.nid)
        node2_dofs = dof_mapper.get_dof_indices(self.node2.nid)
        element_displacements = np.concatenate(
            [u_global_vector[node1_dofs], u_global_vector[node2_dofs]]
        )
        local_displacements = self.T @ element_displacements

        axial_strain = (
            local_displacements[self.AXIAL_DOFS[1]]
            - local_displacements[self.AXIAL_DOFS[0]]
        ) / self.L
        axial_stress = self.mat.E * axial_strain
        axial_force = axial_stress * self.prop.A

        twist_rate = (
            local_displacements[self.TORSIONAL_DOFS[1]]
            - local_displacements[self.TORSIONAL_DOFS[0]]
        ) / self.L
        torque = self.mat.G * self.prop.J * twist_rate

        outer_fiber_distance = getattr(self.prop, "C", 0.0)
        torsional_stress = 0.0
        if self.prop.J > 1.0e-12:
            torsional_stress = torque * outer_fiber_distance / self.prop.J

        return {
            "axial_force": axial_force,
            "axial_stress": axial_stress,
            "axial_strain": axial_strain,
            "torque": torque,
            "torsional_stress": torsional_stress,
            "twist_rate": twist_rate,
        }
