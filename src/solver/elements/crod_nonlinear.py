"""Corotational CROD formulation for geometric nonlinear analysis.

The axial law follows the MSC Nastran SOL 400 references in cases 94--122:
engineering strain based on the reference length, linear MAT1 stress, and
reference PROD area.  Geometry enters through the current chord direction and
the consistent differential (material plus geometric) stiffness.
"""

from dataclasses import dataclass

import numpy as np

from src.preprocessing.data_model import IsotropicMaterial, Node, RodProperty


@dataclass(frozen=True)
class CRODNonlinearState:
    """Element state evaluated at one total-displacement vector."""

    tangent_stiffness: np.ndarray
    internal_force: np.ndarray
    current_length: float
    current_direction: np.ndarray
    axial_strain: float
    axial_stress: float
    axial_force: float
    twist_rate: float
    torque: float
    torsional_stress: float
    strain_energy: float

    def as_result_dict(self) -> dict:
        return {
            "type": "CROD",
            "axial_force": self.axial_force,
            "axial_stress": self.axial_stress,
            "axial_strain": self.axial_strain,
            "torque": self.torque,
            "torsional_stress": self.torsional_stress,
            "twist_rate": self.twist_rate,
            "current_length": self.current_length,
            "strain_energy": self.strain_energy,
        }


class CRODNonlinear:
    """Two-node CROD with finite chord rotation and elastic axial strain.

    Translational terms use a total-Lagrangian reference length with a
    corotational current axis.  The derivative of the internal force is

    ``EA/L0 * n*n.T + N/l * (I - n*n.T)``.

    CROD torsion remains uncoupled and is measured about the reference axis.
    This reproduces cases 116/117, including the fact that axial extension
    does not change ``GJ/L0``. Finite-rotation beam kinematics are implemented
    separately by the nonlinear CBAR formulation and are not part of CROD.
    """

    DOFS_PER_NODE = 6
    ELEMENT_DOFS = 12
    _CURRENT_LENGTH_TOLERANCE = 1.0e-12

    def __init__(
        self,
        node1: Node,
        node2: Node,
        prop: RodProperty,
        mat: IsotropicMaterial,
        *,
        large_displacements: bool = True,
    ):
        self.node1 = node1
        self.node2 = node2
        self.prop = prop
        self.mat = mat
        self.large_displacements = bool(large_displacements)

        self.reference_position_1 = np.asarray(node1.get_coords(), dtype=float)
        self.reference_position_2 = np.asarray(node2.get_coords(), dtype=float)
        reference_chord = self.reference_position_2 - self.reference_position_1
        self.reference_length = float(np.linalg.norm(reference_chord))
        if self.reference_length < self._CURRENT_LENGTH_TOLERANCE:
            raise ValueError(
                f"CROD {node1.nid}-{node2.nid} has zero reference length."
            )
        self.reference_direction = reference_chord / self.reference_length

    def evaluate(self, element_displacements_basic: np.ndarray) -> CRODNonlinearState:
        """Return internal force and consistent tangent in the basic frame."""
        displacement = np.asarray(element_displacements_basic, dtype=float)
        if displacement.shape != (self.ELEMENT_DOFS,):
            raise ValueError(
                "CROD nonlinear evaluation requires a 12-component element "
                f"displacement vector, received shape {displacement.shape}."
            )
        if not np.all(np.isfinite(displacement)):
            raise ValueError("CROD nonlinear displacement contains NaN/Inf values.")

        translation_1 = displacement[0:3]
        rotation_1 = displacement[3:6]
        translation_2 = displacement[6:9]
        rotation_2 = displacement[9:12]

        current_chord = (
            self.reference_position_2
            + translation_2
            - self.reference_position_1
            - translation_1
        )
        actual_current_length = float(np.linalg.norm(current_chord))

        if self.large_displacements:
            if actual_current_length < self._CURRENT_LENGTH_TOLERANCE:
                raise ValueError(
                    f"CROD {self.node1.nid}-{self.node2.nid} has collapsed current length."
                )
            current_length = actual_current_length
            current_direction = current_chord / current_length
            axial_extension = current_length - self.reference_length
        else:
            current_length = actual_current_length
            current_direction = self.reference_direction
            axial_extension = float(
                self.reference_direction @ (translation_2 - translation_1)
            )

        axial_strain = axial_extension / self.reference_length
        axial_stress = self.mat.E * axial_strain
        axial_force = self.prop.A * axial_stress

        direction_outer = np.outer(current_direction, current_direction)
        axial_material_tangent = (
            self.mat.E * self.prop.A / self.reference_length
        ) * direction_outer
        if self.large_displacements:
            projector = np.eye(3) - direction_outer
            axial_geometric_tangent = (axial_force / current_length) * projector
        else:
            axial_geometric_tangent = np.zeros((3, 3))
        translational_tangent = axial_material_tangent + axial_geometric_tangent

        internal_force = np.zeros(self.ELEMENT_DOFS)
        tangent = np.zeros((self.ELEMENT_DOFS, self.ELEMENT_DOFS))
        force_vector = axial_force * current_direction
        internal_force[0:3] = -force_vector
        internal_force[6:9] = force_vector
        self._assemble_two_node_block(tangent, (0, 6), translational_tangent)

        relative_rotation = rotation_2 - rotation_1
        twist_angle = float(relative_rotation @ self.reference_direction)
        twist_rate = twist_angle / self.reference_length
        torque = self.mat.G * self.prop.J * twist_rate
        torsional_vector = torque * self.reference_direction
        internal_force[3:6] = -torsional_vector
        internal_force[9:12] = torsional_vector

        torsional_tangent = (
            self.mat.G * self.prop.J / self.reference_length
        ) * np.outer(self.reference_direction, self.reference_direction)
        self._assemble_two_node_block(tangent, (3, 9), torsional_tangent)

        torsional_stress = 0.0
        if self.prop.J > self._CURRENT_LENGTH_TOLERANCE:
            torsional_stress = torque * self.prop.C / self.prop.J

        axial_energy = 0.5 * (
            self.mat.E * self.prop.A / self.reference_length
        ) * axial_extension**2
        torsional_energy = 0.5 * (
            self.mat.G * self.prop.J / self.reference_length
        ) * twist_angle**2

        return CRODNonlinearState(
            tangent_stiffness=tangent,
            internal_force=internal_force,
            current_length=current_length,
            current_direction=current_direction,
            axial_strain=axial_strain,
            axial_stress=axial_stress,
            axial_force=axial_force,
            twist_rate=twist_rate,
            torque=torque,
            torsional_stress=torsional_stress,
            strain_energy=axial_energy + torsional_energy,
        )

    @staticmethod
    def _assemble_two_node_block(
        matrix: np.ndarray, starts: tuple[int, int], block: np.ndarray
    ) -> None:
        first, second = starts
        matrix[first:first + 3, first:first + 3] += block
        matrix[first:first + 3, second:second + 3] -= block
        matrix[second:second + 3, first:first + 3] -= block
        matrix[second:second + 3, second:second + 3] += block

