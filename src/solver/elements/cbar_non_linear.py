"""Corotational CBAR for geometrically nonlinear static analysis.

``PBARN1`` selects a linear/cubic Euler law (``LC``) or a one-point
linear-shear law (``LS``) over the elastic ``MAT1/PBAR`` section.  A
continuously rotating element frame removes rigid translation and rotation;
only the logarithmic axial stretch and the relative end rotations enter the
section law.  This is the SOL 400 ``PBARN1,SECT=S,BEH2=BAR`` scope: large
displacements and rotations with small elastic section strains.
"""

from dataclasses import dataclass, replace

import numpy as np
from scipy.spatial.transform import Rotation

from src.preprocessing.data_model import BarProperty, IsotropicMaterial, Node

from .bar_element import BarElementSolver
from ..transformations import CoordinateTransformation


@dataclass(frozen=True)
class CBARNonlinearState:
    """CBAR state at one total nodal-displacement vector."""

    tangent_stiffness: np.ndarray
    internal_force: np.ndarray
    current_length: float
    current_axes: np.ndarray
    basic_deformations: np.ndarray
    local_deformations: np.ndarray
    local_actions: np.ndarray
    axial_strain: float
    axial_stress: float
    axial_force: float
    twist_rate: float
    torque: float
    torsional_stress: float
    forces_n1: dict
    forces_n2: dict
    stress_points_end_a: dict
    stress_points_end_b: dict
    max_total_stress: float
    min_total_stress: float
    max_element_stress: float
    min_element_stress: float
    strain_energy: float

    def as_result_dict(self) -> dict:
        return {
            "type": "CBAR",
            "axial_force": self.axial_force,
            "axial_stress": self.axial_stress,
            "axial_strain": self.axial_strain,
            "torque": self.torque,
            "torsional_stress": self.torsional_stress,
            "twist_rate": self.twist_rate,
            "forces_n1": dict(self.forces_n1),
            "forces_n2": dict(self.forces_n2),
            "stress_points": dict(self.stress_points_end_a),
            "stress_points_end_a": dict(self.stress_points_end_a),
            "stress_points_end_b": dict(self.stress_points_end_b),
            "max_total_stress": self.max_total_stress,
            "min_total_stress": self.min_total_stress,
            "max_element_stress": self.max_element_stress,
            "min_element_stress": self.min_element_stress,
            "current_length": self.current_length,
            "strain_energy": self.strain_energy,
        }


class CBARNonlinear:
    """Two-node three-dimensional corotational CBAR.

    Nodal rotations are total rotation vectors in basic axes.  The nodal
    triads are obtained with the exponential map and compared with an element
    frame whose first axis is the current chord.  The PBARN1-selected section
    stiffness then acts on six deformational measures embedded in a 12-DOF
    local vector.  The tangent is the numerical differential of the complete
    spatial internal-force mapping, including frame and geometric-force terms.
    It need not be symmetric when total rotation vectors are used as the
    Newton coordinates.
    """

    DOFS_PER_NODE = 6
    ELEMENT_DOFS = 12
    _LENGTH_TOLERANCE = 1.0e-12
    _FRAME_TOLERANCE = 1.0e-10

    def __init__(
        self,
        node1: Node,
        node2: Node,
        prop: BarProperty,
        mat: IsotropicMaterial,
        orientation,
        *,
        wa=None,
        wb=None,
        pa: int = 0,
        pb: int = 0,
        large_displacements: bool = True,
    ):
        self.node1 = node1
        self.node2 = node2
        self.prop = prop
        self.mat = mat
        self.large_displacements = bool(large_displacements)
        self.wa = np.zeros(3) if wa is None else np.asarray(wa, dtype=float)
        self.wb = np.zeros(3) if wb is None else np.asarray(wb, dtype=float)
        self.pa = int(pa or 0)
        self.pb = int(pb or 0)

        self.reference_grid_a = np.asarray(node1.get_coords(), dtype=float)
        self.reference_grid_b = np.asarray(node2.get_coords(), dtype=float)
        orientation = np.asarray(orientation, dtype=float)
        if orientation.shape != (3,) or not np.all(np.isfinite(orientation)):
            raise ValueError("CBAR nonlinear orientation requires three finite values.")
        if self.wa.shape != (3,) or self.wb.shape != (3,):
            raise ValueError("CBAR nonlinear offsets require three values per end.")

        reference_end_a = self.reference_grid_a + self.wa
        reference_end_b = self.reference_grid_b + self.wb
        reference_chord = reference_end_b - reference_end_a
        self.reference_length = float(np.linalg.norm(reference_chord))
        if self.reference_length < self._LENGTH_TOLERANCE:
            raise ValueError(
                f"CBAR {node1.nid}-{node2.nid} has zero reference length."
            )

        ex, ey, ez = CoordinateTransformation.compute_local_axes(
            reference_end_a, reference_end_b, orientation_vector=orientation
        )
        self.reference_axes = np.column_stack((ex, ey, ez))

        # PBARN1 selects the interpolation used by the SOL 400 element.  LC
        # is the linear/cubic (Euler--Bernoulli) variant and deliberately
        # ignores the PBAR shear-area factors.  LS is a one-point
        # linear-shear element; it is not the closed-form Timoshenko matrix
        # used by the conventional linear CBAR path.
        euler_property = replace(prop, K1=0.0, K2=0.0)
        linear_solver = BarElementSolver(
            node1,
            node2,
            euler_property,
            mat,
            orientation=orientation,
            pa=self.pa,
            pb=self.pb,
        )
        if prop.nonlinear_integration == "LS":
            self.local_stiffness = self._linear_shear_local_stiffness()
        else:
            self.local_stiffness = linear_solver.get_local_stiffness_matrix()
        self._linear_global_stiffness = self._global_stiffness_from_local(
            self.local_stiffness
        )

    @staticmethod
    def _add_block(stiffness: np.ndarray, block: np.ndarray, indices) -> None:
        stiffness[np.ix_(indices, indices)] += block

    def _apply_end_releases(self, stiffness: np.ndarray) -> np.ndarray:
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
        condensed = k_rr - k_rq @ np.linalg.pinv(k_qq, rcond=1.0e-12) @ k_rq.T
        result = np.zeros_like(stiffness)
        result[np.ix_(retained, retained)] = condensed
        return 0.5 * (result + result.T)

    def _linear_shear_local_stiffness(self) -> np.ndarray:
        """Return MSC's PBARN1 ``INT2=LS`` one-point element matrix.

        Transverse displacement and section rotation are linear.  Curvature
        is therefore constant and the shear strain is evaluated at the
        element midpoint.  A zero/blank PBAR shear factor is interpreted as
        the SOL 400 unit default, as isolated by validation cases 159/160.
        """
        length = self.reference_length
        stiffness = np.zeros((12, 12))
        pair = np.array([[1.0, -1.0], [-1.0, 1.0]]) / length
        self._add_block(
            stiffness,
            self.mat.E * self.prop.A * pair,
            BarElementSolver.IDX_AXIAL,
        )
        self._add_block(
            stiffness,
            self.mat.G * self.prop.J * pair,
            BarElementSolver.IDX_TORSION,
        )

        curvature_xy = np.array([0.0, -1.0 / length, 0.0, 1.0 / length])
        shear_xy = np.array([-1.0 / length, -0.5, 1.0 / length, -0.5])
        rotation_signs = np.array([1.0, -1.0, 1.0, -1.0])
        curvature_xz = rotation_signs * curvature_xy
        shear_xz = rotation_signs * shear_xy
        shear_factor_1 = self.prop.K1 if self.prop.K1 > 0.0 else 1.0
        shear_factor_2 = self.prop.K2 if self.prop.K2 > 0.0 else 1.0

        bending_xy = length * (
            self.mat.E * self.prop.I1 * np.outer(curvature_xy, curvature_xy)
            + shear_factor_1
            * self.prop.A
            * self.mat.G
            * np.outer(shear_xy, shear_xy)
        )
        bending_xz = length * (
            self.mat.E * self.prop.I2 * np.outer(curvature_xz, curvature_xz)
            + shear_factor_2
            * self.prop.A
            * self.mat.G
            * np.outer(shear_xz, shear_xz)
        )
        self._add_block(
            stiffness, bending_xy, BarElementSolver.IDX_BEND_XY
        )
        self._add_block(
            stiffness, bending_xz, BarElementSolver.IDX_BEND_XZ
        )

        if self.prop.I12 != 0.0:
            coupling = (
                self.mat.E
                * self.prop.I12
                * length
                * np.outer(curvature_xy, curvature_xz)
            )
            xy = BarElementSolver.IDX_BEND_XY
            xz = BarElementSolver.IDX_BEND_XZ
            stiffness[np.ix_(xy, xz)] += coupling
            stiffness[np.ix_(xz, xy)] += coupling.T

        return self._apply_end_releases(stiffness)

    def _global_stiffness_from_local(self, local_stiffness: np.ndarray) -> np.ndarray:
        local_transform = self._block_rotation(self.reference_axes).T
        offset_transform = np.eye(12)
        offset_transform[0:3, 3:6] = -self._skew(self.wa)
        offset_transform[6:9, 9:12] = -self._skew(self.wb)
        k_reference_ends = local_transform.T @ local_stiffness @ local_transform
        return offset_transform.T @ k_reference_ends @ offset_transform

    @staticmethod
    def _block_rotation(axes: np.ndarray) -> np.ndarray:
        """Map local end forces/moments to basic axes."""
        transform = np.zeros((12, 12))
        for start in (0, 3, 6, 9):
            transform[start:start + 3, start:start + 3] = axes
        return transform

    @staticmethod
    def _rotation_matrix(rotation_vector: np.ndarray) -> np.ndarray:
        return Rotation.from_rotvec(np.asarray(rotation_vector, dtype=float)).as_matrix()

    @staticmethod
    def _rotation_vector(rotation_matrix: np.ndarray) -> np.ndarray:
        return Rotation.from_matrix(rotation_matrix).as_rotvec()

    def _corotated_axes(
        self,
        current_chord: np.ndarray,
        nodal_axes_a: np.ndarray,
        nodal_axes_b: np.ndarray,
    ) -> np.ndarray:
        """Build a chord-following frame with the mean physical end twist.

        Bending rotations must not contaminate the section roll.  The frame is
        therefore obtained by the shortest rotation of the reference chord,
        followed by the mean nodal twist about the current chord.  This is the
        standard swing--twist construction used by corotational beams.
        """
        reference_x = self.reference_axes[:, 0]
        current_x = current_chord / np.linalg.norm(current_chord)
        cross = np.cross(reference_x, current_x)
        sine = float(np.linalg.norm(cross))
        cosine = float(np.clip(reference_x @ current_x, -1.0, 1.0))
        if cosine < -1.0 + self._FRAME_TOLERANCE:
            swing = Rotation.from_rotvec(np.pi * self.reference_axes[:, 1]).as_matrix()
        elif sine < np.finfo(float).eps:
            swing = np.eye(3)
        else:
            cross_matrix = self._skew(cross)
            swing = (
                np.eye(3)
                + cross_matrix
                + (cross_matrix @ cross_matrix) / (1.0 + cosine)
            )

        zero_twist_axes = swing @ self.reference_axes
        twist_a = self._rotation_vector(zero_twist_axes.T @ nodal_axes_a)[0]
        twist_b = self._rotation_vector(zero_twist_axes.T @ nodal_axes_b)[0]
        mean_twist = 0.5 * (twist_a + twist_b)
        twist = Rotation.from_rotvec([mean_twist, 0.0, 0.0]).as_matrix()
        return zero_twist_axes @ twist

    def _nonlinear_force_data(self, displacement: np.ndarray) -> dict:
        translation_a = displacement[0:3]
        rotation_a = displacement[3:6]
        translation_b = displacement[6:9]
        rotation_b = displacement[9:12]

        rotation_matrix_a = self._rotation_matrix(rotation_a)
        rotation_matrix_b = self._rotation_matrix(rotation_b)
        offset_a = rotation_matrix_a @ self.wa
        offset_b = rotation_matrix_b @ self.wb
        end_a = self.reference_grid_a + translation_a + offset_a
        end_b = self.reference_grid_b + translation_b + offset_b
        chord = end_b - end_a
        current_length = float(np.linalg.norm(chord))
        if current_length < self._LENGTH_TOLERANCE:
            raise ValueError(
                f"CBAR {self.node1.nid}-{self.node2.nid} has collapsed current length."
            )

        nodal_axes_a = rotation_matrix_a @ self.reference_axes
        nodal_axes_b = rotation_matrix_b @ self.reference_axes
        current_axes = self._corotated_axes(chord, nodal_axes_a, nodal_axes_b)
        local_rotation_a = self._rotation_vector(current_axes.T @ nodal_axes_a)
        local_rotation_b = self._rotation_vector(current_axes.T @ nodal_axes_b)

        axial_measure = self.reference_length * np.log(
            current_length / self.reference_length
        )
        local_deformations = np.zeros(12)
        local_deformations[0] = -0.5 * axial_measure
        local_deformations[6] = 0.5 * axial_measure
        local_deformations[3:6] = local_rotation_a
        local_deformations[9:12] = local_rotation_b
        local_actions = self.local_stiffness @ local_deformations

        basic_deformations = np.array(
            [
                axial_measure,
                local_rotation_b[0] - local_rotation_a[0],
                local_rotation_a[1],
                local_rotation_a[2],
                local_rotation_b[1],
                local_rotation_b[2],
            ]
        )

        force_basic = self._block_rotation(current_axes) @ local_actions
        if np.any(self.wa):
            force_basic[3:6] += np.cross(offset_a, force_basic[0:3])
        if np.any(self.wb):
            force_basic[9:12] += np.cross(offset_b, force_basic[6:9])

        return {
            "internal_force": force_basic,
            "current_length": current_length,
            "current_axes": current_axes,
            "basic_deformations": basic_deformations,
            "local_deformations": local_deformations,
            "local_actions": local_actions,
            "strain_energy": 0.5 * float(
                local_deformations @ self.local_stiffness @ local_deformations
            ),
        }

    def evaluate_internal_force(self, element_displacements_basic: np.ndarray) -> np.ndarray:
        """Return the total element internal force in basic axes."""
        displacement = self._validated_displacement(element_displacements_basic)
        if not self.large_displacements:
            return self._linear_global_stiffness @ displacement
        data = self._nonlinear_force_data(displacement)
        if np.linalg.norm(data["local_deformations"], ord=np.inf) < 1.0e-12:
            return np.zeros(12)
        return data["internal_force"]

    def _differentiation_steps(self, displacement: np.ndarray) -> np.ndarray:
        translation_scale = max(self.reference_length, 1.0)
        steps = np.zeros(12)
        for column in range(12):
            local_component = column % 6
            scale = translation_scale if local_component < 3 else 1.0
            scale = max(scale, abs(float(displacement[column])))
            steps[column] = 2.0e-6 * scale
        return steps

    def _numerical_tangent(self, displacement: np.ndarray) -> np.ndarray:
        steps = self._differentiation_steps(displacement)
        tangent = np.zeros((12, 12))
        for column, step in enumerate(steps):
            plus = displacement.copy()
            minus = displacement.copy()
            plus[column] += step
            minus[column] -= step
            tangent[:, column] = (
                self._nonlinear_force_data(plus)["internal_force"]
                - self._nonlinear_force_data(minus)["internal_force"]
            ) / (2.0 * step)
        return tangent

    @staticmethod
    def _force_dict(values: np.ndarray) -> dict:
        return dict(zip(("Fx", "Fy", "Fz", "Mx", "My", "Mz"), values))

    def _stress_recovery(self, data: dict) -> dict:
        local_actions = data["local_actions"]
        extension = data["basic_deformations"][0]
        axial_strain = extension / self.reference_length
        axial_stress = self.mat.E * axial_strain
        axial_force = 0.5 * (local_actions[6] - local_actions[0])
        twist = data["basic_deformations"][1]
        twist_rate = twist / self.reference_length
        torque = 0.5 * (local_actions[9] - local_actions[3])

        points = self.prop.stress_points or [(0.0, 0.0)]
        inertia = np.array(
            [[self.prop.I1, self.prop.I12], [self.prop.I12, self.prop.I2]],
            dtype=float,
        )
        end_moments = (
            (-local_actions[5], -local_actions[4]),
            (local_actions[11], local_actions[10]),
        )
        stress_sets = []
        all_stresses = []
        maximum_torsional_stress = 0.0
        for moment_z, moment_y in end_moments:
            coefficients = np.linalg.solve(inertia, np.array([-moment_z, moment_y]))
            result = {}
            for index, (y, z) in enumerate(points):
                bending = float(coefficients @ np.array([y, z]))
                total = axial_stress + bending
                radius = float(np.hypot(y, z))
                torsional = torque * radius / self.prop.J if self.prop.J else 0.0
                maximum_torsional_stress = max(
                    maximum_torsional_stress, abs(torsional)
                )
                result[f"P{index}"] = {
                    "sigma_total": total,
                    "sigma_bending": bending,
                    "sigma_axial": axial_stress,
                    "tau_xy": torsional,
                    "y": y,
                    "z": z,
                }
                all_stresses.append(total)
            stress_sets.append(result)

        end_a_stresses = [value["sigma_total"] for value in stress_sets[0].values()]
        return {
            "axial_strain": axial_strain,
            "axial_stress": axial_stress,
            "axial_force": axial_force,
            "twist_rate": twist_rate,
            "torque": torque,
            "torsional_stress": maximum_torsional_stress,
            "forces_n1": self._force_dict(local_actions[:6]),
            "forces_n2": self._force_dict(local_actions[6:]),
            "stress_points_end_a": stress_sets[0],
            "stress_points_end_b": stress_sets[1],
            "max_total_stress": max(end_a_stresses),
            "min_total_stress": min(end_a_stresses),
            "max_element_stress": max(all_stresses),
            "min_element_stress": min(all_stresses),
        }

    def evaluate(self, element_displacements_basic: np.ndarray) -> CBARNonlinearState:
        displacement = self._validated_displacement(element_displacements_basic)
        if self.large_displacements:
            data = self._nonlinear_force_data(displacement)
            internal_force = data["internal_force"].copy()
            tangent = self._numerical_tangent(displacement)
            if np.linalg.norm(data["local_deformations"], ord=np.inf) < 1.0e-12:
                internal_force.fill(0.0)
            data["internal_force"] = internal_force
        else:
            data = self._linear_force_data(displacement)
            tangent = self._linear_global_stiffness.copy()
        recovered = self._stress_recovery(data)
        return CBARNonlinearState(
            tangent_stiffness=tangent,
            internal_force=data["internal_force"],
            current_length=data["current_length"],
            current_axes=data["current_axes"],
            basic_deformations=data["basic_deformations"],
            local_deformations=data["local_deformations"],
            local_actions=data["local_actions"],
            strain_energy=data["strain_energy"],
            **recovered,
        )

    def _linear_force_data(self, displacement: np.ndarray) -> dict:
        # _block_rotation expects local-to-basic axes.  Its transpose therefore
        # maps basic vectors to the local reference frame.
        local_transform = self._block_rotation(self.reference_axes).T
        offset_transform = np.eye(12)
        offset_transform[0:3, 3:6] = -self._skew(self.wa)
        offset_transform[6:9, 9:12] = -self._skew(self.wb)
        local_deformations = local_transform @ offset_transform @ displacement
        local_actions = self.local_stiffness @ local_deformations
        current_chord = (
            self.reference_grid_b + displacement[6:9] + self.wb
            - self.reference_grid_a - displacement[0:3] - self.wa
        )
        current_length = float(np.linalg.norm(current_chord))
        basic_deformations = np.array(
            [
                local_deformations[6] - local_deformations[0],
                local_deformations[9] - local_deformations[3],
                local_deformations[4],
                local_deformations[5],
                local_deformations[10],
                local_deformations[11],
            ]
        )
        return {
            "internal_force": self._linear_global_stiffness @ displacement,
            "current_length": current_length,
            "current_axes": self.reference_axes,
            "basic_deformations": basic_deformations,
            "local_deformations": local_deformations,
            "local_actions": local_actions,
            "strain_energy": 0.5 * float(
                displacement @ self._linear_global_stiffness @ displacement
            ),
        }

    @staticmethod
    def _skew(vector: np.ndarray) -> np.ndarray:
        x, y, z = vector
        return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])

    def _validated_displacement(self, displacement) -> np.ndarray:
        vector = np.asarray(displacement, dtype=float)
        if vector.shape != (self.ELEMENT_DOFS,):
            raise ValueError(
                "CBAR nonlinear evaluation requires a 12-component element "
                f"displacement vector, received shape {vector.shape}."
            )
        if not np.all(np.isfinite(vector)):
            raise ValueError("CBAR nonlinear displacement contains NaN/Inf values.")
        return vector

    def get_tangent_stiffness(self, displacement: np.ndarray) -> np.ndarray:
        return self.evaluate(displacement).tangent_stiffness

    def get_internal_force(self, displacement: np.ndarray) -> np.ndarray:
        return self.evaluate_internal_force(displacement)
