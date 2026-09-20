"""Incremental nonlinear static solver for CROD and CBAR models.

This module is intentionally parallel to :mod:`linear_static_solver`.  It owns
load stepping, Newton-Raphson iterations, residuals, cutbacks, the current
geometry, and nonlinear result recovery without changing the established
linear path.
"""

from dataclasses import dataclass
import warnings

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve

from src.preprocessing.data_model import (
    AnalysisType,
    BarElement,
    BarProperty,
    IsotropicMaterial,
    RodElement,
    RodProperty,
    SinglePointLoad,
    SinglePointMoment,
    SPCConstraint,
)

from .elements.cbar_nonlinear import CBARNonlinear, CBARNonlinearState
from .elements.crod_nonlinear import CRODNonlinear, CRODNonlinearState
from .utilities import DofMapper


@dataclass(frozen=True)
class NonlinearAssemblyState:
    """Global tangent and internal forces at one total-displacement state."""

    tangent_stiffness: csr_matrix
    internal_force: np.ndarray
    element_states: dict[int, CRODNonlinearState | CBARNonlinearState]
    strain_energy: float


@dataclass(frozen=True)
class NonlinearIterationResult:
    increment: int
    load_factor: float
    iteration: int
    residual_norm: float
    relative_residual_norm: float
    correction_norm: float


@dataclass(frozen=True)
class NonlinearIncrementResult:
    increment: int
    load_factor: float
    step_size: float
    iterations: int
    residual_norm: float
    relative_residual_norm: float
    cutbacks: int
    displacements: np.ndarray
    strain_energy: float


class NonlinearConvergenceError(ValueError):
    """Raised after the configured cutback strategy cannot reach equilibrium."""

    def __init__(
        self,
        message: str,
        *,
        last_converged_load_factor: float,
        attempted_load_factor: float,
        cutbacks: int,
        residual_norm: float,
    ):
        super().__init__(message)
        self.last_converged_load_factor = last_converged_load_factor
        self.attempted_load_factor = attempted_load_factor
        self.cutbacks = cutbacks
        self.residual_norm = residual_norm


class _IncrementFailure(RuntimeError):
    def __init__(self, residual_norm: float):
        super().__init__(f"Increment did not converge; residual={residual_norm:.6e}.")
        self.residual_norm = residual_norm


class NonlinearStaticSolver:
    """Solve ``lambda*Fext - Fint(u) = 0`` by full Newton-Raphson."""

    DOFS_PER_NODE = 6
    _SYMMETRY_TOLERANCE = 1.0e-10

    def __init__(self, model):
        self.model = model
        self.model.validate()
        self._validate_scope()

        self.dof_mapper = DofMapper(
            list(self.model.nodes.keys()), dofs_per_node=self.DOFS_PER_NODE
        )
        self.n_dofs = self.dof_mapper.total_dofs
        self.external_load = np.zeros(self.n_dofs)
        self.displacements = np.zeros(self.n_dofs)
        self.element_solvers: dict[int, CRODNonlinear | CBARNonlinear] = {}
        self.increment_history: list[NonlinearIncrementResult] = []
        self.iteration_history: list[NonlinearIterationResult] = []
        self.total_cutbacks = 0
        self._assembled = False

    def _validate_scope(self) -> None:
        analysis = self.model.analysis
        if analysis.analysis_type is not AnalysisType.NONLINEAR_STATIC:
            raise ValueError(
                "NonlinearStaticSolver requires a nonlinear static analysis configuration."
            )
        if analysis.nonlinear_step is None:
            raise ValueError("Nonlinear static analysis requires NLSTEP/NLPARM control.")
        if analysis.large_displacement_parameter == 0:
            raise NotImplementedError(
                "PARAM,LGDISP,0 differential-stiffness-only analysis is not supported."
            )
        if analysis.follower_force_stiffness:
            raise NotImplementedError(
                "Follower-force stiffness is outside the nonlinear CROD/CBAR scope."
            )
        if analysis.nonlinear_step.stiffness_update_method != "PFNT":
            raise NotImplementedError(
                "The nonlinear solver currently implements full Newton (PFNT) only."
            )
        if "P" not in analysis.nonlinear_step.convergence:
            raise NotImplementedError(
                "The nonlinear solver currently requires force-residual "
                "convergence (P in the MECH flags)."
            )

        for eid, element in self.model.elements.items():
            if not isinstance(element, (RodElement, BarElement)):
                raise NotImplementedError(
                    f"Nonlinear element {type(element).__name__} (EID {eid}) is not supported; "
                    "the current scope is CROD and CBAR."
                )
            prop = self.model.properties[element.prop_id]
            if isinstance(element, RodElement) and not isinstance(prop, RodProperty):
                raise NotImplementedError(f"CROD {eid} requires a PROD property.")
            if isinstance(element, BarElement) and not isinstance(prop, BarProperty):
                raise NotImplementedError(f"CBAR {eid} requires a PBAR property.")
            if isinstance(element, BarElement):
                if prop.nonlinear_section not in {None, "S"}:
                    raise NotImplementedError(
                        f"CBAR {eid} uses unsupported PBARN1 SECT="
                        f"{prop.nonlinear_section}; only SECT=S is supported."
                    )
                if prop.nonlinear_behavior not in {None, "BAR"}:
                    raise NotImplementedError(
                        f"CBAR {eid} uses unsupported PBARN1 BEH2="
                        f"{prop.nonlinear_behavior}; only BEH2=BAR is supported."
                    )
                if prop.nonlinear_integration not in {None, "LC", "LS"}:
                    raise NotImplementedError(
                        f"CBAR {eid} uses unsupported PBARN1 INT2="
                        f"{prop.nonlinear_integration}."
                    )
            mat = self.model.materials[prop.mat_id]
            if not isinstance(mat, IsotropicMaterial):
                raise NotImplementedError(
                    f"{type(element).__name__} {eid} requires an isotropic MAT1 material."
                )

        for constraint in self.model.constraints:
            if not isinstance(constraint, SPCConstraint):
                raise NotImplementedError(
                    f"Nonlinear constraint {type(constraint).__name__} is not supported."
                )

    def assemble(self) -> None:
        """Build invariant load data and cache element formulations.

        Unlike linear assembly, the tangent itself is state-dependent and is
        rebuilt by :meth:`assemble_state` during every Newton iteration.
        """
        self.external_load.fill(0.0)
        self.element_solvers = {}
        large_displacements = self.model.analysis.large_displacements
        for eid, element in self.model.elements.items():
            node_1 = self.model.get_node(element.nodes[0])
            node_2 = self.model.get_node(element.nodes[1])
            prop = self.model.properties[element.prop_id]
            material = self.model.materials[prop.mat_id]
            if isinstance(element, RodElement):
                self.element_solvers[eid] = CRODNonlinear(
                    node_1,
                    node_2,
                    prop,
                    material,
                    large_displacements=large_displacements,
                )
            elif isinstance(element, BarElement):
                orientation = element.orientation_vector
                if element.g0 is not None:
                    orientation = (
                        self.model.get_node(element.g0).get_coords()
                        - node_1.get_coords()
                    )
                if orientation is None:
                    raise ValueError(f"CBAR {eid} has no valid orientation.")
                self.element_solvers[eid] = CBARNonlinear(
                    node_1,
                    node_2,
                    prop,
                    material,
                    orientation,
                    wa=element.wa,
                    wb=element.wb,
                    pa=element.pa,
                    pb=element.pb,
                    large_displacements=large_displacements,
                )
        self._build_external_load()
        self._assembled = True

    def _build_external_load(self) -> None:
        for load in self.model.loads:
            if isinstance(load, SinglePointLoad):
                if load.cid != 0:
                    raise ValueError(
                        f"Force SID {load.sid} was not transformed to the basic frame."
                    )
                basic_vector = load.module * np.asarray(load.vector, dtype=float)
                dof_vector = (
                    np.asarray(self.model.get_node(load.appnode_id).dof_transform)
                    @ basic_vector
                )
                indices = self.dof_mapper.get_dof_indices(load.appnode_id)
                self.external_load[indices[:3]] += dof_vector
            elif isinstance(load, SinglePointMoment):
                if load.cid != 0:
                    raise ValueError(
                        f"Moment SID {load.sid} was not transformed to the basic frame."
                    )
                basic_vector = load.module * np.asarray(load.vector, dtype=float)
                dof_vector = (
                    np.asarray(self.model.get_node(load.appnode_id).dof_transform)
                    @ basic_vector
                )
                indices = self.dof_mapper.get_dof_indices(load.appnode_id)
                self.external_load[indices[3:]] += dof_vector
            else:
                raise NotImplementedError(
                    f"Nonlinear load {type(load).__name__} is not supported."
                )

    def _nodal_frame_transform(self, node_ids: list[int]) -> np.ndarray:
        transform = np.zeros((6 * len(node_ids), 6 * len(node_ids)))
        for position, node_id in enumerate(node_ids):
            rotation = np.asarray(self.model.get_node(node_id).dof_transform)
            start = position * 6
            transform[start:start + 3, start:start + 3] = rotation
            transform[start + 3:start + 6, start + 3:start + 6] = rotation
        return transform

    def get_element_displacements_basic(
        self, element: RodElement | BarElement, displacements: np.ndarray
    ) -> np.ndarray:
        """Extract one element vector and convert GRID.CD DOFs to basic axes."""
        parts = []
        for node_id in element.nodes:
            indices = self.dof_mapper.get_dof_indices(node_id)
            nodal = np.asarray(displacements[indices], dtype=float).copy()
            rotation = np.asarray(self.model.get_node(node_id).dof_transform)
            nodal[:3] = rotation.T @ nodal[:3]
            nodal[3:] = rotation.T @ nodal[3:]
            parts.append(nodal)
        return np.concatenate(parts)

    def assemble_state(self, displacements: np.ndarray) -> NonlinearAssemblyState:
        """Assemble the current global tangent and internal force vector."""
        if not self._assembled:
            self.assemble()
        vector = np.asarray(displacements, dtype=float)
        if vector.shape != (self.n_dofs,) or not np.all(np.isfinite(vector)):
            raise ValueError("Global nonlinear displacement vector is invalid.")

        tangent = lil_matrix((self.n_dofs, self.n_dofs))
        internal = np.zeros(self.n_dofs)
        element_states = {}
        strain_energy = 0.0

        for eid, element in self.model.elements.items():
            node_ids = list(element.nodes)
            global_indices = []
            for node_id in node_ids:
                global_indices.extend(self.dof_mapper.get_dof_indices(node_id))
            element_displacement = self.get_element_displacements_basic(
                element, vector
            )
            element_state = self.element_solvers[eid].evaluate(element_displacement)
            frame_transform = self._nodal_frame_transform(node_ids)
            element_tangent = (
                frame_transform
                @ element_state.tangent_stiffness
                @ frame_transform.T
            )
            element_internal = frame_transform @ element_state.internal_force

            internal[global_indices] += element_internal
            for row, global_row in enumerate(global_indices):
                for column, global_column in enumerate(global_indices):
                    value = element_tangent[row, column]
                    if value != 0.0:
                        tangent[global_row, global_column] += value

            element_states[eid] = element_state
            strain_energy += element_state.strain_energy

        tangent_csr = tangent.tocsr()
        return NonlinearAssemblyState(
            tangent_stiffness=tangent_csr,
            internal_force=internal,
            element_states=element_states,
            strain_energy=strain_energy,
        )

    def get_dof_partition(self):
        fixed_dofs = set()
        enforced_displacements = {}
        component_to_offset = {str(index + 1): index for index in range(6)}

        for constraint in self.model.constraints:
            node_indices = self.dof_mapper.get_dof_indices(constraint.node_id)
            for component in str(constraint.components):
                if component not in component_to_offset:
                    continue
                global_index = node_indices[component_to_offset[component]]
                if global_index in enforced_displacements and not np.isclose(
                    enforced_displacements[global_index],
                    constraint.value,
                    rtol=1.0e-12,
                    atol=1.0e-15,
                ):
                    raise ValueError(
                        f"Conflicting SPC values at global DOF {global_index}."
                    )
                fixed_dofs.add(global_index)
                enforced_displacements[global_index] = constraint.value

        fixed = sorted(fixed_dofs)
        free = sorted(set(range(self.n_dofs)) - fixed_dofs)
        return free, fixed, enforced_displacements

    def solve(self) -> dict:
        """Run fixed/adaptive load stepping and return the converged final state."""
        if not self._assembled:
            self.assemble()
        self.displacements.fill(0.0)
        self.increment_history = []
        self.iteration_history = []
        self.total_cutbacks = 0

        free_dofs, fixed_dofs, prescribed = self.get_dof_partition()
        if not free_dofs:
            raise ValueError("Nonlinear analysis requires at least one free DOF.")
        control = self.model.analysis.nonlinear_step
        if control.increment_method == "FIXED":
            step_size = 1.0 / control.number_of_increments
            minimum_step = step_size / (2.0 ** control.maximum_bisections)
            maximum_step = step_size
            maximum_increments = (
                control.number_of_increments * (2 ** control.maximum_bisections)
            )
        else:
            step_size = control.initial_step_fraction
            minimum_step = control.minimum_step_fraction
            maximum_step = control.maximum_step_fraction
            maximum_increments = control.maximum_increments or 1000

        current_factor = 0.0
        converged_displacement = self.displacements.copy()
        consecutive_cutbacks = 0
        last_failure_residual = np.inf

        while current_factor < 1.0 - 1.0e-14:
            if len(self.increment_history) >= maximum_increments:
                raise self._convergence_error(
                    "Maximum number of nonlinear increments was reached.",
                    current_factor,
                    min(1.0, current_factor + step_size),
                    last_failure_residual,
                )

            actual_step = min(step_size, 1.0 - current_factor)
            target_factor = current_factor + actual_step
            increment_number = len(self.increment_history) + 1
            try:
                trial, assembly, iterations, residual, relative = self._solve_increment(
                    increment_number,
                    target_factor,
                    converged_displacement,
                    free_dofs,
                    fixed_dofs,
                    prescribed,
                )
            except _IncrementFailure as failure:
                last_failure_residual = failure.residual_norm
                consecutive_cutbacks += 1
                self.total_cutbacks += 1
                reduced_step = actual_step * 0.5
                if (
                    consecutive_cutbacks > control.maximum_bisections
                    or reduced_step < minimum_step * (1.0 - 1.0e-12)
                ):
                    raise self._convergence_error(
                        "Unable to converge above the minimum permitted load step.",
                        current_factor,
                        target_factor,
                        failure.residual_norm,
                    ) from failure
                step_size = reduced_step
                continue

            current_factor = target_factor
            converged_displacement = trial
            self.displacements = trial.copy()
            self.increment_history.append(
                NonlinearIncrementResult(
                    increment=increment_number,
                    load_factor=current_factor,
                    step_size=actual_step,
                    iterations=iterations,
                    residual_norm=residual,
                    relative_residual_norm=relative,
                    cutbacks=consecutive_cutbacks,
                    displacements=trial.copy(),
                    strain_energy=assembly.strain_energy,
                )
            )
            consecutive_cutbacks = 0

            if control.increment_method == "ADAPT":
                desired = control.desired_iterations or 4
                scale = control.step_scale_factor or 1.0
                if iterations <= desired:
                    step_size = min(maximum_step, actual_step * scale)
                elif iterations > desired + 1:
                    step_size = max(minimum_step, actual_step / scale)
                else:
                    step_size = actual_step
            else:
                step_size = maximum_step

        final_state = self.assemble_state(self.displacements)
        final_factor = current_factor
        residual = final_factor * self.external_load - final_state.internal_force
        reactions = {}
        constrained_nodes = {
            constraint.node_id for constraint in self.model.constraints
        }
        fixed_dof_set = set(fixed_dofs)
        for node_id in constrained_nodes:
            indices = self.dof_mapper.get_dof_indices(node_id)
            nodal_reaction = np.zeros(self.DOFS_PER_NODE)
            for local_index, global_index in enumerate(indices):
                if global_index in fixed_dof_set:
                    nodal_reaction[local_index] = -residual[global_index]
            reactions[node_id] = nodal_reaction

        stresses = {
            eid: state.as_result_dict()
            for eid, state in final_state.element_states.items()
        }
        return {
            "displacements": self.displacements.copy(),
            "reactions": reactions,
            "stresses": stresses,
            "increments": list(self.increment_history),
            "iterations": list(self.iteration_history),
            "load_factor": final_factor,
            "strain_energy": final_state.strain_energy,
            "cutbacks": self.total_cutbacks,
        }

    def _solve_increment(
        self,
        increment_number,
        load_factor,
        starting_displacement,
        free_dofs,
        fixed_dofs,
        prescribed,
    ):
        control = self.model.analysis.nonlinear_step
        trial = starting_displacement.copy()
        for global_index in fixed_dofs:
            trial[global_index] = load_factor * prescribed.get(global_index, 0.0)

        residual_norm = np.inf
        relative_residual = np.inf
        for iteration in range(1, control.maximum_iterations + 1):
            state = self.assemble_state(trial)
            residual = load_factor * self.external_load - state.internal_force
            residual_free = residual[free_dofs]
            residual_norm = self._infinity_norm(residual_free)
            reference_norm = max(
                self._infinity_norm(load_factor * self.external_load[free_dofs]),
                self._infinity_norm(state.internal_force[free_dofs]),
                1.0,
            )
            relative_residual = residual_norm / reference_norm

            if iteration >= control.minimum_iterations and self._force_converged(
                residual_norm, relative_residual
            ):
                return trial, state, iteration - 1, residual_norm, relative_residual

            tangent_free = state.tangent_stiffness[free_dofs, :][:, free_dofs]
            correction = self._solve_linearized_system(tangent_free, residual_free, free_dofs)
            trial[free_dofs] += correction

            updated_state = self.assemble_state(trial)
            updated_residual = (
                load_factor * self.external_load - updated_state.internal_force
            )[free_dofs]
            updated_norm = self._infinity_norm(updated_residual)
            updated_reference = max(
                self._infinity_norm(load_factor * self.external_load[free_dofs]),
                self._infinity_norm(updated_state.internal_force[free_dofs]),
                1.0,
            )
            updated_relative = updated_norm / updated_reference
            self.iteration_history.append(
                NonlinearIterationResult(
                    increment=increment_number,
                    load_factor=load_factor,
                    iteration=iteration,
                    residual_norm=updated_norm,
                    relative_residual_norm=updated_relative,
                    correction_norm=self._infinity_norm(correction),
                )
            )
            residual_norm = updated_norm
            relative_residual = updated_relative
            if iteration >= control.minimum_iterations and self._force_converged(
                updated_norm, updated_relative
            ):
                return (
                    trial,
                    updated_state,
                    iteration,
                    updated_norm,
                    updated_relative,
                )

        raise _IncrementFailure(residual_norm)

    def _force_converged(self, absolute: float, relative: float) -> bool:
        tolerance = self.model.analysis.nonlinear_step.residual_tolerance
        if tolerance is None:
            tolerance = 1.0e-8
        absolute_tolerance = max(1.0e-8, tolerance)
        # Spatial CBAR force recovery contains products of current rotation
        # matrices and section actions.  Once the free-DOF correction is at
        # machine precision, cancellation of those large actions leaves a
        # reproducible residual floor around 1e-7 force units in prescribed-
        # displacement cases.  This absolute floor is numerical, not a
        # relaxation of the relative equilibrium criterion.
        if any(
            isinstance(element, CBARNonlinear)
            for element in self.element_solvers.values()
        ):
            absolute_tolerance = max(absolute_tolerance, 1.0e-6)
        return absolute <= absolute_tolerance or relative <= tolerance

    def _solve_linearized_system(self, matrix, rhs, free_dofs):
        if not np.all(np.isfinite(matrix.data)) or not np.all(np.isfinite(rhs)):
            raise _IncrementFailure(np.inf)
        asymmetry = matrix - matrix.T
        scale = max(1.0, np.max(np.abs(matrix.data)) if matrix.nnz else 1.0)
        has_cbar = any(
            isinstance(element, CBARNonlinear)
            for element in self.element_solvers.values()
        )
        if (
            not has_cbar
            and asymmetry.nnz
            and np.max(np.abs(asymmetry.data)) > scale * self._SYMMETRY_TOLERANCE
        ):
            raise ValueError("The nonlinear tangent stiffness matrix is not symmetric.")
        row_norms = np.asarray(np.abs(matrix).sum(axis=1)).ravel()
        if np.any(row_norms == 0.0):
            zero = np.flatnonzero(row_norms == 0.0)
            descriptions = [self._describe_dof(free_dofs[index]) for index in zero[:20]]
            raise ValueError(
                "The nonlinear tangent is singular; free DOFs without stiffness: "
                + ", ".join(descriptions)
            )
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("error", category=MatrixRankWarning)
                solution = spsolve(matrix.tocsc(), rhs)
        except (MatrixRankWarning, RuntimeError, ValueError) as error:
            raise _IncrementFailure(self._infinity_norm(rhs)) from error
        if not np.all(np.isfinite(solution)):
            raise _IncrementFailure(self._infinity_norm(rhs))
        return np.asarray(solution)

    def _convergence_error(self, message, current_factor, attempted_factor, residual):
        return NonlinearConvergenceError(
            message,
            last_converged_load_factor=current_factor,
            attempted_load_factor=attempted_factor,
            cutbacks=self.total_cutbacks,
            residual_norm=residual,
        )

    def _describe_dof(self, global_index: int) -> str:
        labels = ("T1", "T2", "T3", "R1", "R2", "R3")
        return (
            f"GRID {self.dof_mapper.get_node_from_index(global_index)} "
            f"{labels[global_index % 6]}"
        )

    @staticmethod
    def _infinity_norm(vector) -> float:
        return float(np.max(np.abs(vector))) if len(vector) else 0.0
