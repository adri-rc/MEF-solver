"""Format-independent data structures for the finite-element model."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional

import numpy as np

@dataclass
class Node:
    """A GRID entry with both its input and basic-coordinate positions.

    ``x1``--``x3`` are expressed in ``coord_id`` (GRID.CP), while ``x``--``z``
    are the corresponding Cartesian coordinates in the basic frame. The
    displacement frame is identified by ``cd_id`` (GRID.CD).
    """

    nid: int
    coord_id: int
    x1: float
    x2: float
    x3: float
    cd_id: int
    ps: str
    seid: int
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    # Maps vectors expressed in the basic/global frame to the GRID displacement
    # frame (CD). It is the identity for CD=0. The solver uses the same frame
    # for nodal translations and rotations, exactly as Nastran does.
    dof_transform: np.ndarray = field(
        default_factory=lambda: np.eye(3), repr=False, compare=False
    )

    def __post_init__(self):
        if self.nid <= 0:
            raise ValueError("Node ID must be a positive integer.")
        raw = np.asarray([self.x1, self.x2, self.x3], dtype=float)
        if not np.all(np.isfinite(raw)):
            raise ValueError(f"Node {self.nid} contains non-finite coordinates.")

    def get_coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])
    

@dataclass
class Material:
    """Base class for material definitions."""

    mid: int


@dataclass
class IsotropicMaterial(Material):
    """Linear isotropic material equivalent to Nastran MAT1.

    Any two of Young's modulus ``E``, shear modulus ``G`` and Poisson's ratio
    ``nu`` determine the third. When all three are supplied, they must agree
    within the project's five-percent engineering tolerance.
    """

    E: Optional[float] = None
    nu: Optional[float] = None
    G: Optional[float] = None
    rho: Optional[float] = None

    def __post_init__(self):
        supplied = [value for value in (self.E, self.G, self.nu) if value is not None]
        if len(supplied) < 2:
            raise ValueError(
                f"Review inputs in material '{self.mid}'\n"
                "There must be at least 2 parameters defined between (E, nu, G)."
            )
        if not all(np.isfinite(value) for value in supplied):
            raise ValueError(
                f"Material '{self.mid}' contains non-finite elastic constants."
            )

        if self.E is None:
            self.E = 2 * self.G * (1 + self.nu)
        elif self.G is None:
            if self.nu == -1:
                raise ValueError(
                    f"Review inputs in material '{self.mid}'\n"
                    "Poisson's ratio cannot be -1."
                )
            self.G = self.E / (2 * (1 + self.nu))
        elif self.nu is None:
            if self.G == 0:
                raise ValueError(
                    f"Review inputs in material '{self.mid}'\n"
                    "G cannot be 0 to calculate nu."
                )
            self.nu = (self.E / (2 * self.G)) - 1
        else:
            E_calculated = 2 * self.G * (1 + self.nu)
            if self.E == 0.0:
                raise ValueError(
                    f"Review inputs in material '{self.mid}'\n"
                    "E must be greater than zero."
                )
            error = abs(self.E - E_calculated) / abs(self.E)
            if error > 0.05:
                raise ValueError(
                    f"Material ID '{self.mid}' is physically inconsistent. \n"
                    f"        Calculated E based on G/nu: {E_calculated:.2e}\n"
                    f"        User input E:               {self.E:.2e}\n"
                    f"        Difference is > 5% ({error:.1%}). Review material inputs."
                )

        if self.E <= 0.0 or self.G <= 0.0:
            raise ValueError(f"Material '{self.mid}' requires E > 0 and G > 0.")
        if not (-1.0 < self.nu <= 0.5):
            raise ValueError(f"Material '{self.mid}' requires -1 < nu <= 0.5.")
        if self.rho is not None and (
            not np.isfinite(self.rho) or self.rho < 0.0
        ):
            raise ValueError(
                f"Material '{self.mid}' requires a finite, non-negative density."
            )

    def __repr__(self):
        return (
            f"IsotropicMaterial with id:{self.mid}\n"
            f"(E={self.E:.2e}, nu={self.nu:.4f}, "
            f"G={self.G:.2e}, rho={self.rho})"
        )


@dataclass
class Property:
    """Base class for element properties."""

    pid: int

@dataclass
class RodProperty(Property):
    """PROD section properties for a CROD element."""

    mat_id: int
    A: float
    J: float
    C: float
    NSM: float

    def __post_init__(self):
        values = (self.A, self.J, self.C, self.NSM)
        if not all(np.isfinite(value) for value in values):
            raise ValueError(f"RodProperty {self.pid} contains non-finite values.")
        if self.A <= 0.0:
            raise ValueError(f"RodProperty {self.pid} requires A > 0.")
        if self.J < 0.0 or self.NSM < 0.0:
            raise ValueError(f"RodProperty {self.pid} requires J and NSM to be non-negative.")


@dataclass
class BarProperty(Property):
    """PBAR section properties for a CBAR element."""

    mat_id: int
    A: float
    I1: float
    I2: float
    J: float
    I12: float = 0.0
    K1: float = 0.0
    K2: float = 0.0
    NSM: float = 0.0

    # PBAR recovery points C--F, stored as (y, z) pairs in section axes.
    stress_points: List[tuple] = field(default_factory=list)
    # Optional SOL 400 PBARN1 extension associated with this PBAR. Keeping
    # these fields on the format-independent property lets the nonlinear
    # solver distinguish advanced beam behavior without coupling it to
    # pyNastran (which rejects PBARN1 in the supported version).
    nonlinear_section: Optional[str] = None
    nonlinear_behavior: Optional[str] = None
    nonlinear_integration: Optional[str] = None

    def __post_init__(self):
        for attribute in (
            "nonlinear_section",
            "nonlinear_behavior",
            "nonlinear_integration",
        ):
            value = getattr(self, attribute)
            if value is not None:
                setattr(self, attribute, str(value).upper())
        values = (
            self.A,
            self.I1,
            self.I2,
            self.J,
            self.I12,
            self.K1,
            self.K2,
            self.NSM,
        )
        if not all(np.isfinite(value) for value in values):
            raise ValueError(f"BarProperty {self.pid} contains non-finite values.")
        if self.A <= 0.0 or self.I1 <= 0.0 or self.I2 <= 0.0:
            raise ValueError(
                f"BarProperty {self.pid} requires A, I1 and I2 to be positive."
            )
        if self.J < 0.0 or self.K1 < 0.0 or self.K2 < 0.0 or self.NSM < 0.0:
            raise ValueError(
                f"BarProperty {self.pid} requires J, K1, K2 and NSM to be "
                "non-negative."
            )
        # A positive-definite centroidal inertia tensor requires a strict inequality.
        if (self.I1 * self.I2) <= (self.I12**2):
            raise ValueError(
                f"Inconsistent inertial properties in BarProperty with id: {self.pid}"
            )


@dataclass
class ShellProperty(Property):
    """PSHELL thickness, material assignments and stiffness factors."""

    mat_id: Optional[int]
    t: float
    mat_id_bending: Optional[int] = None
    mat_id_shear: Optional[int] = None
    bending_inertia_ratio: float = 1.0
    shear_thickness_ratio: float = 5.0 / 6.0

    def __post_init__(self):
        if not np.isfinite(self.t) or self.t <= 0.0:
            raise ValueError(f"ShellProperty {self.pid} requires a finite thickness t > 0.")
        if self.mat_id == -1:
            self.mat_id = None
        if self.mat_id_bending == -1:
            self.mat_id_bending = None
        if self.mat_id_shear == -1:
            self.mat_id_shear = None
        if self.mat_id_shear is not None and self.mat_id_bending is None:
            raise ValueError(
                f"ShellProperty {self.pid} cannot define MID3 without MID2."
            )
        if self.mat_id is None and self.mat_id_bending is None:
            raise ValueError(
                f"ShellProperty {self.pid} requires MID1 or MID2."
            )
        if (
            not np.isfinite(self.bending_inertia_ratio)
            or self.bending_inertia_ratio <= 0.0
        ):
            raise ValueError(f"ShellProperty {self.pid} requires 12I/T^3 > 0.")
        if (
            not np.isfinite(self.shear_thickness_ratio)
            or self.shear_thickness_ratio <= 0.0
        ):
            raise ValueError(f"ShellProperty {self.pid} requires TS/T > 0.")


@dataclass
class Element:
    """Base class for elements."""
    eid: int
    prop_id: int
    nodes: List[int]

@dataclass
class RodElement(Element):
    """Two-node CROD element with axial and torsional response."""

    pass

@dataclass
class BarElement(Element):
    """Two-node CBAR element, including orientation, offsets and releases."""

    # When both are present, the solver uses the G0 node for orientation.
    orientation_vector: Optional[List[float]] = None
    g0: Optional[int] = None

    offt: str = "GGG"
    wa: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    wb: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])

    # Integers encode released end DOFs; zero means full continuity.
    pa: int = 0
    pb: int = 0

@dataclass
class CQUAD4Element(Element):
    """Four-node CQUAD4 with a material-axis angle in degrees."""

    theta: float = 0.0

@dataclass
class CTRIA3Element(Element):
    """Three-node CTRIA3 with a material-axis angle in degrees."""

    theta: float = 0.0


@dataclass
class LoadCase:
    """Base class for loads grouped by set ID."""

    sid: int


@dataclass
class SinglePointLoad(LoadCase):
    """Concentrated FORCE with magnitude and direction components."""

    appnode_id: int
    module: float
    vector: List[float]
    cid: int = 0

    def __post_init__(self):
        vector = np.asarray(self.vector, dtype=float)
        if vector.shape != (3,) or not np.all(np.isfinite(vector)):
            raise ValueError("FORCE direction must contain three finite components.")
        if not np.isfinite(self.module):
            raise ValueError("FORCE magnitude must be finite.")


@dataclass
class SinglePointMoment(LoadCase):
    """Concentrated MOMENT with magnitude and direction components."""

    appnode_id: int
    module: float
    vector: List[float]
    cid: int = 0

    def __post_init__(self):
        vector = np.asarray(self.vector, dtype=float)
        if vector.shape != (3,) or not np.all(np.isfinite(vector)):
            raise ValueError("MOMENT direction must contain three finite components.")
        if not np.isfinite(self.module):
            raise ValueError("MOMENT magnitude must be finite.")


@dataclass
class Constraint:
    """Base class for boundary conditions grouped by set ID."""

    sid: int


@dataclass
class SPCConstraint(Constraint):
    """SPC or SPC1 constraint on one or more nodal components."""

    node_id: int
    components: str
    value: float = 0.0

    def __post_init__(self):
        components = str(self.components)
        if not components or any(component not in "123456" for component in components):
            raise ValueError(f"Invalid SPC components '{self.components}'.")
        if len(set(components)) != len(components):
            raise ValueError(f"Repeated SPC components in '{self.components}'.")
        if not np.isfinite(self.value):
            raise ValueError("SPC enforced displacement must be finite.")


class AnalysisType(str, Enum):
    """Analysis families understood by the project orchestration layer."""

    LINEAR_STATIC = "linear_static"
    NONLINEAR_STATIC = "nonlinear_static"


@dataclass
class NonlinearStepControl:
    """Format-independent subset of Nastran NLSTEP/NLPARM controls.

    Values that Nastran derives from ``CTRLDEF`` are left as ``None`` unless
    explicitly supplied, so the model does not silently embed defaults tied
    to one MSC Nastran version.
    """

    sid: int
    source: str = "NLSTEP"
    total_time: float = 1.0
    control_definition: Optional[str] = None
    increment_method: str = "FIXED"
    number_of_increments: Optional[int] = 50
    initial_step_fraction: Optional[float] = None
    minimum_step_fraction: Optional[float] = None
    maximum_step_fraction: Optional[float] = None
    desired_iterations: Optional[int] = None
    step_scale_factor: Optional[float] = None
    output_intervals: int = 0
    maximum_increments: Optional[int] = None
    maximum_iterations: int = 10
    minimum_iterations: int = 1
    maximum_bisections: int = 10
    convergence: str = "PV"
    displacement_tolerance: Optional[float] = None
    residual_tolerance: Optional[float] = None
    work_tolerance: Optional[float] = None
    stiffness_update_method: str = "PFNT"
    stiffness_update_interval: Optional[int] = None

    def __post_init__(self):
        self.source = self.source.upper()
        self.increment_method = self.increment_method.upper()
        self.convergence = self.convergence.upper()
        self.stiffness_update_method = self.stiffness_update_method.upper()
        if self.control_definition is not None:
            self.control_definition = self.control_definition.upper()

        if self.sid <= 0:
            raise ValueError("Nonlinear control ID must be a positive integer.")
        if self.source not in {"NLSTEP", "NLPARM"}:
            raise ValueError(f"Unsupported nonlinear control source: {self.source}")
        if not np.isfinite(self.total_time) or self.total_time <= 0.0:
            raise ValueError("Nonlinear step total_time must be finite and positive.")
        if self.control_definition not in {None, "QLINEAR", "MILDLY", "SEVERELY"}:
            raise ValueError(
                f"Unsupported NLSTEP control definition: {self.control_definition}"
            )
        if self.increment_method not in {"FIXED", "ADAPT"}:
            raise ValueError(
                f"Unsupported nonlinear increment method: {self.increment_method}"
            )
        if self.increment_method == "FIXED":
            if self.number_of_increments is None or self.number_of_increments <= 0:
                raise ValueError("FIXED nonlinear stepping requires NINC > 0.")
        else:
            fractions = (
                self.initial_step_fraction,
                self.minimum_step_fraction,
                self.maximum_step_fraction,
            )
            if any(
                value is None or not np.isfinite(value) or value <= 0.0
                for value in fractions
            ):
                raise ValueError(
                    "ADAPT stepping requires positive DTINITF, DTMINF and DTMAXF."
                )
            if not (
                self.minimum_step_fraction
                <= self.initial_step_fraction
                <= self.maximum_step_fraction
            ):
                raise ValueError(
                    "ADAPT step fractions must satisfy "
                    "DTMINF <= DTINITF <= DTMAXF."
                )
            if self.desired_iterations is not None and self.desired_iterations <= 0:
                raise ValueError("ADAPT desired_iterations must be positive.")
            if self.step_scale_factor is not None and (
                not np.isfinite(self.step_scale_factor) or self.step_scale_factor <= 0.0
            ):
                raise ValueError("ADAPT step_scale_factor must be finite and positive.")
            if self.maximum_increments is not None and self.maximum_increments <= 0:
                raise ValueError("ADAPT maximum_increments must be positive.")
        if self.output_intervals < -1:
            raise ValueError(
                "Nonlinear output_intervals must be -1, 0 or a positive integer."
            )
        if self.minimum_iterations < 0 or self.maximum_iterations <= 0:
            raise ValueError("Nonlinear iteration limits must be non-negative/positive.")
        if self.minimum_iterations > self.maximum_iterations:
            raise ValueError("MINITER cannot exceed MAXITER.")
        if self.maximum_bisections < 0:
            raise ValueError("MAXBIS must be non-negative in this project.")
        if not self.convergence or not set(self.convergence) <= set("UPWVNA"):
            raise ValueError(f"Unsupported nonlinear convergence flags: {self.convergence}")
        for name, value in (
            ("displacement_tolerance", self.displacement_tolerance),
            ("residual_tolerance", self.residual_tolerance),
            ("work_tolerance", self.work_tolerance),
        ):
            if value is not None and (not np.isfinite(value) or value <= 0.0):
                raise ValueError(f"{name} must be finite and positive when supplied.")
        if self.stiffness_update_method not in {"PFNT", "ITER", "AUTO", "SEMI"}:
            raise ValueError(
                f"Unsupported stiffness update method: {self.stiffness_update_method}"
            )
        if (
            self.stiffness_update_method == "ITER"
            and (
                self.stiffness_update_interval is None
                or self.stiffness_update_interval <= 0
            )
        ):
            raise ValueError("ITER stiffness updates require KSTEP > 0.")


@dataclass
class AnalysisConfiguration:
    """Analysis selection and global geometric/load-stiffness switches."""

    solution_sequence: Optional[int] = None
    analysis_type: AnalysisType = AnalysisType.LINEAR_STATIC
    analysis_keyword: Optional[str] = None
    large_displacement_parameter: Optional[int] = None
    large_displacements: bool = False
    follower_force_stiffness: bool = False
    nonlinear_step: Optional[NonlinearStepControl] = None


class FEMModel:
    """Container for the normalized physical model consumed by the solvers."""

    def __init__(self):
        self.nodes: Dict[int, Node] = {}
        self.materials: Dict[int, Material] = {}
        self.properties: Dict[int, Property] = {}
        self.elements: Dict[int, Element] = {}

        # Lists preserve repeated SIDs after LOAD/SPCADD expansion.
        self.constraints: List[Constraint] = []
        self.loads: List[LoadCase] = []

        self.k6rot: float = 0.0
        self.active_subcase_id: Optional[int] = None
        self.active_load_sid: Optional[int] = None
        self.active_spc_sid: Optional[int] = None
        # Logical Bulk Data entry counts retained for an F06-like model summary.
        # They are deliberately separate from normalized loads and constraints,
        # which may expand or combine the source entries used by the solver.
        self.source_entry_counts: Dict[str, int] = {}
        self.analysis = AnalysisConfiguration()

    def validate(self, require_elements: bool = True) -> None:
        """Validate references and finite geometry before matrix assembly."""
        if require_elements and not self.nodes:
            raise ValueError("The FEM model does not contain any GRID nodes.")
        if require_elements and not self.elements:
            raise ValueError("The FEM model does not contain supported finite elements.")
        for node in self.nodes.values():
            coords = node.get_coords()
            if not np.all(np.isfinite(coords)):
                raise ValueError(f"Node {node.nid} contains non-finite global coordinates.")
            transform = np.asarray(node.dof_transform, dtype=float)
            if transform.shape != (3, 3) or not np.all(np.isfinite(transform)):
                raise ValueError(f"Node {node.nid} has an invalid CD transformation.")
            if not np.allclose(transform @ transform.T, np.eye(3), atol=1e-10):
                raise ValueError(f"Node {node.nid} CD transformation is not orthonormal.")
        for element in self.elements.values():
            if element.prop_id not in self.properties:
                raise ValueError(f"Element {element.eid} references missing property {element.prop_id}.")
            missing_nodes = [nid for nid in element.nodes if nid not in self.nodes]
            if missing_nodes:
                raise ValueError(f"Element {element.eid} references missing nodes {missing_nodes}.")
            prop = self.properties[element.prop_id]
            for attribute in ("mat_id", "mat_id_bending", "mat_id_shear"):
                mid = getattr(prop, attribute, None)
                if mid is not None and mid not in self.materials:
                    raise ValueError(f"Property {prop.pid} references missing material {mid}.")
        for load in self.loads:
            if load.appnode_id not in self.nodes:
                raise ValueError(f"Load SID {load.sid} references missing node {load.appnode_id}.")
        for constraint in self.constraints:
            if constraint.node_id not in self.nodes:
                raise ValueError(
                    f"Constraint SID {constraint.sid} references missing node {constraint.node_id}."
                )

    def get_node(self, node_id: int) -> Node:
        """Return a node by ID, raising ``KeyError`` when it is absent."""
        if node_id not in self.nodes:
            raise KeyError(
                f"[FEMModel ERROR] Node ID: '{node_id}' does not exist in the model."
            )
        return self.nodes[node_id]

