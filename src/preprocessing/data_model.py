from dataclasses import dataclass, field
from enum import Enum
from typing import List, Dict, Optional
import numpy as np

"""
@file data_model.py
@brief Definition of data structures for the FEM model.
@details This module contains the classes that represent the physical 
    entities of the problem (nodes, elements, materials) 
    DECOUPLED from the input format.
"""

@dataclass
class Node:
    """
    @class Node
    @brief Generic node definition where raw data from input file can be saved and general 
        coordinates in global sistem are also saved.
    
    @details
        This class stores both the “raw” data from the input file (such as a Nastran GRID card)
        and the coordinates processed in the global Cartesian system,
        which are necessary for matrix assembly.

    @var nid ID único del nodo (Entero > 0).
    @var coord_id ID del sistema de coordenadas (CP) donde se definen x1, x2, x3. (0 = Global Cartesiano).
    @var x1 Primera coordenada en el sistema local CP.
    @var x2 Segunda coordenada en el sistema local CP.
    @var x3 Tercera coordenada en el sistema local CP.
    @var cd_id ID del sistema de coordenadas de salida (CD) para los desplazamientos.
    @var ps Restricción permanente de punto único (Permanent SPC), ej: "123".
    @var seid ID del super-elemento (SEID).
    @var x Coordenada X procesada en el Sistema Cartesiano Global (Calculada).
    @var y Coordenada Y procesada en el Sistema Cartesiano Global (Calculada).
    @var z Coordenada Z procesada en el Sistema Cartesiano Global (Calculada).
    """
    nid: int # node id 
    # Raw data fom input file
    coord_id: int # Identification number of coordinate system where node is defined.
    x1: float
    x2: float
    x3: float
    cd_id: int
    ps: str # Permanent Single Point Constraint
    seid: int # Superelement ID
    
    # Processed coordinates (Global Cartesian) - Calculated after reading input file
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    # Maps vectors expressed in the basic/global frame to the GRID displacement
    # frame (CD).  It is the identity for CD=0.  The solver uses the same frame
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
    

# region Materials classes
@dataclass
class Material:
    """
    @class Material
    @brief Base class for properties constituting materials.
    """
    mid: int #: ID material

# region IsotropicMaterial(Material)
@dataclass
class IsotropicMaterial(Material):
    """
    @class IsotropicMaterial
    @brief Linear isotropic elastic material (equivalent to NASTRAN MAT1).
    @details 
    This class implements the elastic constant triangulation logic.
    The triplet (E, G, nu) is guaranteed to be complete and consistent.    
    Allows the class to be instantiated with ANY pair of constants:
        1. (E, nu) -> Calculates G
        2. (E, G) -> Calculates nu
        3. (G, nu) -> Calculates E
    @var E Módulo de Young's Module [Pa].
    @var nu Coeficiente de Poisson [-].
    @var G Shear Module [Pa].
    @var rho density [kg/m^3].

    """
    E: Optional[float] = None   #: Young´s Module
    nu: Optional[float] = None  #: Poisson Coef
    G: Optional[float] = None   #: Shear module (optional, calculado si es 0)
    rho: Optional[float] = None #: density (if needed in future)

    def __post_init__(self):
            """
            Triangulation logic to calculate the missing constant.
            """
            # How many elastic parameters have been defined is count (excluding rho).
            params = [self.E is not None, self.nu is not None, self.G is not None]
            inputs_count = sum(params)

            # Test 1: There must be at least 2 parameters.
            if inputs_count < 2:
                raise ValueError(f"Review inputs in material '{self.mid}'\nThere must be at least 2 parameters defined between (E, nu, G).")

            # --- CALC ---
                        
            supplied = [value for value in (self.E, self.G, self.nu) if value is not None]
            if not all(np.isfinite(value) for value in supplied):
                raise ValueError(f"Material '{self.mid}' contains non-finite elastic constants.")

            if self.E is None: # CASE 1: Missing E -> (G, nu) Known
                self.E = 2 * self.G * (1 + self.nu) # Formula: E = 2G(1 + nu)

            
            elif self.G is None: # CASE 2: Missing G -> (E, nu) Known
                if self.nu == -1: # Avoid theoretical division by zero
                    raise ValueError(f"Review inputs in material '{self.mid}'\nPoisson's ratio cannot be -1.")
                self.G = self.E / (2 * (1 + self.nu)) # Formula: G = E / (2(1 + nu))

            
            elif self.nu is None: # CASE 3: Missing nu -> (E, G) Known
                
                if self.G == 0:
                    raise ValueError(f"Review inputs in material '{self.mid}'\nG cannot be 0 to calculate nu.")
                self.nu = (self.E / (2 * self.G)) - 1 # Formula: nu = (E / 2G) - 1

            else: # CASE 4: The three parameters (E, G and nu) were given -> Consistency validation
                E_calculated = 2 * self.G * (1 + self.nu) # We verify whether the three values given comply with physics.
                if self.E == 0.0:
                    raise ValueError(f"Review inputs in material '{self.mid}'\nE must be greater than zero.")
                error = abs(self.E - E_calculated) / abs(self.E)
                
                
                if error > 0.05: # If the error is greater than 5%, it´s notified (engineering tolerance).
                    raise ValueError(
                        f"Material ID '{self.mid}' is physically inconsistent. \n"
                        f"        Calculated E based on G/nu: {E_calculated:.2e}\n"
                        f"        User input E:               {self.E:.2e}\n"
                        f"        Difference is > 5% ({error:.1%}). Review material inputs."
                    )

            if self.E <= 0.0 or self.G <= 0.0:
                raise ValueError(f"Material '{self.mid}' requires E > 0 and G > 0.")
            if not (-1.0 < self.nu <= 0.5):
                raise ValueError(
                    f"Material '{self.mid}' requires -1 < nu <= 0.5."
                )
            if self.rho is not None and (not np.isfinite(self.rho) or self.rho < 0.0):
                raise ValueError(f"Material '{self.mid}' requires a finite, non-negative density.")
                    

    def __repr__(self):
        # Usamos 'or 0' por seguridad en el f-string, aunque post_init asegura valores
        return (f"IsotropicMaterial with id:{self.mid}\n"
                f"(E={self.E:.2e}, nu={self.nu:.4f}, G={self.G:.2e}, rho={self.rho})")
# endregion IsotropicMaterial(Material)

# endregion Material Classes

# region Property Classes
@dataclass
class Property:
    """
    @class Property
    @brief Base class for any elemen property.
    """
    pid: int #: ID property
    

# region RodProperty(Property)
@dataclass
class RodProperty(Property):
    """
    @class RodProperty
    @brief Corresponds to PROD in NASTRAN.
    """
    mat_id: int  #: ID from asociated material.
    A: float     #: Transverse section area.
    J: float     #: Torsion constant (Polar inertia for circulars coordinates)
    C: float     #: Coefficient for calculating torsional stresses
    NSM: float   #: Non-structural mass per unit length

    def __post_init__(self):
        values = (self.A, self.J, self.C, self.NSM)
        if not all(np.isfinite(value) for value in values):
            raise ValueError(f"RodProperty {self.pid} contains non-finite values.")
        if self.A <= 0.0:
            raise ValueError(f"RodProperty {self.pid} requires A > 0.")
        if self.J < 0.0 or self.NSM < 0.0:
            raise ValueError(f"RodProperty {self.pid} requires J and NSM to be non-negative.")
# endregion RodProperty(Property)

# region BarProperty(Property)
@dataclass
class BarProperty(Property):
    '''
    @class BarProperty
    @brief Geometric properties for the beam element (PBAR).
    @details Defines the area, inertia, and torsional constant.
    '''
    mat_id: int         #: ID from asociated material.
    A: float            #: Transverse section area.
    I1: float           #: Moment of inertia of area around plane 1 (local z) - XY bending
    I2: float           #: Moment of inertia of area around plane 2 (and local) - Bending XZ
    J: float            #: Torsion constant
    I12: float = 0.0    #: Producto de inercia (Iyz). Default 0 (sección simétrica)
    K1: float = 0.0     #: Area factor for shear, plane 1. Default 0 (Euler-Bernoulli)
    K2: float = 0.0     #: Area factor for shear, plane 2. Default 0 (Euler-Bernoulli)
    NSM: float = 0.0    #: Non-structural mass per unit length
    
    # Stress Recovery Points
    # List of tuples [(y_c, z_c), (y_d, z_d), (y_e, z_e), (y_f, z_f)]
    # If the list is empty, it will be calculated only at the centroid (0,0) or default endpoints.
    stress_points: List[tuple] = field(default_factory=list)
    # Optional SOL 400 PBARN1 extension associated with this PBAR.  Keeping
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
        values = (self.A, self.I1, self.I2, self.J, self.I12, self.K1, self.K2, self.NSM)
        if not all(np.isfinite(value) for value in values):
            raise ValueError(f"BarProperty {self.pid} contains non-finite values.")
        if self.A <= 0.0 or self.I1 <= 0.0 or self.I2 <= 0.0:
            raise ValueError(f"BarProperty {self.pid} requires A, I1 and I2 to be positive.")
        if self.J < 0.0 or self.K1 < 0.0 or self.K2 < 0.0 or self.NSM < 0.0:
            raise ValueError(f"BarProperty {self.pid} requires J, K1, K2 and NSM to be non-negative.")
        # A positive-definite centroidal inertia tensor requires a strict inequality.
        if (self.I1 * self.I2) <= (self.I12 ** 2):
            raise ValueError(f"Inconsistent inertial properties in BarProperty with id: {self.pid}")
# endregion BarProperty(Property)

# region ShellProperty(Property)
@dataclass
class ShellProperty(Property):
    """
    @class ShellProperty
    @brief Properties for Shell elements (PSHELL).
    @details Defines thickness and material IDs for membrane and bending.
    """
    mat_id: Optional[int] #: Material ID for membrane (blank MID1 disables membrane)
    t: float             #: Thickness
    mat_id_bending: Optional[int] = None #: Material ID for bending (blank MID2 disables bending)
    mat_id_shear: Optional[int] = None #: Material ID for transverse-shear flexibility (MID3)
    bending_inertia_ratio: float = 1.0 #: PSHELL 12I/T^3 factor
    shear_thickness_ratio: float = 5.0 / 6.0 #: PSHELL TS/T factor

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
            raise ValueError(
                f"ShellProperty {self.pid} requires 12I/T^3 > 0."
            )
        if (
            not np.isfinite(self.shear_thickness_ratio)
            or self.shear_thickness_ratio <= 0.0
        ):
            raise ValueError(f"ShellProperty {self.pid} requires TS/T > 0.")
# endregion ShellProperty(Property)
# endregion Property Classes


# region Element Classes
@dataclass
class Element:
    """Base class for elements."""
    eid: int
    prop_id: int
    nodes: List[int]

@dataclass
class RodElement(Element):
    """
    @brief 2-node bar (CROD).
    @details Defines a tension-compression-torsion element"""
    pass

@dataclass
class BarElement(Element):
    """
    @class BarElement
    @brief 2-node prismatic beam element (CBAR).
    @details 
    Supports traction, torsion, and bending in two planes.
    Includes support for vector or node orientation, 
    offsets (eccentricities), and releases (pin flags).
    """
    # --- Orientation ---
    # Can be defined by vector (v) or by auxiliary node (g0).
    # Both are saved and solver will decide which one to use (priority to g0 if it exists).
    orientation_vector: Optional[List[float]] = None 
    g0: Optional[int] = None

    # --- Offsets
    offt: str = 'GGG' # offset type: 'GGG', 'B', 'LOC'. Default Global.
    wa: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0]) # Offset Node A
    wb: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0]) # Offset Node B

    # --- Releases (Joints/Pin Flags) ---
    # Integers representing the released DOFs (e.g., 456).
    pa: int = 0 # Default: 0 = Clamped (fully fixed), complete continuity between elements.
    pb: int = 0 # Default: 0 = Clamped (fully fixed), complete continuity between elements.
    



@dataclass
class CQUAD4Element(Element):
    """
    @brief 4-node Quadrilateral Plate Element (CQUAD4).
    @details Stores connectivity and material angle. The active solver uses
    the validated planar Nastran reconstruction.
    """
    theta: float = 0.0 # Material coordinate system orientation

@dataclass
class CTRIA3Element(Element):
    """
    @brief 3-node Triangular Shell Element (CTRIA3).
    @details Supports membrane and bending stiffness (Reissner-Mindlin).
    """
    theta: float = 0.0 # Material coordinate system orientation
# endregion Element Classes


# region LOADs Classes
@dataclass
class LoadCase:
    """
    @class LoadCase
    @brief Abstract base class for any load type (point, distributed, thermal).
    """
    sid: int  #: Set ID (SID) para agrupar casos de carga

# region SinglePointLoad(LoadCase) (FORCE)
@dataclass
class SinglePointLoad(LoadCase):
    """
    @class SinglePointLoad
    @brief Represents a static point load (FORCE card).
    @details Defines a force applied to a specific node, 
        with a magnitude and direction vector referenced to a coordinate system.
    """
    appnode_id: int      #: Application node ID (nid)
    module: float        #: Force module (F, scale factor)
    vector: List[float]  #: Director vector [N1, N2, N3]
    cid: int = 0    #: Coordinate system ID (CID) where the vector is defined. 0=Global.

    def __post_init__(self):
        vector = np.asarray(self.vector, dtype=float)
        if vector.shape != (3,) or not np.all(np.isfinite(vector)):
            raise ValueError("FORCE direction must contain three finite components.")
        if not np.isfinite(self.module):
            raise ValueError("FORCE magnitude must be finite.")
# endregion SinglePointLoad(LoadCase) (FORCE)

# region SinglePointMoment(LoadCase) (MOMENT)
@dataclass
class SinglePointMoment(LoadCase):
    """
    @class SinglePointMoment
    @brief Point-in-time applied to a node (MOMENT card).
    @details affects the degrees of rotational freedom. (For beams(CBAR) or plates (CQUAD4) elements.)
    """
    appnode_id: int      #: Application node ID (nid)
    module: float        #: Moment module (M, scale factor)
    vector: List[float]  #: Director vector [N1, N2, N3]
    cid: int = 0         #: Coordinate system ID (CID) where the vector is defined. 0=Global.

    def __post_init__(self):
        vector = np.asarray(self.vector, dtype=float)
        if vector.shape != (3,) or not np.all(np.isfinite(vector)):
            raise ValueError("MOMENT direction must contain three finite components.")
        if not np.isfinite(self.module):
            raise ValueError("MOMENT magnitude must be finite.")
# endregion SinglePointMoment(LoadCase) (MOMENT)

# endregion LOADS Classes

# region Constraints
@dataclass
class Constraint:
    """
    @class ConstraintS
    @brief Base class for boundary conditions and kinematic constraints.
    """
    sid: int  #: Set ID (SID) para agrupar restricciones (Load Case)

# region SPCConstraint
@dataclass
class SPCConstraint(Constraint):
    """
    @class SPCConstraint
    @brief Restricción de Punto Único (Tarjetas SPC y SPC1).
    @details Fija grados de libertad específicos de un nodo a un valor constante.
    """
    node_id: int        #: Grid ID (Gi)
    components: str     #: Componentes restringidas (ej. "123456" o "3")
    value: float = 0.0  #: Enforced displacement (Di). Default=0.0 (empotramiento).

    def __post_init__(self):
        components = str(self.components)
        if not components or any(component not in "123456" for component in components):
            raise ValueError(f"Invalid SPC components '{self.components}'.")
        if len(set(components)) != len(components):
            raise ValueError(f"Repeated SPC components in '{self.components}'.")
        if not np.isfinite(self.value):
            raise ValueError("SPC enforced displacement must be finite.")
# endregion SPCConstraint

# endregion Constraints


# region Analysis configuration
class AnalysisType(str, Enum):
    """Analysis families understood by the project orchestration layer."""

    LINEAR_STATIC = "linear_static"
    NONLINEAR_STATIC = "nonlinear_static"


@dataclass
class NonlinearStepControl:
    """Format-independent subset of Nastran NLSTEP/NLPARM controls.

    Values that Nastran derives from ``CTRLDEF`` are left as ``None`` unless
    they were explicitly supplied.  This prevents the future nonlinear solver
    from silently baking one MSC-version-dependent default into the model.
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
            if any(value is None or not np.isfinite(value) or value <= 0.0 for value in fractions):
                raise ValueError("ADAPT stepping requires positive DTINITF, DTMINF and DTMAXF.")
            if not (
                self.minimum_step_fraction
                <= self.initial_step_fraction
                <= self.maximum_step_fraction
            ):
                raise ValueError("ADAPT step fractions must satisfy DTMINF <= DTINITF <= DTMAXF.")
            if self.desired_iterations is not None and self.desired_iterations <= 0:
                raise ValueError("ADAPT desired_iterations must be positive.")
            if self.step_scale_factor is not None and (
                not np.isfinite(self.step_scale_factor) or self.step_scale_factor <= 0.0
            ):
                raise ValueError("ADAPT step_scale_factor must be finite and positive.")
            if self.maximum_increments is not None and self.maximum_increments <= 0:
                raise ValueError("ADAPT maximum_increments must be positive.")
        if self.output_intervals < -1:
            raise ValueError("Nonlinear output_intervals must be -1, 0 or a positive integer.")
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
            and (self.stiffness_update_interval is None or self.stiffness_update_interval <= 0)
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


# endregion Analysis configuration

# region FEMModel
class FEMModel:
    """
    @class FEMModel
    @brief Main container (“Database”) of the physical model.
    @details Stores all parsed and processed entities, 
        acting as the single source for the Solver.
    """
    def __init__(self):
        self.nodes: Dict[int, Node] = {}            #: Diccionario de Nodos {ID: Node}
        self.materials: Dict[int, Material] = {}    #: Diccionario de Materiales {ID: Material}
        self.properties: Dict[int, Property] = {}   #: Propiedades de sección {ID: Property}
        self.elements: Dict[int, Element] = {}      #: Elementos finitos {ID: Element}

                # Future user coordinate systems will be implemented
                #self.coords: Dict[int, CoordSystem] = {}    #: Sistemas de Coordenadas definidos por usuario
        
        # Loads and restrictions are stored in lists because different subcase can count with the same load_id or SPC_id multiple times,
        # therefor in a dictionary they can be overwritten.
        self.constraints: List[Constraint] = []     
        self.loads: List[LoadCase] = []

        # --- PARÁMETROS GLOBALES DE NASTRAN ---
        self.k6rot: float = 0.0                     #: Factor de rigidez de drilling (PARAM, K6ROT)
        self.active_subcase_id: Optional[int] = None
        self.active_load_sid: Optional[int] = None
        self.active_spc_sid: Optional[int] = None
        # Logical Bulk Data entry counts retained for an F06-like model summary.
        # They are deliberately separate from normalized loads/constraints,
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
        """
        @brief Retrieves a Node object by its ID.
        @param node_id ID of the node to search for.
        @return Class Node with the node id requested.
        @throws KeyError if the node does not exist in the model.
        """

        if node_id not in self.nodes:
            raise KeyError(f"[FEMModel ERROR] Node ID: '{node_id}' does not exist in the model.")
        
        return self.nodes[node_id]
    
        '''def get_coord_system(self, cid: int) -> Optional[CoordSystem]:
            """
            @brief Recupera un sistema de coordenadas.
            @return Objeto CoordSystem o None si es el global (CID=0).
            """
            if cid == 0:
                return None
            return self.coords[cid]'''

# endregion FEMModel















