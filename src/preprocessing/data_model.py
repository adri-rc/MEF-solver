from dataclasses import dataclass, field
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
                error = abs(self.E - E_calculated) / self.E
                
                
                if error > 0.05: # If the error is greater than 5%, it´s notified (engineering tolerance).
                    raise ValueError(
                        f"Material ID '{self.mid}' is physically inconsistent. \n"
                        f"        Calculated E based on G/nu: {E_calculated:.2e}\n"
                        f"        User input E:               {self.E:.2e}\n"
                        f"        Difference is > 5% ({error:.1%}). Review material inputs."
                    )
                    

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
    
    def __post_init__(self):
        # physical consistency validation (Teorema de ejes paralelos/Mohr)
        # I1 * I2 must be greater than I12^2 for the section to be physically real.
        if (self.I1 * self.I2) < (self.I12 ** 2):
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
    mat_id: int          #: Material ID for Membrane
    t: float             #: Thickness
    mat_id_bending: int = -1 #: Material ID for Bending (Default: same as membrane if -1)
    mat_id_shear: int = -1   #: Material ID for Transverse Shear (MID3) (Default: same as MID2)

    def __post_init__(self):
        if self.mat_id_bending == -1:
            self.mat_id_bending = self.mat_id
        if self.mat_id_shear == -1:
            self.mat_id_shear = self.mat_id_bending
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
    @details Supports membrane and bending stiffness using Reissner-Mindlin and isoparametric formulation.
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
# endregion SPCConstraint

# endregion Constraints

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















