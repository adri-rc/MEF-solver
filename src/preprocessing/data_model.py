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
class NODE:
    """
    @class Node
    @brief Generic node definition where raw data from input file can be saved and general 
        coordinates in global sistem are also saved.
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
    

# region MATERIALS CLASES
@dataclass
class MATERIAL:
    """
    @class MATERIAL
    @brief Base class for properties constituting materials.
    """
    mid: int #: ID material

# region ISOTROPIC_MATERIAL(MATERIAL)
@dataclass
class ISOTROPIC_MATERIAL(MATERIAL):
    """
    @class ISOTROPIC_MATERIAL
    @brief Linear isotropic elastic material (equivalent to NASTRAN MAT1).
    
    Allows the class to be instantiated with ANY pair of constants:
        1. (E, nu) -> Calculates G
        2. (E, G) -> Calculates nu
        3. (G, nu) -> Calculates E

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
                    print(f"Review inputs in material '{self.mid}'\n"
                          "The parameters entered are not consistent. "
                        f"E theoretical: {E_calculated:.2e}, E introduced: {self.E:.2e}")

    def __repr__(self):
        # Usamos 'or 0' por seguridad en el f-string, aunque post_init asegura valores
        return (f"IsotropicMaterial with id:{self.mid}\n"
                f"(E={self.E:.2e}, nu={self.nu:.4f}, G={self.G:.2e}, rho={self.rho})")
# endregion ISOTROPIC_MATERIAL(MATERIAL)

# endregion MATERIAL CLASES

# region PROPERTY CLASES
@dataclass
class PROPERTY:
    """
    @class PPROPERTY
    @brief Base class for any elemen property.
    """
    pid: int #: ID property
    

# region ROD_PROPERTY(PPROPERTY)
@dataclass
class ROD_PROPERTY(PROPERTY):
    """
    @class ROD_PROPERTY
    @brief Corresponds to PROD in NASTRAN.
    """
    mat_id: int   #: ID from asociated material.
    A: float = 0.0    #: Transverse section area.
    J: float = 0.0    #: Torsion constant (Polar inertia for circulars coordinates)
    C: float = 0.0    #: Coefficient for calculating torsional stresses
    NSM: float = 0.0  #: Non-structural mass per unit length
# endregion ROD_PROPERTY(PPROPERTY)
# endregion PROPERTY CLASES


# region ELEMENT CLASES
@dataclass
class ELEMENT:
    """Base class for elements."""
    eid: int
    prop_id: int
    nodes: List[int]

@dataclass
class ROD_ELEMENT(ELEMENT):
    """
    @brief 2-node bar (CROD).
    @details Defines a tension-compression-torsion element"""
    pass

@dataclass
class QUADELEMENT(ELEMENT):
    """
    @brief 4-node Quadrilateral Plate Element (CQUAD4).
    @details ¿Supports flat tension and plate bending.?
    """
    pass

# endregion ELEMENT CLASES


# region LOADs CLASES
@dataclass
class LOAD_CASE:
    """
    @class LOAD_CASE
    @brief Abstract base class for any load type (point, distributed, thermal).
    """
    sid: int  #: Set ID (SID) para agrupar casos de carga

# region SINGLE_POINT_LOAD(LOAD_CASE) (FORCE)
@dataclass
class SINGLE_POINT_LOAD(LOAD_CASE):
    """
    @class SINGLE_POINT_LOAD
    @brief Represents a static point load (FORCE card).
    @details Defines a force applied to a specific node, 
        with a magnitude and direction vector referenced to a coordinate system.
    """
    appnode_id: int      #: Application node ID (nid)
    module: float        #: Force module (F, scale factor)
    vector: List[float]  #: Director vector [N1, N2, N3]
    cid: int = 0    #: Coordinate system ID (CID) where the vector is defined. 0=Global.
# endregion SINGLE_POINT_LOAD(LOAD_CASE) (FORCE)

# region SINGLE_MOMENT_LOAD(LOAD_CASE) (MOMENT)
@dataclass
class SINGLE_MOMENT_LOAD(LOAD_CASE):
    """
    @class SINGLE_MOMENT_LOAD
    @brief Point-in-time applied to a node (MOMENT card).
    @details affects the degrees of rotational freedom. (For beams(CBAR) or plates (CQUAD4) elements.)
    """
    appnode_id: int      #: Application node ID (nid)
    module: float        #: Moment module (M, scale factor)
    vector: List[float]  #: Director vector [N1, N2, N3]
    cid: int = 0         #: Coordinate system ID (CID) where the vector is defined. 0=Global.
# endregion SINGLE_MOMENT_LOAD(LOAD_CASE) (MOMENT)

# endregion LOADS CLASES

# region CONSTRAINTs
@dataclass
class CONSTRAINT:
    """
    @class CONSTRAINTS
    @brief Base class for boundary conditions and kinematic constraints.
    """
    sid: int  #: Set ID (SID) para agrupar restricciones (Load Case)

# region SPC_CONSTRAINT
@dataclass
class SPC_CONSTRAINT(CONSTRAINT):
    """
    @class SPC_CONSTRAINT
    @brief Restricción de Punto Único (Tarjetas SPC y SPC1).
    @details Fija grados de libertad específicos de un nodo a un valor constante.
    """
    node_id: int        #: Grid ID (Gi)
    components: str     #: Componentes restringidas (ej. "123456" o "3")
    value: float = 0.0  #: Enforced displacement (Di). Default=0.0 (empotramiento).
# endregion SPC_CONSTRAINT

# endregion CONSTRAINTs

# region FEM_MODEL
class FEM_MODEL:
    """
    @class FEM_MODEL
    @brief Main container (“Database”) of the physical model.
    @details Stores all parsed and processed entities, 
        acting as the single source for the Solver.
    """
    def __init__(self):
        self.nodes: Dict[int, NODE] = {}            #: Diccionario de Nodos {ID: NODE}
        self.materials: Dict[int, MATERIAL] = {}    #: Diccionario de Materiales {ID: MATERIAL}
        self.properties: Dict[int, PROPERTY] = {}   #: Propiedades de sección {ID: PROPERTY}
        self.elements: Dict[int, ELEMENT] = {}      #: Elementos finitos {ID: ELEMENT}

                # Future user coordinate systems will be implemented
                #self.coords: Dict[int, CoordSystem] = {}    #: Sistemas de Coordenadas definidos por usuario
        
        # Loads and restrictions are stored in lists because different sucase can count with the same load_id or SPC_id multiple times,
        # therefor in a dictionary they can be overwritten.
        self.constraints: List[CONSTRAINT] = []     
        self.loads: List[LOAD_CASE] = []

    def get_node(self, node_id: int) -> NODE:
        """
        @brief Retrieves a NODE object by its ID.
        @param node_id ID of the node to search for.
        @return Class NODE.
        @throws KeyError if the node does not exist in the model.
        """
        return self.nodes[node_id]
    
        '''def get_coord_system(self, cid: int) -> Optional[CoordSystem]:
            """
            @brief Recupera un sistema de coordenadas.
            @return Objeto CoordSystem o None si es el global (CID=0).
            """
            if cid == 0:
                return None
            return self.coords[cid]'''

# endregion FEM_MODEL















