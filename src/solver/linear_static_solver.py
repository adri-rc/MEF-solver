
"""
@file linear_static_solver.py
@brief Solver for Linear Static Analysis (Ku = F).
@details 
Orchestrates the Finite Element Analysis workflow:
1. DOF Mapping (Nodes -> Matrix Indices).
2. Global Stiffness Matrix Assembly.
3. Load Vector Assembly.
4. System Partitioning (Free vs. Constrained DOFs).
5. Numerical Solution (Sparse Solver).
Design Philosophy:
- FAIL FAST: The solver explicitly raises errors for any unsupported entity 
  (Element, Load, Constraint) preventing silent errors or physical inconsistencies.
- SEMANTIC ERRORS: Uses NotImplementedError for missing features and ValueError for invalid data.
"""

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import warnings
from scipy.sparse.linalg import spsolve, MatrixRankWarning

from .utilities import DofMapper

# --- ELEMENT SOLVERS ---
from .elements.rod_element_v2 import RodElementSolverV2
from .elements.bar_element import BarElementSolver
from .elements.ctria3_element import CTRIA3ElementSolver


# --- DATA MODEL ---
from src.preprocessing.data_model import (
    FEMModel, 
    RodElement,
    BarElement,
    CTRIA3Element,
    SinglePointLoad,  
    SinglePointMoment,
    SPCConstraint
)

class LinearStaticSolver:
    """
    @class LinearStaticSolver
    @brief Solves the linear static equilibrium equation: [K]{u} = {F}.
    """
    DOFS_PER_NODE = 6 # Assume 6 DOFs per node for a general Nastran-type solver.

    def __init__(self, model: FEMModel):
        """
        @brief Initialize the solver with a parsed FEM Model.
        @param model The database containing nodes, elements, materials, loads, etc.
        """
        self.model = model

        # 1. Initialize DOF Mapping
        # We convert keys to a list to ensure deterministic ordering.
        node_ids = list(self.model.nodes.keys())
        self.dof_mapper = DofMapper(node_ids, dofs_per_node=self.DOFS_PER_NODE)
        
        self.n_dofs = self.dof_mapper.total_dofs
        
        # 2. Initialize Matrices
        # 'lil_matrix' (list of list) is efficient for the incremental construction (Assembly phase)
        self.K_global = lil_matrix((self.n_dofs, self.n_dofs)) 
        self.F_global = np.zeros(self.n_dofs)
        
        # 3. Output Storage
        self.displacements = np.zeros(self.n_dofs)

    # ================= HELPER METHODS =================

    def _get_global_dof_indices(self, node_or_list) -> list[int]:
        """
        @brief Retrieves the global matrix indices corresponding to specific nodes.
        @param node_or_list Single Node ID (int) or a list of Node IDs.
        @return Flat list of 0-based global indices.
        """
        if isinstance(node_or_list, list):
            indices = []
            for nid in node_or_list:
                indices.extend(self.dof_mapper.get_dof_indices(nid))
            return indices
        else:
            return self.dof_mapper.get_dof_indices(node_or_list)

   
    # ================= ASSEMBLY PHASE =================

    def assemble(self):
        """
        @brief Master method to build the Global Stiffness Matrix (K) and Load Vector (F).
        """
        print(f"--- STARTING MODEL ASSEMBLY (Total DOFs: {self.n_dofs}) ---") 
        self.build_stiffness_matrix()
        self.build_load_vector()        
        print("--- MODEL ASSEMBLY COMPLETED ---\n")
    
    def build_stiffness_matrix(self):
        """
        @brief Iterates over all elements to compute and aggregate their stiffness.
        @throws NotImplementedError if an unsupported element type is encountered.
        """
        print(f"[Solver] Assembling Stiffness Matrix...")
        
        count = 0
        for eid, element in self.model.elements.items():
            
            k_elem_global = None
            node_ids = []

            if isinstance(element, RodElement):
                k_elem_global, node_ids = self._solve_rod_element(element)
            
            elif isinstance(element, BarElement):
                k_elem_global, node_ids = self._solve_bar_element(element)

            elif isinstance(element, CTRIA3Element):
                k_elem_global, node_ids = self._solve_ctria3_element(element)

            
            # Placeholder for future elements (Beam, Plate, Solid)
            # FUTURE: elif isinstance(element, BeamElement): ...

            else:
                # STRICT MODE: Fail fast if the model contains unknown elements.
                elem_type = type(element).__name__
                raise NotImplementedError(
                    f"[Solver Error] Element ID {eid} is of type '{elem_type}', "
                    "which is not currently implemented in this solver."
                )
            
            # Assemble if the element returned stifness
            if k_elem_global is not None:
                self._add_element_stiffness(node_ids, k_elem_global)
                count += 1
        
        print(f"[Solver] Assembled {count} elements.")
    
    def _solve_rod_element(self, element: RodElement):
        """
        @brief Computes the 12x12 stiffness matrix for a CROD element in global coordinates.
        """
        n1 = self.model.get_node(element.nodes[0])
        n2 = self.model.get_node(element.nodes[1])
        prop = self.model.properties[element.prop_id]
        mat = self.model.materials[prop.mat_id]
        
        # Delegate physics calculation to the specific element solver
        elem_solver = RodElementSolverV2(n1, n2, prop, mat)
        
        return elem_solver.get_stiffness_matrix(), [n1.nid, n2.nid]

    def _solve_bar_element(self, element: BarElement) -> tuple:
        """
        @brief Computes the 12x12 stiffness matrix for a CBAR element in global coordinates.
        """
        n1 = self.model.get_node(element.nodes[0])
        n2 = self.model.get_node(element.nodes[1])
        prop = self.model.properties[element.prop_id]
        mat = self.model.materials[prop.mat_id]
        
        # Validate Orientation (Crucial for Beam Physics)
        orientation = element.orientation_vector
        if not orientation or np.linalg.norm(orientation) < 1e-6:
             # Fallback to Y-axis if missing (though Reader usually checks this)
             orientation = [0.0, 1.0, 0.0]

        solver = BarElementSolver(n1, n2, prop, mat, orientation=orientation)
        return solver.get_stiffness_matrix(), [n1.nid, n2.nid]

    def _solve_ctria3_element(self, element: CTRIA3Element) -> tuple:
        """
        @brief Computes the 18x18 stiffness matrix for a CTRIA3 element in global coordinates.
        """
        n1 = self.model.get_node(element.nodes[0])
        n2 = self.model.get_node(element.nodes[1])
        n3 = self.model.get_node(element.nodes[2])
        prop = self.model.properties[element.prop_id]
        mat = self.model.materials[prop.mat_id]
        
        # Check for Bending Material (if property supports it)
        mat_bending = None
        # Check if attribute exists (it is a ShellProperty) and points to a valid material
        if hasattr(prop, 'mat_id_bending') and prop.mat_id_bending in self.model.materials:
             mat_bending = self.model.materials[prop.mat_id_bending]

        solver = CTRIA3ElementSolver(n1, n2, n3, prop, mat, mat_bending=mat_bending)
        
        return solver.get_stiffness_matrix(), [n1.nid, n2.nid, n3.nid]

    
    # --- CORE ASSEMBLY & SOLVE LOGIC ---

    def _add_element_stiffness(self, node_ids: list, k_elem: np.ndarray):
        """
        @brief Maps and adds the local element matrix into the global system matrix.
        """
        global_indices = self._get_global_dof_indices(node_ids)

        # Iterative assembly
        rows = len(global_indices)
        for i in range(rows):
            gi = global_indices[i]
            for j in range(rows):
                gj = global_indices[j]
                val = k_elem[i, j]
                if val != 0.0:
                    self.K_global[gi, gj] += val
    
    def build_load_vector(self):
        """
        @brief Aggregates nodal forces and moments into the global force vector (F).
        @throws NotImplementedError if an unsupported load type is found.
        @throws ValueError if the load definition is invalid (e.g. unsupported Coordinate System).
        """
        print("[Solver] Assembling Load Vector...")
        
        count_forces = 0
        count_moments = 0

        for load in self.model.loads:

            # --- CASE A: POINT FORCE ---
            if isinstance(load, SinglePointLoad):
                
                nid = load.appnode_id
                scale_factor = load.module
                vec = np.array(load.vector)
                cid = load.cid
                
                # Validation: Coordinate System
                if cid != 0:
                    raise ValueError(
                        f"[Solver Error] Force on Node {nid} uses CID={load.cid}. "
                        "This solver currently supports Global coordinates (CID=0) only."
                    )
                
                # Calculation
                force_vec = scale_factor * vec
                
                # Assembly: Force affects Translational DOFs (Tx, Ty, Tz) -> Indices 0, 1, 2
                dof_indices = self._get_global_dof_indices(nid)
                self.F_global[dof_indices[0]] += force_vec[0] # Fx
                self.F_global[dof_indices[1]] += force_vec[1] # Fy
                self.F_global[dof_indices[2]] += force_vec[2] # Fz
                
                count_forces += 1
            
            # --- CASE B: POINT MOMENT ---
            elif isinstance(load, SinglePointMoment):
                nid = load.appnode_id
                scale = load.module
                vec = np.array(load.vector)
                
                if load.cid != 0:
                     raise ValueError(
                        f"[Solver Error] Moment on Node {nid} uses CID={load.cid}. "
                        "This solver currently supports Global coordinates (CID=0) only."
                    )
                
                # Calculation
                moment_vec = scale * vec
                
                # Assembly: Moment affects Rotational DOFs (Rx, Ry, Rz) -> Indices 3, 4, 5
                dof_indices = self._get_global_dof_indices(nid) 
                self.F_global[dof_indices[3]] += moment_vec[0] # Mx
                self.F_global[dof_indices[4]] += moment_vec[1] # My
                self.F_global[dof_indices[5]] += moment_vec[2] # Mz
                
                count_moments += 1
            
            # --- CASE C: UNSUPPORTED ---
            else:
                load_type = type(load).__name__
                raise NotImplementedError(
                    f"[Solver Error] The load type '{load_type}' exists in the model "
                    f"but logic to process it is not yet implemented in 'build_load_vector'."
                )

        print(f"[Solver] Applied forces: {count_forces} Forces, {count_moments} Moments.")
    
    # ================= SOLUTION PHASE =================
    def solve(self)-> np.ndarray:
        """
        @brief Partitions the system and solves for unknown displacements.
        @details
        Solves [K_ff]{u_f} = {F_f} - [K_fs]{u_s}
        This approach handles both zero (clamped) and non-zero (enforced) displacements.
        @return A dictionary containing:
            - 'displacements': np.ndarray (Total Displacement Vector)
            - 'reactions': dict {NodeID: np.ndarray(6)}
            - 'stresses': dict {ElemID: dict(axial_force, stress...)}
        """
        print("[Solver] Solving system...")
        
        # 1. Identify Boundary Conditions
        free_dofs, fixed_dofs, u_s_dict = self._get_dof_partition()
        
        if not free_dofs:
            raise ValueError("Model is fully constrained or empty. Nothing to solve.")

        # 2. Optimization: Convert LIL (List of List) to CSC (Compressed Sparse Column)
        # CSC is efficient for column slicing (needed for K_fs) and linear algebra.
        K_csc = self.K_global.tocsc()
        
        # 3. Partitioning
        # K_ff: Interaction between Free DOFs (The system to solve)
        K_ff = K_csc[free_dofs, :][:, free_dofs]

        # K_fs: Interaction between Free and Fixed DOFs
        # (Necessary to calculate reaction forces from enforced displacements)
        K_fs = K_csc[free_dofs, :][:, fixed_dofs]
        
        # 4. Partition Force Vector
        F_f = self.F_global[free_dofs]
        
        # 5. Build Enforced Displacement Vector (u_s)
        u_s = np.zeros(len(fixed_dofs))
        for i, global_idx in enumerate(fixed_dofs):
            if global_idx in u_s_dict:
                u_s[i] = u_s_dict[global_idx]
        
        # 6. Calculate Effective Load Vector (RHS)
        # Equation: K_ff * u_f + K_fs * u_s = F_f  =>  K_ff * u_f = F_f - K_fs * u_s
        # If u_s is all zeros, the term K_fs * u_s vanishes (standard clamping).
        rhs = F_f - (K_fs @ u_s) # subtract the effect of the imposed displacement.
        
        # 7. Numerical Solution
        print(f"[Solver] Solving linear system of size {K_ff.shape}...")

        # SECURITY BLOCK FOR SINGULAR MATRIX
        try:
            # We convert MatrixRankWarning into a real exception.
            with warnings.catch_warnings():
                warnings.filterwarnings('error', category=MatrixRankWarning)
                u_f = spsolve(K_ff, rhs)
                
            # Double check: If Scipy returned NaNs without warning (happens in some backends)
            if np.any(np.isnan(u_f)):
                raise ValueError("Solver returned NaNs.")

        except (MatrixRankWarning, RuntimeError, ValueError) as e:
            raise ValueError(
                f"[Solver Error] The Stiffness Matrix is singular (unstable structure). "
                f"Check for mechanisms or insufficient boundary conditions.\nDetails: {str(e)}"
            )
        
        # 8. Reconstruct Full Vector
        self.displacements[free_dofs] = u_f
        self.displacements[fixed_dofs] = u_s
        
        max_disp = np.max(np.abs(self.displacements))
        print(f"[Solver] Solution done. Max displacement: {max_disp:.4e}")

        # 9. Post-Process (Reactions and Stresses)
        return self.reactions_and_stresses()

    def reactions_and_stresses(self) -> dict:
        """
        @brief Calculates reaction forces and element stresses based on the solution.
        @details
        Reactions: R = K_global * u_total - F_external
        Stresses: Delegates to each element solver (e.g. RodElementSolverV2.calculate_stresses)

        @return Dictionary with structure:
            {
                "displacements": np.ndarray,
                "reactions": { NodeID: np.ndarray(6) },
                "stresses": { 
                    ElemID: { 
                    'type': str('CROD', 'CBAR')
                    'axial_force': float, 
                    ...
                            }
                            }
            }
        """
        print("[Solver] Post-processing results (Reactions & Stresses)...")
        
        results = {
            "displacements": self.displacements,
            "reactions": {},
            "stresses": {}
        }

        # --- A. REACTION FORCES ---
        # R = K * u - F
        # Only relevant at constrained DOFs. At free DOFs, R should be ~0 (equilibrium check).
        if self.model.constraints:
            K_csr = self.K_global.tocsr()
            F_internal = K_csr.dot(self.displacements)
            R_total = F_internal - self.F_global
            
            # Filter and store reactions only for constrained nodes
            # Optimization: We loop through unique constrained nodes to avoid duplicate work
            constrained_nodes = set(c.node_id for c in self.model.constraints if isinstance(c, SPCConstraint))
            
            for node_id in constrained_nodes:
                indices = self._get_global_dof_indices(node_id)
                # Store the reaction vector [Rx, Ry, Rz, Mx, My, Mz]
                results["reactions"][node_id] = R_total[indices]

        # --- B. ELEMENT STRESSES ---
        for eid, element in self.model.elements.items():
            
            # ---------------------------------------------------------
            # 0. FAIL-FAST CHECK
            # ---------------------------------------------------------
            # If the element is not supported -> abort
            if not isinstance(element, (RodElement, BarElement, CTRIA3Element)):

                elem_type = type(element).__name__
                raise NotImplementedError(
                    f"[Error] Element Type '{elem_type}' (EID: {eid}) is not supported for stress recovery. "
                    "Analysis stopped to prevent misleading results."
                )
            
            # ---------------------------------------------------------
            # 1. CASO CROD (Rod Element)
            # ---------------------------------------------------------
            if isinstance(element, RodElement):
                # TODO: Optimization - Retrieve cached ElementSolver instead of re-instantiating.
                # Currently, geometry and T matrix are re-computed here.
                n1 = self.model.get_node(element.nodes[0])
                n2 = self.model.get_node(element.nodes[1])
                prop = self.model.properties[element.prop_id]
                mat = self.model.materials[prop.mat_id]
                
                elem_solver = RodElementSolverV2(n1, n2, prop, mat)
                
                # Call the specific stress recovery method
                stress_data = elem_solver.calculate_stresses(self.displacements, self.dof_mapper)
                stress_data['type'] = 'CROD'  # Etiqueta identificativa

                results["stresses"][eid] = stress_data
            
            # ---------------------------------------------------------
            # 2. CASO CBAR (Bar Element)
            # ---------------------------------------------------------
            elif isinstance(element, BarElement):
                n1 = self.model.get_node(element.nodes[0])
                n2 = self.model.get_node(element.nodes[1])
                prop = self.model.properties[element.prop_id]
                mat = self.model.materials[prop.mat_id]
                
                orientation = element.orientation_vector if element.orientation_vector else [0.0, 1.0, 0.0]

                elem_solver = BarElementSolver(n1, n2, prop, mat, orientation=orientation)

                stress_data = elem_solver.calculate_stresses(self.displacements, self.dof_mapper)
                stress_data['type'] = 'CBAR'  # Etiqueta identificativa

                results["stresses"][eid] = stress_data

            # ---------------------------------------------------------
            # 3. CASO CTRIA3 (Shell Element)
            # ---------------------------------------------------------
            elif isinstance(element, CTRIA3Element):
                n1 = self.model.get_node(element.nodes[0])
                n2 = self.model.get_node(element.nodes[1])
                n3 = self.model.get_node(element.nodes[2])
                prop = self.model.properties[element.prop_id]
                mat = self.model.materials[prop.mat_id]
                
                # Check for Bending Material (if property supports it)
                mat_bending = None
                # Check if attribute exists (it is a ShellProperty) and points to a valid material
                if hasattr(prop, 'mat_id_bending') and prop.mat_id_bending in self.model.materials:
                     mat_bending = self.model.materials[prop.mat_id_bending]

                elem_solver = CTRIA3ElementSolver(n1, n2, n3, prop, mat, mat_bending=mat_bending)
                
                stress_data = elem_solver.calculate_stresses(self.displacements, self.dof_mapper)
                stress_data['type'] = 'CTRIA3' 

                results["stresses"][eid] = stress_data

        
        return results 

    def _get_dof_partition(self):
        """
        @brief Separate DOFs into Free (f) and Restricted (s).
        @return (free_dofs, fixed_dofs, enforced_displacements).
        @raises NotImplementedError If it encounters an unsupported constraint type (e.g., MPC).
        """
        fixed_dofs = set()
        enforced_displacements = {}
        str_to_index = {'1':0, '2':1, '3':2, '4':3, '5':4, '6':5}

        for constraint in self.model.constraints:
            
            # SPC (Single Point Constraint)
            if isinstance(constraint, SPCConstraint):
                # 1. Obtain global base indexes for the node (e.g., [60, 61... 65])
                node_global_indices = self._get_global_dof_indices(constraint.node_id)
                
                # 2. Parse SPC components. Example: “123” -> indices 0, 1, 2
                comps_str = str(constraint.components)
                
                for char in comps_str:
                    if char in str_to_index:
                        local_idx = str_to_index[char]     # 0 for '1'
                        global_idx = node_global_indices[local_idx]
                        
                        fixed_dofs.add(global_idx)
                        
                        # If there is imposed displacement (not zero)
                        if constraint.value != 0.0:
                            enforced_displacements[global_idx] = constraint.value
            
            # ERROR - UNKNOWN CONSTRAIN 
            else:
                constraint_type = type(constraint).__name__
                raise NotImplementedError(
                    f"[Solver Error] The constraint type '{constraint_type}' defined in the model is not supported. "
                    f"(Only SPCConstraint is implemented)."
                )
        
        # Calculate free DOFs (Total - Fixed)
        all_dofs = set(range(self.n_dofs))
        free_dofs = sorted(list(all_dofs - fixed_dofs))
        fixed_dofs = sorted(list(fixed_dofs))
        
        return free_dofs, fixed_dofs, enforced_displacements