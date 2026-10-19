
"""Assembly, solution, and result recovery for linear static analysis."""

import warnings

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve

from src.preprocessing.data_model import (
    BarElement,
    CQUAD4Element,
    CTRIA3Element,
    FEMModel,
    RodElement,
    SinglePointLoad,
    SinglePointMoment,
    SPCConstraint,
)

from .elements.bar_element import BarElementSolver
from .elements.cquad4_nastran_element import CQUAD4NastranElementSolver
from .elements.ctria3_element import CTRIA3ElementSolver
from .elements.rod_element import RodElementSolver
from .utilities import DofMapper


class LinearStaticSolver:
    """Solve the linear static equilibrium equation ``K u = F``."""

    DOFS_PER_NODE = 6

    def __init__(self, model: FEMModel):
        """Initialize a solver for a parsed finite-element model."""
        self.model = model
        self.k6rot = getattr(model, "k6rot", 0.0)

        node_ids = list(self.model.nodes.keys())
        self.dof_mapper = DofMapper(node_ids, dofs_per_node=self.DOFS_PER_NODE)
        self.n_dofs = self.dof_mapper.total_dofs

        # LIL supports efficient incremental assembly.
        self.K_global = lil_matrix((self.n_dofs, self.n_dofs))
        self.F_global = np.zeros(self.n_dofs)
        self.displacements = np.zeros(self.n_dofs)

    def _get_global_dof_indices(self, node_or_list) -> list[int]:
        """Return the global matrix indices for one node or a list of nodes."""
        if isinstance(node_or_list, list):
            indices = []
            for nid in node_or_list:
                indices.extend(self.dof_mapper.get_dof_indices(nid))
            return indices
        return self.dof_mapper.get_dof_indices(node_or_list)

    def assemble(self):
        """Build the global stiffness matrix and load vector once."""
        print(f"--- STARTING MODEL ASSEMBLY (Total DOFs: {self.n_dofs}) ---")
        self.model.validate()
        self.build_stiffness_matrix()
        self.build_load_vector()
        print("--- MODEL ASSEMBLY COMPLETED ---\n")
    
    def build_stiffness_matrix(self):
        """Compute and assemble every supported element stiffness matrix."""
        print("[Solver] Assembling Stiffness Matrix...")
        
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
            elif isinstance(element, CQUAD4Element):
                k_elem_global, node_ids = self._solve_cquad4_element(element)

            else:
                elem_type = type(element).__name__
                raise NotImplementedError(
                    f"[Solver Error] Element ID {eid} is of type '{elem_type}', "
                    "which is not currently implemented in this solver."
                )
            
            # Element solvers currently always return a stiffness matrix.
            if k_elem_global is not None:
                self._add_element_stiffness(node_ids, k_elem_global)
                count += 1
        
        print(f"[Solver] Assembled {count} elements.")
    
    def _solve_rod_element(self, element: RodElement):
        """Return one CROD global stiffness matrix and its node IDs."""
        n1 = self.model.get_node(element.nodes[0])
        n2 = self.model.get_node(element.nodes[1])
        prop = self.model.properties[element.prop_id]
        mat = self.model.materials[prop.mat_id]
        
        elem_solver = RodElementSolver(n1, n2, prop, mat)
        
        return elem_solver.get_stiffness_matrix(), [n1.nid, n2.nid]

    def _solve_bar_element(self, element: BarElement) -> tuple:
        """Return one CBAR global stiffness matrix and its node IDs."""
        n1 = self.model.get_node(element.nodes[0])
        n2 = self.model.get_node(element.nodes[1])
        prop = self.model.properties[element.prop_id]
        mat = self.model.materials[prop.mat_id]
        
        orientation = element.orientation_vector
        if element.g0 is not None:
            orientation = self.model.get_node(element.g0).get_coords() - n1.get_coords()
        if orientation is None or np.linalg.norm(orientation) < 1e-6:
             # The reader normally supplies a valid orientation.
             orientation = [0.0, 1.0, 0.0]

        solver = BarElementSolver(
            n1, n2, prop, mat, orientation=orientation,
            wa=element.wa, wb=element.wb, pa=element.pa, pb=element.pb,
        )
        return solver.get_stiffness_matrix(), [n1.nid, n2.nid]

    def _solve_ctria3_element(self, element: CTRIA3Element) -> tuple:
        """Return one CTRIA3 global stiffness matrix and its node IDs."""
        n1 = self.model.get_node(element.nodes[0])
        n2 = self.model.get_node(element.nodes[1])
        n3 = self.model.get_node(element.nodes[2])
        prop = self.model.properties[element.prop_id]
        mat = self.model.materials.get(prop.mat_id)
        
        mat_bending = None
        if hasattr(prop, 'mat_id_bending') and prop.mat_id_bending in self.model.materials:
             mat_bending = self.model.materials[prop.mat_id_bending]
        mat_shear = None
        if hasattr(prop, 'mat_id_shear') and prop.mat_id_shear in self.model.materials:
             mat_shear = self.model.materials[prop.mat_id_shear]

        solver = CTRIA3ElementSolver(
            n1,
            n2,
            n3,
            prop,
            mat,
            mat_bending=mat_bending,
            mat_shear=mat_shear,
            k6rot=self.k6rot,
            theta=element.theta,
        )
        
        return solver.get_stiffness_matrix(), [n1.nid, n2.nid, n3.nid]

    def _solve_cquad4_element(self, element: CQUAD4Element) -> tuple:
        """Build one CQUAD4 with the active Nastran reconstruction."""
        n1 = self.model.get_node(element.nodes[0])
        n2 = self.model.get_node(element.nodes[1])
        n3 = self.model.get_node(element.nodes[2])
        n4 = self.model.get_node(element.nodes[3])
        prop = self.model.properties[element.prop_id]
        
        # Membrane material (MID1).
        mat_membrane = self.model.materials.get(prop.mat_id)

        # Bending material (MID2).
        mat_bending = None
        if hasattr(prop, 'mat_id_bending') and prop.mat_id_bending in self.model.materials:
            mat_bending = self.model.materials[prop.mat_id_bending]

        # Transverse-shear material (MID3).
        mat_shear = None
        if hasattr(prop, 'mat_id_shear') and prop.mat_id_shear in self.model.materials:
            mat_shear = self.model.materials[prop.mat_id_shear]

        solver = CQUAD4NastranElementSolver(
            n1, n2, n3, n4, 
            prop, 
            mat_membrane=mat_membrane, 
            mat_bending=mat_bending,
            mat_shear=mat_shear,
            k6rot=self.k6rot,
            theta=element.theta,
        )
        return solver.get_stiffness_matrix(), [n1.nid, n2.nid, n3.nid, n4.nid]

    def _add_element_stiffness(self, node_ids: list, k_elem: np.ndarray):
        """Transform and scatter an element matrix into the global matrix."""
        global_indices = self._get_global_dof_indices(node_ids)
        nodal_transform = self._get_nodal_frame_transform(node_ids)
        k_elem = nodal_transform @ k_elem @ nodal_transform.T

        # Sparse incremental assembly.
        rows = len(global_indices)
        for i in range(rows):
            gi = global_indices[i]
            for j in range(rows):
                gj = global_indices[j]
                val = k_elem[i, j]
                if val != 0.0:
                    self.K_global[gi, gj] += val

    def _get_nodal_frame_transform(self, node_ids: list[int]) -> np.ndarray:
        """Map basic/global element vectors to the DOF frame of each GRID."""
        transform = np.zeros((6 * len(node_ids), 6 * len(node_ids)))
        for position, node_id in enumerate(node_ids):
            rotation = np.asarray(self.model.get_node(node_id).dof_transform)
            start = 6 * position
            transform[start:start+3, start:start+3] = rotation
            transform[start+3:start+6, start+3:start+6] = rotation
        return transform

    def _displacements_in_basic(self) -> np.ndarray:
        result = np.zeros_like(self.displacements)
        for node_id in self.model.nodes:
            indices = self._get_global_dof_indices(node_id)
            rotation = np.asarray(self.model.get_node(node_id).dof_transform)
            result[indices[:3]] = rotation.T @ self.displacements[indices[:3]]
            result[indices[3:]] = rotation.T @ self.displacements[indices[3:]]
        return result
    
    def build_load_vector(self):
        """Assemble supported nodal forces and moments into ``F_global``."""
        print("[Solver] Assembling Load Vector...")
        
        count_forces = 0
        count_moments = 0

        for load in self.model.loads:

            if isinstance(load, SinglePointLoad):
                
                nid = load.appnode_id
                scale_factor = load.module
                vec = np.array(load.vector)
                cid = load.cid
                
                if cid != 0:
                    raise ValueError(
                        f"[Solver Error] Force on Node {nid} uses CID={load.cid}. "
                        "This solver currently supports Global coordinates (CID=0) only."
                    )
                
                force_vec = scale_factor * vec
                force_vec = np.asarray(self.model.get_node(nid).dof_transform) @ force_vec
                
                # Forces act on translational DOFs Tx, Ty, and Tz.
                dof_indices = self._get_global_dof_indices(nid)
                self.F_global[dof_indices[0]] += force_vec[0] # Fx
                self.F_global[dof_indices[1]] += force_vec[1] # Fy
                self.F_global[dof_indices[2]] += force_vec[2] # Fz
                
                count_forces += 1
            
            elif isinstance(load, SinglePointMoment):
                nid = load.appnode_id
                scale = load.module
                vec = np.array(load.vector)
                
                if load.cid != 0:
                     raise ValueError(
                        f"[Solver Error] Moment on Node {nid} uses CID={load.cid}. "
                        "This solver currently supports Global coordinates (CID=0) only."
                    )
                
                moment_vec = scale * vec
                moment_vec = np.asarray(self.model.get_node(nid).dof_transform) @ moment_vec
                
                # Moments act on rotational DOFs Rx, Ry, and Rz.
                dof_indices = self._get_global_dof_indices(nid) 
                self.F_global[dof_indices[3]] += moment_vec[0] # Mx
                self.F_global[dof_indices[4]] += moment_vec[1] # My
                self.F_global[dof_indices[5]] += moment_vec[2] # Mz
                
                count_moments += 1
            
            else:
                load_type = type(load).__name__
                raise NotImplementedError(
                    f"[Solver Error] The load type '{load_type}' exists in the model "
                    f"but logic to process it is not yet implemented in 'build_load_vector'."
                )

        print(f"[Solver] Applied forces: {count_forces} Forces, {count_moments} Moments.")
    
    def solve(self) -> dict:
        """Solve ``K_ff u_f = F_f - K_fs u_s`` and recover results."""
        print("[Solver] Solving system...")
        
        free_dofs, fixed_dofs, u_s_dict = self._get_dof_partition()
        
        # CSC supports efficient column slicing and sparse factorization.
        K_csc = self.K_global.tocsc()
        
        # Partition the free/free and free/fixed blocks.
        K_ff = K_csc[free_dofs, :][:, free_dofs]

        K_fs = K_csc[free_dofs, :][:, fixed_dofs]
        
        F_f = self.F_global[free_dofs]
        
        u_s = np.zeros(len(fixed_dofs))
        for i, global_idx in enumerate(fixed_dofs):
            if global_idx in u_s_dict:
                u_s[i] = u_s_dict[global_idx]
        
        if not free_dofs:
            self.displacements[fixed_dofs] = u_s
            print("[Solver] Model is fully constrained; calculating reactions only.")
            return self.reactions_and_stresses()

        rhs = F_f - (K_fs @ u_s)
        
        if not np.all(np.isfinite(K_ff.data)) or not np.all(np.isfinite(rhs)):
            raise ValueError("[Solver Error] K_ff or the right-hand side contains NaN/Inf values.")
        asymmetry = K_ff - K_ff.T
        if asymmetry.nnz and np.max(np.abs(asymmetry.data)) > max(1.0, np.max(np.abs(K_ff.data))) * 1e-10:
            raise ValueError("[Solver Error] The assembled stiffness matrix is not symmetric.")
        row_norms = np.asarray(np.abs(K_ff).sum(axis=1)).ravel()
        zero_rows = np.flatnonzero(row_norms == 0.0)
        if zero_rows.size:
            global_zero_dofs = [free_dofs[index] for index in zero_rows]
            raise ValueError(
                "[Solver Error] The stiffness matrix is singular; free DOFs without stiffness were found: "
                f"{self._describe_dofs(global_zero_dofs)}"
            )

        print(f"[Solver] Solving linear system of size {K_ff.shape}...")

        try:
            # Promote singularity warnings to errors so mechanisms cannot pass silently.
            with warnings.catch_warnings():
                warnings.filterwarnings('error', category=MatrixRankWarning)
                u_f = spsolve(K_ff, rhs)
                
            # Some SciPy backends can return non-finite values without a warning.
            if not np.all(np.isfinite(u_f)):
                raise ValueError("Solver returned NaN/Inf values.")

        except (MatrixRankWarning, RuntimeError, ValueError) as e:
            raise ValueError(
                f"[Solver Error] The Stiffness Matrix is singular (unstable structure). "
                f"Check for mechanisms or insufficient boundary conditions.\nDetails: {str(e)}"
            )
        
        self.displacements[free_dofs] = u_f
        self.displacements[fixed_dofs] = u_s
        
        max_disp = np.max(np.abs(self.displacements))
        print(f"[Solver] Solution done. Max displacement: {max_disp:.4e}")

        return self.reactions_and_stresses()

    def _describe_dofs(self, indices: list[int]) -> list[str]:
        labels = ("T1", "T2", "T3", "R1", "R2", "R3")
        descriptions = []
        for index in indices[:20]:
            node_id = self.dof_mapper.index_to_node[index // self.DOFS_PER_NODE]
            descriptions.append(f"GRID {node_id} {labels[index % 6]}")
        if len(indices) > 20:
            descriptions.append(f"... and {len(indices) - 20} more")
        return descriptions

    def reactions_and_stresses(self) -> dict:
        """Recover reactions and element results from the solved displacement vector."""
        print("[Solver] Post-processing results (Reactions & Stresses)...")
        
        results = {
            "displacements": self.displacements.copy(),
            "reactions": {},
            "stresses": {}
        }

        # Reactions are relevant at constrained DOFs; free-DOF residuals should vanish.
        if self.model.constraints:
            K_csr = self.K_global.tocsr()
            F_internal = K_csr.dot(self.displacements)
            R_total = F_internal - self.F_global
            
            constrained_nodes = set(c.node_id for c in self.model.constraints if isinstance(c, SPCConstraint))
            
            for node_id in constrained_nodes:
                indices = self._get_global_dof_indices(node_id)
                # Reaction vector: [Fx, Fy, Fz, Mx, My, Mz].
                results["reactions"][node_id] = R_total[indices]

        displacements_basic = self._displacements_in_basic()
        for eid, element in self.model.elements.items():
            
            if not isinstance(element, (RodElement, BarElement, CTRIA3Element, CQUAD4Element)):

                elem_type = type(element).__name__
                raise NotImplementedError(
                    f"[Error] Element Type '{elem_type}' (EID: {eid}) is not supported for stress recovery. "
                    "Analysis stopped to prevent misleading results."
                )
            
            # CROD
            if isinstance(element, RodElement):
                n1 = self.model.get_node(element.nodes[0])
                n2 = self.model.get_node(element.nodes[1])
                prop = self.model.properties[element.prop_id]
                mat = self.model.materials[prop.mat_id]
                
                elem_solver = RodElementSolver(n1, n2, prop, mat)
                
                # Call the specific stress recovery method
                stress_data = elem_solver.calculate_stresses(displacements_basic, self.dof_mapper)
                stress_data['type'] = 'CROD'

                results["stresses"][eid] = stress_data
            
            # CBAR
            elif isinstance(element, BarElement):
                n1 = self.model.get_node(element.nodes[0])
                n2 = self.model.get_node(element.nodes[1])
                prop = self.model.properties[element.prop_id]
                mat = self.model.materials[prop.mat_id]
                
                orientation = element.orientation_vector
                if element.g0 is not None:
                    orientation = self.model.get_node(element.g0).get_coords() - n1.get_coords()
                if orientation is None:
                    orientation = [0.0, 1.0, 0.0]

                elem_solver = BarElementSolver(
                    n1, n2, prop, mat, orientation=orientation,
                    wa=element.wa, wb=element.wb, pa=element.pa, pb=element.pb,
                )

                stress_data = elem_solver.calculate_stresses(displacements_basic, self.dof_mapper)
                stress_data['type'] = 'CBAR'

                results["stresses"][eid] = stress_data

            # CTRIA3
            elif isinstance(element, CTRIA3Element):
                n1 = self.model.get_node(element.nodes[0])
                n2 = self.model.get_node(element.nodes[1])
                n3 = self.model.get_node(element.nodes[2])
                prop = self.model.properties[element.prop_id]
                mat = self.model.materials.get(prop.mat_id)
                
                mat_bending = None
                if hasattr(prop, 'mat_id_bending') and prop.mat_id_bending in self.model.materials:
                     mat_bending = self.model.materials[prop.mat_id_bending]
                
                mat_shear = None
                if hasattr(prop, 'mat_id_shear') and prop.mat_id_shear in self.model.materials:
                    mat_shear = self.model.materials[prop.mat_id_shear]
                
                elem_solver = CTRIA3ElementSolver(
                    n1,
                    n2,
                    n3,
                    prop,
                    mat,
                    mat_bending=mat_bending,
                    mat_shear=mat_shear,
                    k6rot=self.k6rot,
                    theta=element.theta,
                )

                stress_data = elem_solver.calculate_stresses(displacements_basic, self.dof_mapper)
                stress_data['type'] = 'CTRIA3' 

                results["stresses"][eid] = stress_data

            elif isinstance(element, CQUAD4Element):
                nodes = [self.model.get_node(nid) for nid in element.nodes]
                prop = self.model.properties[element.prop_id]
                mat_membrane = self.model.materials.get(prop.mat_id)
                mat_bending = self.model.materials.get(prop.mat_id_bending)
                mat_shear = self.model.materials.get(prop.mat_id_shear)
                elem_solver = CQUAD4NastranElementSolver(
                    *nodes, prop,
                    mat_membrane=mat_membrane,
                    mat_bending=mat_bending,
                    mat_shear=mat_shear,
                    k6rot=self.k6rot,
                    theta=element.theta,
                )
                stress_data = elem_solver.calculate_stresses(displacements_basic, self.dof_mapper)
                stress_data['type'] = 'CQUAD4'
                results["stresses"][eid] = stress_data

        
        return results 

    def _get_dof_partition(self):
        """Return free DOFs, fixed DOFs, and enforced displacement values."""
        fixed_dofs = set()
        enforced_displacements = {}
        str_to_index = {'1':0, '2':1, '3':2, '4':3, '5':4, '6':5}

        for constraint in self.model.constraints:
            
            # SPC (Single Point Constraint)
            if isinstance(constraint, SPCConstraint):
                node_global_indices = self._get_global_dof_indices(constraint.node_id)
                
                comps_str = str(constraint.components)
                
                for char in comps_str:
                    if char in str_to_index:
                        local_idx = str_to_index[char]     # 0 for '1'
                        global_idx = node_global_indices[local_idx]
                        
                        fixed_dofs.add(global_idx)
                        
                        if global_idx in enforced_displacements and not np.isclose(
                            enforced_displacements[global_idx], constraint.value,
                            rtol=1e-12, atol=1e-15,
                        ):
                            raise ValueError(
                                "Conflicting SPC values for "
                                f"{self._describe_dofs([global_idx])[0]}: "
                                f"{enforced_displacements[global_idx]} and {constraint.value}."
                            )
                        enforced_displacements[global_idx] = constraint.value
            
            else:
                constraint_type = type(constraint).__name__
                raise NotImplementedError(
                    f"[Solver Error] The constraint type '{constraint_type}' defined in the model is not supported. "
                    f"(Only SPCConstraint is implemented)."
                )
        
        all_dofs = set(range(self.n_dofs))
        free_dofs = sorted(list(all_dofs - fixed_dofs))
        fixed_dofs = sorted(list(fixed_dofs))
        
        return free_dofs, fixed_dofs, enforced_displacements
