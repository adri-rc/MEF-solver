"""Small utilities shared by the finite-element solvers."""


class DofMapper:
    """Map arbitrary GRID IDs to contiguous zero-based matrix indices."""

    def __init__(self, node_ids: list[int], dofs_per_node: int = 6):
        self.dofs_per_node = dofs_per_node
        self.node_to_index: dict[int, int] = {}
        self.index_to_node: dict[int, int] = {}
        
        if len(node_ids) != len(set(node_ids)):
            raise ValueError(f"Duplicate Node IDs found in input list: {node_ids}")
                
        # To ensure that the stiffness matrix is always constructed the same way,
        # regardless of the order in which you read the input file. 
        sorted_ids = sorted(node_ids) 
        
        for i, nid in enumerate(sorted_ids):
            self.node_to_index[nid] = i
            self.index_to_node[i] = nid
            
        self.total_dofs = len(node_ids) * dofs_per_node

    def get_dof_indices(self, node_id: int) -> list[int]:
        """Return the global DOF indices assigned to ``node_id``."""
        if node_id not in self.node_to_index:
            raise ValueError(f"Node ID {node_id} is not mapped.")
            
        start_index = self.node_to_index[node_id] * self.dofs_per_node
        return list(range(start_index, start_index + self.dofs_per_node))
    
    def get_node_from_index(self, global_index: int) -> int:
        """Return the GRID ID associated with a global matrix index."""
        if not (0 <= global_index < self.total_dofs):
            raise IndexError(f"Index {global_index} out of bounds (Total: {self.total_dofs})")
        
        internal_index = global_index // self.dofs_per_node
        return self.index_to_node[internal_index]
