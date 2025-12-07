
"""
@file utilities.py
@brief Helper classes for the Finite Element Solver.
"""

class DofMapper:
    """
    @class DofMapper
    @brief Maps arbitrary Node IDs to continuous Global Matrix Indices.
    @details
    Nastran allows arbitrary IDs (e.g., 10, 50, 200). The solver matrices 
    need continuous 0-based indexing (0, 1, 2...).
    
    This class handles the translation:
    Global_Index = (Node_Internal_Index * DOFs_per_node) + DOF_Component
    """
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
        """
        @brief Returns the range of global indices for a given node.
        @param node_id The User ID of the node (e.g., 101).
        @return List of integers (e.g., [12, 13, 14, 15, 16, 17])
        """
        if node_id not in self.node_to_index:
            raise ValueError(f"Node ID {node_id} is not mapped.")
            
        start_index = self.node_to_index[node_id] * self.dofs_per_node
        return list(range(start_index, start_index + self.dofs_per_node))
    
    def get_node_from_index(self, global_index: int) -> int:
        """
        @brief Reverse lookup: Given a matrix row/col index, find the Node ID.
        @details Useful for debugging or mapping results back to nodes.
        """
        if not (0 <= global_index < self.total_dofs):
            raise IndexError(f"Index {global_index} out of bounds (Total: {self.total_dofs})")
        
        internal_index = global_index // self.dofs_per_node # integer division (//)
        return self.index_to_node[internal_index]