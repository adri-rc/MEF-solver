
"""
@file test_utilities.py
@brief Unit tests for solver helper utilities.
@details
Focuses on DofMapper logic:
    - Verifies correct index calculation for sequential and non-sequential IDs.
    - Checks error handling for duplicate nodes or missing IDs.
    - Validates total DOF count.
"""

import pytest
from src.solver.utilities import DofMapper

class TestDofMapper:
    """
    @class TestDofMapper
    @brief Test suite for the Degree of Freedom mapping logic.
    """

    def test_basic_mapping_sequential(self):
        """
        @test Basic Sequential IDs.
        @brief Nodes 1, 2, 3 with 6 DOFs each.
        """
        nodes = [1, 2, 3]
        mapper = DofMapper(nodes, dofs_per_node=6)
        
        assert mapper.total_dofs == 18 # 3 * 6
        assert mapper.get_dof_indices(1) == [0, 1, 2, 3, 4, 5] # Node 1 -> Indices [0, 1, 2, 3, 4, 5]
        assert mapper.get_dof_indices(2) == [6, 7, 8, 9, 10, 11] # Node 2 -> Indices [6, 7, 8, 9, 10, 11]

    def test_sorting_non_sequential(self):
        """
        @test Sorting Logic.
        @brief Input: [100, 5]. The mapper MUST sort them internally.
        @details
        Matrix needs to be: Node 5 (rows 0-5), Node 100 (rows 6-11).
        """
        nodes = [502, 100, 5]
        mapper = DofMapper(nodes, dofs_per_node=6)

        indices_5 = mapper.get_dof_indices(5) # Node 5 should come first in the matrix
        assert indices_5[0] == 0
        indices_100 = mapper.get_dof_indices(100) # Node 100 should come second
        assert indices_100[0] == 6
        indices_502 = mapper.get_dof_indices(502) # Node 502 should come second
        assert indices_502[2] == 14

    def test_custom_dofs_per_node(self):
        """
        @test Custom DOF count (e.g., 2D Heat Transfer).
        @brief 1 DOF per node.
        """
        nodes = [10, 20]
        mapper = DofMapper(nodes, dofs_per_node=1)
        
        assert mapper.total_dofs == 2
        assert mapper.get_dof_indices(10) == [0]
        assert mapper.get_dof_indices(20) == [1]

    def test_error_missing_node(self):
        """
        @test Lookup Error.
        @brief Requesting a node that wasn't in the initialization list.
        """
        mapper = DofMapper([1, 2])
        with pytest.raises(ValueError) as excinfo:
            mapper.get_dof_indices(999)
        assert "not mapped" in str(excinfo.value)
    
    def test_error_duplicate_nodes(self):
        """
        @test Initialization Error.
        @brief Should fail if duplicate node IDs are provided.
        """
        duplicate_nodes = [1, 2, 1] # Node 1 appears twice        
        with pytest.raises(ValueError) as excinfo:
            DofMapper(duplicate_nodes)
        assert "Duplicate Node IDs" in str(excinfo.value)

    def test_reverse_lookup(self):
        """
        @test Get Node from Matrix Index.
        @brief For post-processing results.
        """
        # Nodes: 10, 20. DOFs: 6
        # Node 10: 0-5
        # Node 20: 6-11
        mapper = DofMapper([20, 10], dofs_per_node=6)
        
        # Index 0 belongs to Node 10 (because of sorting)
        assert mapper.get_node_from_index(0) == 10
        assert mapper.get_node_from_index(5) == 10
        
        # Index 7 belongs to Node 20
        assert mapper.get_node_from_index(7) == 20
        
        # Error out of bounds
        with pytest.raises(IndexError)as excinfo:
            mapper.get_node_from_index(99)
        error_msg = str(excinfo.value)
        assert "out of bounds" in error_msg
        assert "99" in error_msg