"""Tests for global degree-of-freedom indexing utilities."""

import pytest

from src.solver.utilities import DofMapper


class TestDofMapper:
    def test_basic_mapping_sequential(self):
        nodes = [1, 2, 3]
        mapper = DofMapper(nodes, dofs_per_node=6)

        assert mapper.total_dofs == 18
        assert mapper.get_dof_indices(1) == [0, 1, 2, 3, 4, 5]
        assert mapper.get_dof_indices(2) == [6, 7, 8, 9, 10, 11]

    def test_sorting_non_sequential(self):
        nodes = [502, 100, 5]
        mapper = DofMapper(nodes, dofs_per_node=6)

        indices_5 = mapper.get_dof_indices(5)
        assert indices_5[0] == 0
        indices_100 = mapper.get_dof_indices(100)
        assert indices_100[0] == 6
        indices_502 = mapper.get_dof_indices(502)
        assert indices_502[2] == 14

    def test_custom_dofs_per_node(self):
        nodes = [10, 20]
        mapper = DofMapper(nodes, dofs_per_node=1)

        assert mapper.total_dofs == 2
        assert mapper.get_dof_indices(10) == [0]
        assert mapper.get_dof_indices(20) == [1]

    def test_error_missing_node(self):
        mapper = DofMapper([1, 2])
        with pytest.raises(ValueError) as excinfo:
            mapper.get_dof_indices(999)
        assert "not mapped" in str(excinfo.value)

    def test_error_duplicate_nodes(self):
        duplicate_nodes = [1, 2, 1]
        with pytest.raises(ValueError) as excinfo:
            DofMapper(duplicate_nodes)
        assert "Duplicate Node IDs" in str(excinfo.value)

    def test_reverse_lookup(self):
        mapper = DofMapper([20, 10], dofs_per_node=6)

        assert mapper.get_node_from_index(0) == 10
        assert mapper.get_node_from_index(5) == 10

        assert mapper.get_node_from_index(7) == 20

        with pytest.raises(IndexError) as excinfo:
            mapper.get_node_from_index(99)
        error_msg = str(excinfo.value)
        assert "out of bounds" in error_msg
        assert "99" in error_msg
