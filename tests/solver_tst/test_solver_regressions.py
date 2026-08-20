import numpy as np
import pytest

from src.preprocessing.data_model import (
    FEMModel, Node, IsotropicMaterial, RodProperty, ShellProperty,
    RodElement, CQUAD4Element, SinglePointLoad, SPCConstraint,
)
from src.solver.linear_static_solver import LinearStaticSolver


def _node(nid, x, y=0.0, z=0.0):
    node = Node(nid, 0, x, y, z, 0, "", 0)
    node.x, node.y, node.z = x, y, z
    return node


def test_fully_constrained_model_calculates_reactions():
    model = FEMModel()
    model.nodes = {1: _node(1, 0.0), 2: _node(2, 1.0)}
    model.materials[1] = IsotropicMaterial(1, E=1000.0, nu=0.3)
    model.properties[1] = RodProperty(1, 1, A=1.0, J=1.0, C=0.0, NSM=0.0)
    model.elements[1] = RodElement(1, 1, [1, 2])
    model.constraints = [
        SPCConstraint(1, 1, "123456"), SPCConstraint(1, 2, "123456")
    ]
    model.loads = [SinglePointLoad(1, 2, 10.0, [1.0, 0.0, 0.0])]

    solver = LinearStaticSolver(model)
    solver.assemble()
    results = solver.solve()
    assert results["reactions"][2][0] == pytest.approx(-10.0)


def test_conflicting_prescribed_displacements_are_rejected():
    model = FEMModel()
    model.nodes = {1: _node(1, 0.0), 2: _node(2, 1.0)}
    model.materials[1] = IsotropicMaterial(1, E=1000.0, nu=0.3)
    model.properties[1] = RodProperty(1, 1, A=1.0, J=1.0, C=0.0, NSM=0.0)
    model.elements[1] = RodElement(1, 1, [1, 2])
    model.constraints = [
        SPCConstraint(1, 1, "123456"),
        SPCConstraint(1, 2, "123456", value=0.0),
        SPCConstraint(2, 2, "1", value=0.1),
    ]
    solver = LinearStaticSolver(model)
    solver.assemble()
    with pytest.raises(ValueError, match="Conflicting SPC"):
        solver.solve()


def test_cquad4_is_supported_by_global_postprocessing():
    model = FEMModel()
    model.nodes = {
        1: _node(1, 0.0, 0.0), 2: _node(2, 1.0, 0.0),
        3: _node(3, 1.0, 1.0), 4: _node(4, 0.0, 1.0),
    }
    model.materials[1] = IsotropicMaterial(1, E=1000.0, nu=0.3)
    model.properties[1] = ShellProperty(1, 1, t=0.1)
    model.elements[1] = CQUAD4Element(1, 1, [1, 2, 3, 4])
    model.constraints = [SPCConstraint(1, nid, "123456") for nid in model.nodes]

    solver = LinearStaticSolver(model)
    solver.assemble()
    results = solver.solve()
    assert results["stresses"][1]["type"] == "CQUAD4"
    assert np.isfinite(results["stresses"][1]["von_mises"])
