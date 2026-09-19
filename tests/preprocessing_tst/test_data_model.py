"""Unit tests for the format-independent preprocessing data model."""

import numpy as np
import pytest

from src.preprocessing.data_model import (
    BarProperty,
    BarElement,
    CQUAD4Element,
    CTRIA3Element,
    FEMModel,
    IsotropicMaterial,
    Node,
    RodElement,
    ShellProperty,
    SinglePointLoad,
    SPCConstraint,
)

class TestNode:
    def test_initialization_defaults(self):
        node = Node(
            nid=1,
            coord_id=0,
            x1=10.0,
            x2=20.0,
            x3=30.0,
            cd_id=0,
            ps="",
            seid=0,
        )

        assert node.x == 0.0
        assert node.y == 0.0
        assert node.z == 0.0
        assert node.nid == 1

    def test_coordinate_storage(self):
        node = Node(
            nid=2,
            coord_id=0,
            x1=5.0,
            x2=0.0,
            x3=2.0,
            cd_id=0,
            ps="",
            seid=0,
        )

        node.x = node.x1
        node.y = node.x2
        node.z = node.x3

        coords = node.get_coords()

        assert isinstance(coords, np.ndarray), "The method must return a NumPy array"

        expected = np.array([5.0, 0.0, 2.0])

        assert np.allclose(coords, expected), f"Expected {expected}, got {coords}"


class TestIsotropicMaterial:
    def test_calc_G_from_E_nu(self):
        E_in = 210e9
        nu_in = 0.3

        mat = IsotropicMaterial(mid=1, E=E_in, nu=nu_in)

        expected_G = E_in / (2 * (1 + nu_in))

        assert mat.G == pytest.approx(expected_G, rel=1e-5)
        assert mat.E == E_in
        assert mat.nu == nu_in

    def test_calc_E_from_G_nu(self):
        G_in = 70e9
        nu_in = 0.35

        mat = IsotropicMaterial(mid=2, G=G_in, nu=nu_in)

        expected_E = 2 * G_in * (1 + nu_in)

        assert mat.E == pytest.approx(expected_E, rel=1e-5)

    def test_calc_nu_from_E_G(self):
        E_in = 1000.0
        G_in = 400.0

        mat = IsotropicMaterial(mid=3, E=E_in, G=G_in)

        expected_nu = E_in / (2 * G_in) - 1
        assert mat.nu == pytest.approx(expected_nu, rel=1e-5)

    def test_error_missing_parameters(self):
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=10, E=210e9)
        assert "must be at least 2 parameters" in str(excinfo.value)
        
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=10, G=20e9)
        assert "must be at least 2 parameters" in str(excinfo.value)
        
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=10, nu=0.3)
        assert "must be at least 2 parameters" in str(excinfo.value)

    def test_error_poisson_value_minus_one(self):
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=11, E=210e9, nu=-1.0)

        assert "Poisson's ratio cannot be -1" in str(excinfo.value)

    def test_error_zero_shear_modulus(self):
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=12, E=210e9, G=0.0)

        assert "G cannot be 0" in str(excinfo.value)

    def test_consistency_error(self):
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=20, E=200.0, nu=0.0, G=50.0)

        error_msg = str(excinfo.value)
        assert "physically inconsistent" in error_msg
        assert "Difference is > 5%" in error_msg

    def test_consistency_ok_no_warning(self, capsys: pytest.CaptureFixture[str]):
        IsotropicMaterial(mid=21, E=200.0, nu=0.0, G=100.0)

        captured = capsys.readouterr()

        assert captured.out == ""


class TestFEMModel:
    def test_initialization_empty(self):
        model = FEMModel()

        assert len(model.nodes) == 0
        assert len(model.materials) == 0
        assert len(model.elements) == 0
        assert len(model.properties) == 0
        assert len(model.constraints) == 0
        assert len(model.loads) == 0

    def test_add_and_get_node(self):
        model = FEMModel()

        new_node = Node(
            nid=50,
            coord_id=0,
            x1=10.0,
            x2=20.0,
            x3=30.0,
            cd_id=0,
            ps="",
            seid=0,
        )
        model.nodes[50] = new_node
        retrieved = model.get_node(50)

        assert retrieved == new_node
        assert retrieved.nid == 50
        assert retrieved.coord_id == 0
        assert retrieved.x1 == 10.0
        assert retrieved.x2 == 20.0
        assert retrieved.x3 == 30.0
        assert retrieved.cd_id == 0
        assert retrieved.ps == ""
        assert retrieved.seid == 0

    def test_get_node_not_found(self):
        model = FEMModel()
        model.nodes[1] = Node(
            nid=1,
            coord_id=0,
            x1=0.0,
            x2=0.0,
            x3=0.0,
            cd_id=0,
            ps="",
            seid=0,
        )
        id_searched = 999

        with pytest.raises(KeyError) as excinfo:
            model.get_node(id_searched)
        error_message = str(excinfo.value)
        assert "[FEMModel ERROR]" in error_message
        assert str(id_searched) in error_message
        assert "does not exist" in error_message

    def test_store_materials(self):
        model = FEMModel()

        mat = IsotropicMaterial(mid=100, E=210e9, nu=0.3)
        model.materials[100] = mat

        assert 100 in model.materials
        assert model.materials[100].E == 210e9
        assert len(model.materials) == 1

    def test_store_elements(self):
        model = FEMModel()

        elem = RodElement(eid=10, prop_id=5, nodes=[1, 2])
        model.elements[10] = elem

        assert 10 in model.elements
        assert model.elements[10].nodes == [1, 2]
        assert model.elements[10].prop_id == 5

    def test_store_loads(self):
        model = FEMModel()

        load1 = SinglePointLoad(
            sid=2, appnode_id=5, module=100.0, vector=[0, 0, 1]
        )
        load2 = SinglePointLoad(
            sid=2, appnode_id=6, module=200.0, vector=[0, 1, 0]
        )

        model.loads.append(load1)
        model.loads.append(load2)

        assert len(model.loads) == 2
        assert model.loads[0].module == 100.0
        assert model.loads[1].module == 200.0
        assert model.loads[0].sid == 2
        assert model.loads[1].sid == 2

    def test_store_constraints(self):
        model = FEMModel()

        spc1 = SPCConstraint(sid=1, node_id=10, components="123")
        spc2 = SPCConstraint(sid=1, node_id=11, components="456")

        model.constraints.append(spc1)
        model.constraints.append(spc2)

        assert len(model.constraints) == 2
        assert model.constraints[0].node_id == 10
        assert model.constraints[1].node_id == 11
        assert model.constraints[0].components == "123"
        assert model.constraints[1].components == "456"
class TestBarProperty:
    def test_full_initialization(self):
        stress_points = [
            (0.1, 0.1),
            (0.1, -0.1),
            (-0.1, -0.1),
            (-0.1, 0.1),
        ]

        prop = BarProperty(
            pid=10,
            mat_id=1,
            A=0.05,
            I1=1e-4,
            I2=2e-4,
            J=5e-5,
            NSM=1.0,
            I12=1e-6,
            K1=0.85,
            K2=0.85,
            stress_points=stress_points,
        )

        assert prop.pid == 10
        assert prop.I12 == 1e-6
        assert prop.K1 == 0.85
        assert prop.NSM == 1.0
        assert prop.J == 5e-5
        assert len(prop.stress_points) == 4
        assert prop.stress_points[1] == (0.1, -0.1)
        assert prop.stress_points[3] == (-0.1, 0.1)

    def test_physical_inconsistency(self):
        with pytest.raises(ValueError) as excinfo:
            BarProperty(
                pid=99,
                mat_id=1,
                A=1.0,
                J=1.0,
                I1=2.0,
                I2=2.0,
                I12=10.0,
            )
        assert "Inconsistent inertial properties" in str(excinfo.value)


class TestBarElement:
    def test_full_initialization(self):
        elem = BarElement(
            eid=200,
            prop_id=10,
            nodes=[1, 2],
            orientation_vector=[0.0, 1.0, 0.0],
            g0=None,
            offt="GGG",
            pa=45,
            pb=0,
            wa=[0.0, 0.05, 0.0],
            wb=[0.0, 0.05, 0.0],
        )

        assert elem.eid == 200
        assert elem.pa == 45
        assert elem.pb == 0
        assert elem.wa[1] == 0.05
        assert elem.orientation_vector == [0.0, 1.0, 0.0]


class TestShellProperty:
    def test_initialization(self):
        prop = ShellProperty(
            pid=1,
            mat_id=10,
            t=0.005,
            mat_id_bending=11,
            mat_id_shear=12,
        )
        assert prop.pid == 1
        assert prop.t == 0.005
        assert prop.mat_id == 10
        assert prop.mat_id_bending == 11
        assert prop.mat_id_shear == 12

    def test_blank_mid2_and_mid3_remain_disabled(self):
        prop = ShellProperty(pid=2, mat_id=20, t=0.1)
        assert prop.mat_id == 20
        assert prop.mat_id_bending is None
        assert prop.mat_id_shear is None


class TestCTRIA3Element:
    def test_initialization(self):
        elem = CTRIA3Element(eid=100, prop_id=1, nodes=[1, 2, 3], theta=30.0)
        assert elem.eid == 100
        assert elem.nodes == [1, 2, 3]
        assert elem.prop_id == 1
        assert elem.theta == 30.0

class TestCQUAD4Element:
    def test_initialization(self):
        elem = CQUAD4Element(eid=200, prop_id=2, nodes=[1, 2, 3, 4], theta=45.0)
        assert elem.eid == 200
        assert len(elem.nodes) == 4
        assert elem.nodes == [1, 2, 3, 4]
        assert elem.prop_id == 2
        assert elem.theta == 45.0


