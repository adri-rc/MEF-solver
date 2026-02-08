import pytest
import numpy as np
from src.preprocessing.data_model import (          
    Node,
    IsotropicMaterial, # Material
    FEMModel,
    Property,
    BarProperty,
    RodElement,        # Element
    BarElement,
    SPCConstraint,     # Constraint
    SinglePointLoad,   # LoadCase
    ShellProperty,
    CTRIA3Element
    )


"""
@file test_data_model.py
@brief Unit tests for the main data clasess defined in 'data_model.py'.
"""

class TestNode:
    """
    @class TestNode
    @brief Suite de pruebas para validar la entidad NODE.
    @details
    Followings issues are covered:
        1. Verifica la inicialización
        2. Almacenamiento de datos raw
    """

    def test_initialization_defaults(self):
        """
        @test Validation of default values during initialization.
        @brief Verify that when creating a node, the processed coordinates (x, y, z) start at 0.0.
        """
        node = Node(nid=1, coord_id=0, x1=10.0, x2=20.0, x3=30.0, cd_id=0, ps="", seid=0)

        assert node.x == 0.0
        assert node.y == 0.0
        assert node.z == 0.0
        assert node.nid == 1

    def test_coordinate_storage(self):
        """
        @test Storage and retrieval of coordinates.
        @brief Verify that the get_coords() method returns a correct 
            numpy array after assigning global coordinates.
        """
        node = Node(nid=2, coord_id=0, x1=5.0, x2=0.0, x3=2.0, cd_id=0, ps="", seid=0)
        
        node.x = node.x1
        node.y = node.x2
        node.z = node.x3

        coords = node.get_coords()

        assert isinstance(coords, np.ndarray), "El método debe devolver un numpy array"
        
        expected = np.array([5.0, 0.0, 2.0])
        
        # Nota: Para flotantes es mejor usar np.allclose para evitar errores de precisión decimal
        assert np.allclose(coords, expected), f"Se esperaba {expected}, se obtuvo {coords}"


class TestIsotropicMaterial:
    """
    @class TestIsotropicMaterial
    @brief Test suite for the IsotropicMaterial class.
    @details    
    Followings issues are covered:
        1. Cálculo automático de G dado (E, nu).
        2. Cálculo automático de E dado (G, nu).
        3. Cálculo automático de nu dado (E, G).
        4. Detección de inconsistencias físicas (cuando se dan los 3 parámetros).
        5. Manejo de errores matemáticos (división por cero).
    

    """

    # --- CASOS DE ÉXITO (Cálculos Automáticos) ---

    def test_calc_G_from_E_nu(self):
        """
        @test Cálculo del Módulo de Cortante (G).
        @brief Checks: G = E / (2 * (1 + nu)).
        """
        E_in = 210e9  # 210 GPa (Fe)
        nu_in = 0.3   # 0.3     (Fe)
        
        mat = IsotropicMaterial(mid=1, E=E_in, nu=nu_in)
        
        # Valor teórico esperado
        expected_G = E_in / (2 * (1 + nu_in)) # ~80.7 GPa
        
        assert mat.G == pytest.approx(expected_G, rel=1e-5)
        assert mat.E == E_in
        assert mat.nu == nu_in

    def test_calc_E_from_G_nu(self):
        """
        @test Cálculo del Módulo de Young (E).
        @brief Checks: E = 2 * G * (1 + nu).
        """
        # Caso: Aluminio (aprox)
        G_in = 70e9   # 70 GPa (Al)
        nu_in = 0.35  # 0.35 (Al)
        
        mat = IsotropicMaterial(mid=2, G=G_in, nu=nu_in)
        
        expected_E = 2 * G_in * (1 + nu_in)
        
        assert mat.E == pytest.approx(expected_E, rel=1e-5)

    def test_calc_nu_from_E_G(self):
        """
        @test Cálculo del Coeficiente de Poisson (nu).
        @brief Checks: nu = (E / (2 * G)) - 1.
        """
        
        E_in = 1000.0
        G_in = 400.0  
        
        mat = IsotropicMaterial(mid=3, E=E_in, G=G_in)
        
        expected_nu = E_in/(2*G_in)-1 # -> nu = 0.25
        assert mat.nu == pytest.approx(expected_nu, rel=1e-5)

    # --- CASOS DE ERROR Y VALIDACIÓN (Edge Cases) ---

    def test_error_missing_parameters(self):
        """
        @test Error due to lack of parameters.
        @brief Fail if fewer than 2 elastic constants are provided.
        """
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=10, E=210e9) # Missing nu or G
        assert "must be at least 2 parameters" in str(excinfo.value)
        
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=10, G=20e9) # Missing nu or E    
        assert "must be at least 2 parameters" in str(excinfo.value)
        
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=10, nu=0.3) # Missing E or G
        assert "must be at least 2 parameters" in str(excinfo.value)

    def test_error_poisson_value_minus_one(self):
        """
        @test Mathematical error due to singularity.
        @brief If nu = -1, G's formula involves dividing by zero.
                G = E / (2*(1 + (-1))) -> División por cero
        """
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=11, E=210e9, nu=-1.0)
        
        assert "Poisson's ratio cannot be -1" in str(excinfo.value)

    def test_error_zero_shear_modulus(self):
        """
        @test Error if G=0 when attempting to calculate nu.
        @brief nu = (E/2G) - 1. If G=0 -> Error
        """
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=12, E=210e9, G=0.0)
        
        assert "G cannot be 0" in str(excinfo.value)

    # --- CASOS DE CONSISTENCIA (Warning) ---

    def test_consistency_error(self):
        """
        @test Physical consistency error.
        @brief If the user provides E, G, and nu, but they do not match the formula,
                a ValueError is raised (Fail Fast).
        """
        # Non real material (E=200 y nu=0.0-> Should have G=100.) but G is forced to 50. (50% Error)
    
        with pytest.raises(ValueError) as excinfo:
            IsotropicMaterial(mid=20, E=200.0, nu=0.0, G=50.0)
        
        # Check if expected error is raised
        error_msg = str(excinfo.value)
        assert "physically inconsistent" in error_msg
        assert "Difference is > 5%" in error_msg

    def test_consistency_ok_no_warning(self, capsys: pytest.CaptureFixture[str]):
        """
        @test Silent verification when the data is correct.
        @brief If we provide the three parameters and they are consistent, it should NOT print anything.
        """
        IsotropicMaterial(mid=21, E=200.0, nu=0.0, G=100.0)
        
        captured = capsys.readouterr()
        
        assert captured.out == "" # The output must be blank.


class TestFEMModel:
    """
    @class TestFEMModel
    @brief Test suite for the model's management class.
    @details
    Followings issues are covered:
        1. Inicialización vacía.
        2. Inserción y recuperación de Nodos (get_node).
        3. Manejo de errores al buscar entidades inexistentes.
        4. Almacenamiento correcto de Elementos, Materiales, Cargas y Restricciones.

    """

    def test_initialization_empty(self):
        """
        @test Initial state of the model.
        @brief Checks that when creating a FEMModel, all its containers are empty.
        """
        model = FEMModel()
        
        # Dictionaries and lists must be empty
        assert len(model.nodes) == 0
        assert len(model.materials) == 0
        assert len(model.elements) == 0
        assert len(model.properties) == 0
        assert len(model.constraints) == 0
        assert len(model.loads) == 0

    # --- TEST DE MÉODOS (Método get_node) ---

    def test_add_and_get_node(self):
        """
        @test Test for 'get_node' method.
        @brief Add and retrieve a node. Checks that the get_node(nid) method returns the correct object.
        """
        model = FEMModel()
        
        new_node = Node(nid=50, coord_id=0, x1=10., x2=20., x3=30., cd_id=0, ps="", seid=0)
        model.nodes[50] = new_node # Manual insertion (simulating what the nastran_reader.py does)
        retrieved = model.get_node(50) # class method defined
        
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
        """
        @test Error searching for non-existent node.
        @brief get_node() must raise KeyError if the ID does not exist.
        """
        model = FEMModel()
        model.nodes[1] = Node(nid=1, coord_id=0, x1=0., x2=0., x3=0., cd_id=0, ps="", seid=0)
        id_searched = 999

        with pytest.raises(KeyError) as excinfo:
            model.get_node(id_searched)
        error_message = str(excinfo.value)
        assert "[FEMModel ERROR]" in error_message
        assert str(id_searched) in error_message
        assert "does not exist" in error_message

    # --- TEST DE OTRAS ENTIDADES (Materiales, Elementos, Cargas) ---

    def test_store_materials(self):
        """
        @test Material Storage (Dictionary).
        @brief Check insertion and retrieval by material ID.
        """
        model = FEMModel()
        
        mat = IsotropicMaterial(mid=100, E=210e9, nu=0.3)
        model.materials[100] = mat

        assert 100 in model.materials
        assert model.materials[100].E == 210e9
        assert len(model.materials) == 1

    def test_store_elements(self):
        """
        @test Elements Storage (Dictionary).
        @brief Check insertion and retrieval of elements fields.
        """
        model = FEMModel()
        
        elem = RodElement(eid=10, prop_id=5, nodes=[1, 2])
        model.elements[10] = elem
        
        assert 10 in model.elements
        assert model.elements[10].nodes == [1, 2]
        assert model.elements[10].prop_id == 5
    
    def test_store_loads(self):
        """
        @test Loads Storage (List).
        @brief Check insertion and retrieval of loads fields and that the loads are added sequentially to the list.
        """
        model = FEMModel()

        load1 = SinglePointLoad(sid=2, appnode_id=5, module=100.0, vector=[0,0,1])
        load2 = SinglePointLoad(sid=2, appnode_id=6, module=200.0, vector=[0,1,0])
        
        model.loads.append(load1)
        model.loads.append(load2)

        assert len(model.loads) == 2
        assert model.loads[0].module == 100.0
        assert model.loads[1].module == 200.0
        assert model.loads[0].sid == 2 # Same Set ID allowed
        assert model.loads[1].sid == 2 # Same Set ID allowed

    def test_store_constraints(self):
        """
        @test Constraints Storage (List).
        @brief Verify that constraints are added correctly, allowing multiple entries.
        """
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
    """
    @class TestBarProperty
    @brief Specific tests for bar properties (PBAR).
    @details Checks correct initialization and physical validations (inertia).
    """

    def test_full_initialization(self):
        """
        @test BarProperty inicialization
        @brief Checks that the beam property correctly stores all coefficients
        extended (I12, K1, K2) and the stress recovery points.
        """
        stress_points = [(0.1, 0.1), (0.1, -0.1), (-0.1, -0.1), (-0.1, 0.1)]
        
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
            stress_points=stress_points
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
        """
        @test Checks physics consistency (Mohr's circle).
        @brief Verify that a ValueError is raised if I12^2 > I1*I2.
        """
        # Case: I1*I2 = 4. If I12=10 -> 100 > 4 -> Physically impossible.
        with pytest.raises(ValueError) as excinfo:
            BarProperty(
                pid=99, mat_id=1, A=1.0, J=1.0,
                I1=2.0, I2=2.0, 
                I12=10.0 
            )
        assert "Inconsistent inertial properties" in str(excinfo.value)


class TestBarElement:
    """
    @class TestBarElement
    @brief Tests suite for the beam element entity (CBAR).
    """
    def test_full_initialization(self):
        """
        @test BarProperty inicialization
        @brief Checks the initialization of a CBAR element with advanced options
        (Offset, Orientation, Pin Flags).
        """
        elem = BarElement(
            eid=200,
            prop_id=10,
            nodes=[1, 2],
            orientation_vector=[0.0, 1.0, 0.0],
            g0=None,
            offt='GGG',
            pa=45,          # Pin Flag in node A (free Rx, Ry)
            pb=0,
            wa=[0.0, 0.05, 0.0], # Offset in A
            wb=[0.0, 0.05, 0.0]  # Offset in B
        )

        assert elem.eid == 200
        assert elem.pa == 45
        assert elem.pb == 0
        assert elem.wa[1] == 0.05
        assert elem.orientation_vector == [0.0, 1.0, 0.0]


class TestShellProperty:
    """
    @class TestShellProperty
    @brief Test suite for ShellProperty (PSHELL).
    """
    def test_initialization(self):
        prop = ShellProperty(pid=1, mat_id=10, t=0.005, mat_id_bending=10)
        assert prop.pid == 1
        assert prop.t == 0.005
        assert prop.mat_id == 10
        assert prop.mat_id_bending == 10


class TestCTRIA3Element:
    """
    @class TestCTRIA3Element
    @brief Test suite for CTRIA3Element.
    """
    def test_initialization(self):
        elem = CTRIA3Element(eid=100, prop_id=1, nodes=[1, 2, 3], theta=30.0)
        assert elem.eid == 100
        assert elem.nodes == [1, 2, 3]
        assert elem.prop_id == 1
        assert elem.theta == 30.0





