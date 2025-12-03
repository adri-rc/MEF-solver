import pytest
import sys
import os
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(os.path.dirname(current_dir))
target_dir = os.path.join(root_dir, 'src', 'preprocessing')

if not os.path.exists(target_dir):
    print(f"[ERROR] Path not found: {target_dir}")

sys.path.insert(0, target_dir) # Insert into sys.path

from data_model import (          
    Node,
    IsotropicMaterial, # Material
    FEMModel,
    Property,
    RodElement,        # Element
    SPCConstraint,     # Constraint
    SinglePointLoad,   # LoadCase
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

    def test_consistency_warning(self, capsys: pytest.CaptureFixture[str]):
        """
        @test Physical consistency warning.
        @brief If the user provides E, G, and nu, but they do not match the formula,
                a warning is raised.
        """
        # Non real material (E=200 y nu=0.0-> Should have G=100.) but G is forced to 50. (50% Error)
        IsotropicMaterial(mid=20, E=200.0, nu=0.0, G=50.0)
        
        captured = capsys.readouterr()
        
        assert "Review inputs in material id: '20'" in captured.out
        assert "consistent" in captured.out
        assert "E theoretical" in captured.out

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










