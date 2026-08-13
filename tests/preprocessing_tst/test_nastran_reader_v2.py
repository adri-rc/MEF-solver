import pytest
from src.preprocessing.data_model import FEMModel
from tests.conftest import test_data_dir
from src.preprocessing.nastran_reader_v2 import parse_nastran_input

"""
@file test_nastran_reader_v2.py
@brief Tests for the Nastran reader module using pyNastran (v2).
@details Two levels (Parsing logic is now handled by pyNastran internally):
    1. Card Handling: Integration tests generating temporary files in memory.
    2. File Loading: Tests with real external files.
    3. Error Handling: 
        - Ensuring invalid materials stop execution.
        - Ensuring invalid elements stop execution.
        - Ensuring invalid spcs stop execution.
"""

# region --- Dummy Nastran texts Tests --
class TestNastranIntegrationV2:
    """
    @class TestNastranIntegrationV2
    @brief Tests that create temporary files to test the pyNastran integration.
    """

    def test_read_complete_model_from_temp_file(self, tmp_path):
        """
        @test Complete reading of a dummy.bdf with pyNastran.
        @brief Use the pytest 'tmp_path' fixture to avoid cluttering up the disk.
        """
        
        bdf_content = """
$ Model Test
GRID, 101, , 1.0, 0.0, 0.
GRID, 102, , 2.0, 3.000, 0.
MAT1, 5, 2.1E11, , 0.3

PROD, 10, 5, 0.1
CROD, 20, 10, 101, 102
SPC1, 100, 123456, 101
SPC1, 100, 123, 102
FORCE, 200, 102, , 500.0, 0.0, 1.0, 0.0
        """
        
        # Temp file
        d = tmp_path / "subdir_v2"
        d.mkdir()
        p = d / "temp_model_v2.bdf"
        p.write_text(bdf_content, encoding='utf-8')

        # Run function to test (V2)
        model = parse_nastran_input(str(p)) 

        # Nodes check
        # Note: pyNastran resolves coordinates globally, 
        # so we expect correct values in model.nodes
        assert 101 in model.nodes
        assert 102 in model.nodes
        assert model.nodes[101].x == pytest.approx(1.0)
        assert model.nodes[101].z == pytest.approx(0.0)
        assert model.nodes[102].y == pytest.approx(3.0)
        
        # Materials check
        assert 5 in model.materials
        assert model.materials[5].E == pytest.approx(2.1e11)
        assert model.materials[5].G is not None # G is autocalculated by DataModel

        # Elements check
        assert 20 in model.elements
        assert model.elements[20].nodes == [101, 102]

        # Constraints check
        # pyNastran might read SPC1 order differently, so we check existence
        assert len(model.constraints) == 2
        spc_nodes = [c.node_id for c in model.constraints]
        assert 101 in spc_nodes
        assert 102 in spc_nodes

        # Loads check
        assert len(model.loads) == 1
        assert model.loads[0].module == pytest.approx(500.0)
        assert model.loads[0].vector[1] == pytest.approx(1.0)

    def test_single_include(self, tmp_path):
        """
        @test Reading INCLUDE cards (Single include)
        @brief Checks that pyNastran follows the referenced files correctly.
        """
        # Include_file.bdf
        include_content = "GRID, 99, , 9.9, 9.9, 9.9\nMAT1, 5, 2.1E11,, 0.3"
        include_file = tmp_path / "mesh.bdf"
        include_file.write_text(include_content)

        # main.dat
        parent_content = f"""
GRID, 101, , 1.0, 0.0, 0.
GRID, 102, , 2.0, 3.0, 0.
INCLUDE '{include_file.name}'"""
        parent_file = tmp_path / "master.bdf"
        parent_file.write_text(parent_content)

        model = parse_nastran_input(str(parent_file))

        assert 99 in model.nodes
        assert model.nodes[99].y == pytest.approx(9.9)
        assert 5 in model.materials
        assert model.materials[5].G is not None

        assert 101 in model.nodes
        assert model.nodes[102].y == pytest.approx(3.0)

    def test_multiple_includes(self, tmp_path):
        """
        @test Reading multiple INCLUDE cards
        @brief Checks that the reader follows multiple referenced files.
        """
        # Include_file1.bdf
        include_content1 = "GRID, 99, , 9.9, 9.9, 9.9\nMAT1, 5, 2.1E11,, 0.3"
        include_file1 = tmp_path / "mesh1.bdf"
        include_file1.write_text(include_content1)

        # Include_file2.bdf
        include_content2 = "GRID, 88, , 8., 0., 1.\nMAT1, 8, 2.1E11,, 0.3"
        include_file2 = tmp_path / "mesh2.bdf"
        include_file2.write_text(include_content2)


        # main.dat
        main_content = f"""
GRID, 101, , 1.0, 0.0, 0.
GRID, 102, , 2.0, 3.0, 0.
INCLUDE '{include_file1.name}'
INCLUDE '{include_file2.name}'"""
        main_file = tmp_path / "master.bdf"
        main_file.write_text(main_content)

        model = parse_nastran_input(str(main_file))

        # Check nodes from include 1
        assert 99 in model.nodes
        assert model.nodes[99].y == pytest.approx(9.9)
        
        # Check nodes from include 2
        assert 88 in model.nodes
        assert model.nodes[88].y == pytest.approx(0.)
        
        # Check master nodes
        assert 101 in model.nodes
    
    def test_pbar_stress_points_parsing(self, tmp_path):
        """
        @test PBAR Stress Points (C, D, E, F).
        @brief Verifies that stress recovery points are correctly mapped to BarProperty.
        """
        # Definition:
        # PBAR, PID=1, MID=1, A=10., I1=20., I2=20., J=10.
        # Continuation line (+): C(1.,1.), D(1.,-1.), E(-1.,-1.), F(-1.,1.)
        # Format PBAR fields 9-16 correspond to C1,C2, D1,D2...
        
        bdf_content = """
MAT1, 1, 2.1E11, , 0.3
PBAR, 100, 1, 10.0, 20.0, 20.0, 10.0, 0.0
, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0
        """
        p = tmp_path / "pbar_stress.bdf"
        p.write_text(bdf_content)

        model = parse_nastran_input(str(p))
        
        prop = model.properties[100]
        
        # Check that we have 4 points
        assert len(prop.stress_points) == 4
        # Check values (Point C)
        assert prop.stress_points[0] == pytest.approx((1.0, 1.0))
        # Check values (Point E)
        assert prop.stress_points[2] == pytest.approx((-1.0, -1.0))
      
    def test_ctria3_parsing(self, tmp_path):
        """
        @test CTRIA3 and PSHELL parsing.
        @brief Verifies that shell elements and properties are correctly mapped.
        """
        bdf_content = """
GRID, 1, , 0.0, 0.0, 0.0
GRID, 2, , 1.0, 0.0, 0.0
GRID, 3, , 0.0, 1.0, 0.0
MAT1, 1, 70.E9, , 0.3
PSHELL, 10, 1, 0.005, 1, , 1 
CTRIA3, 100, 10, 1, 2, 3, 30.0
        """
        p = tmp_path / "ctria3_test.bdf"
        p.write_text(bdf_content)

        model = parse_nastran_input(str(p))

        # Check Property
        assert 10 in model.properties
        prop = model.properties[10]
        assert prop.t == 0.005
        assert prop.mat_id == 1
        assert prop.mat_id_bending == 1 # mid2 mapped


        # Check Element
        assert 100 in model.elements
        elem = model.elements[100]
        assert elem.nodes == [1, 2, 3]
        assert elem.prop_id == 10
        assert elem.theta == 30.0

    def test_cquad4_parsing(self, tmp_path):
        """
        @test CQUAD4 and PSHELL parsing.
        @brief Verifies that quadrilateral shell elements and properties are correctly mapped.
        """
        bdf_content = """
GRID, 1, , 0.0, 0.0, 0.0
GRID, 2, , 1.0, 0.0, 0.0
GRID, 3, , 1.0, 1.0, 0.0
GRID, 4, , 0.0, 1.0, 0.0
MAT1, 1, 70.E9, , 0.3
PSHELL, 10, 1, 0.01, 1, , 1 
CQUAD4, 200, 10, 1, 2, 3, 4, 45.0
        """
        p = tmp_path / "cquad4_test.bdf"
        p.write_text(bdf_content)

        model = parse_nastran_input(str(p))

        # Check Property
        assert 10 in model.properties
        prop = model.properties[10]
        assert prop.t == 0.01
        assert prop.mat_id == 1
        assert prop.mat_id_bending == 1

        # Check Element
        assert 200 in model.elements
        elem = model.elements[200]
        assert elem.nodes == [1, 2, 3, 4]
        assert elem.prop_id == 10
        assert elem.theta == 45.0
        assert elem.__class__.__name__ == "CQUAD4Element"

    def test_pshell_blank_mid3_does_not_inherit_mid2(self, tmp_path):
        """A blank MID3 means no physical transverse-shear material."""
        bdf_content = """
GRID, 1, , 0.0, 0.0, 0.0
GRID, 2, , 1.0, 0.0, 0.0
GRID, 3, , 1.0, 1.0, 0.0
GRID, 4, , 0.0, 1.0, 0.0
MAT1, 1, 70.E9, , 0.3
PSHELL, 10, 1, 0.01, 1, 0.75, , 0.5
CQUAD4, 200, 10, 1, 2, 3, 4
        """
        path = tmp_path / "pshell_blank_mid3.bdf"
        path.write_text(bdf_content)

        prop = parse_nastran_input(str(path)).properties[10]

        assert prop.mat_id_bending == 1
        assert prop.mat_id_shear is None
        assert prop.bending_inertia_ratio == pytest.approx(0.75)
        assert prop.shear_thickness_ratio == pytest.approx(0.5)

# endregion --- Dummy Nastran texts Tests --


class TestNastranEntityMapping:
    """
    @brief Tests to verify correct extraction of physical entities.
    """

    @pytest.fixture
    def valid_bdf_path(self, tmp_path):
        """Creates a comprehensive BDF file for positive testing."""
        content = """
$ --- NODES ---
GRID, 1, , 0.0, 0.0, 0.0
GRID, 2, , 1.0, 0.0, 0.0
GRID, 3, , 2.0, 0.0, 0.0
$ --- MATERIALS (E=2.1e11, nu=0.3 -> G approx 8e10) ---
MAT1, 10, 2.1E11, , 0.3
$ --- PROPERTIES ---
PROD, 100, 10, 0.01
$ --- ELEMENTS ---
CROD, 50, 100, 1, 2
CROD, 51, 100, 2, 3
$ --- LOADS ---
$ Force on Node 2, Mag 1000, Dir Y (0,1,0)
FORCE, 10, 2, , 1000., 0., 1., 0.
$ Moment on Node 3, Mag 50, Dir Z (0,0,1)
MOMENT, 10, 3, , 50., 0., 0., 1.
$ --- CONSTRAINTS ---
$ SPC1: Node 1 fixed in 123456
SPC1, 20, 123456, 1
$ SPC: Node 3 fixed in 3 (Z) with displacement 0.005
SPC, 20, 3, 3, 0.005
        """
        p = tmp_path / "valid_model.bdf"
        p.write_text(content, encoding='utf-8')
        return str(p)

    def test_entity_counts(self, valid_bdf_path):
        """
        @test Verify the number of extracted entities.
        """
        model = parse_nastran_input(valid_bdf_path)

        assert len(model.nodes) == 3
        assert len(model.elements) == 2
        assert len(model.materials) == 1
        assert len(model.properties) == 1
        assert len(model.loads) == 2      # 1 Force + 1 Moment
        assert len(model.constraints) == 2 # 1 SPC1 + 1 SPC

    def test_loads_extraction(self, valid_bdf_path):
        """
        @test Verify Force and Moment data mapping.
        """
        model = parse_nastran_input(valid_bdf_path)
        
        # Filter loads by type (implicitly by checking attributes or class)
        # Assuming mapped order matches file order for this simple case
        force_load = model.loads[0] 
        moment_load = model.loads[1]

        # Check Force
        assert force_load.appnode_id == 2
        assert force_load.module == pytest.approx(1000.0)
        assert force_load.vector == [0.0, 1.0, 0.0]

        # Check Moment
        assert moment_load.appnode_id == 3
        assert moment_load.module == pytest.approx(50.0)
        assert moment_load.vector == [0.0, 0.0, 1.0]

    def test_constraints_extraction(self, valid_bdf_path):
        """
        @test Verify SPC (enforced disp) and SPC1 (0 disp) mapping.
        """
        model = parse_nastran_input(valid_bdf_path)
        
        # Check SPC1 (Node 1, comps 123456, value 0)
        spc1 = next(c for c in model.constraints if c.node_id == 1)
        assert "1" in spc1.components and "6" in spc1.components
        assert spc1.value == 0.0

        # Check SPC (Node 3, comp 3, value 0.005)
        spc_custom = next(c for c in model.constraints if c.node_id == 3)
        assert spc_custom.components == "3"
        assert spc_custom.value == pytest.approx(0.005)


class TestNastranErrorHandling:
    """
    @brief Tests for failure scenarios and validations.
    """

    def test_material_inconsistency_failure(self, tmp_path):
        """
        @test Fail Fast on Inconsistent Material.
        @brief Input E, G, nu that violate E = 2G(1+nu). Should raise ValueError.
        """
        # --- CORREGIDO: Sin 'BEGIN BULK' ni 'ENDDATA' para evitar conflicto con punch=True ---
        bad_mat_content = """
MAT1, 1, 100.0, 10.0, 0.0
        """
        p = tmp_path / "bad_material.bdf"
        p.write_text(bad_mat_content)

        with pytest.raises(ValueError) as excinfo:
            parse_nastran_input(str(p))
        
        error_msg = str(excinfo.value)
        assert "physically inconsistent" in error_msg or "Review material inputs" in error_msg

    def test_unsupported_element_warning(self, tmp_path, capsys):
        """
        @test Unsupported Element Detection (Fail Fast).
        @brief If the BDF contains elements we don't map (e.g. CTRIA3), 
               the reader MUST raise a ValueError.
        """
        # --- CORREGIDO: Añadimos MAT1 para que PSHELL sea válido y no falle por CrossRef ---
        content = """
GRID, 1, , 0.,0.,0.
GRID, 2, , 1.,0.,0.
GRID, 3, , 0.,1.,0.
GRID, 4, , 1.,1.,0.
MAT1, 1, 2.1E11, , 0.3
PSOLID, 10, 1
CTETRA, 99, 10, 1, 2, 3, 4

        """
        p = tmp_path / "unsupported_elem.bdf"
        p.write_text(content)

        with pytest.raises(ValueError) as excinfo:
            parse_nastran_input(str(p))
        
        # Verify message is about the element type
        error_msg = str(excinfo.value)
        assert "CTETRA" in error_msg or "FATAL ERROR" in error_msg


    def test_unsupported_load_error(self, tmp_path):
        """
        @test Unsupported Load Detection.
        @brief If the BDF contains loads like GRAV, it must fail.
        """
        content = """
GRID, 1, , 0.,0.,0.
MAT1, 1, 2.1E11, , 0.3
GRAV, 10, 0, 9.81, 0., -1., 0.
        """
        p = tmp_path / "unsupported_load.bdf"
        p.write_text(content)

        with pytest.raises(ValueError) as excinfo:
            parse_nastran_input(str(p))
        
        assert "GRAV" in str(excinfo.value)

    def test_unsupported_mpc_error(self, tmp_path):
        """
        @test Unsupported MPC/RBE Detection.
        @brief If the BDF contains Rigid Elements (RBE2), it must fail.
        """
        content = """
GRID, 1, , 0.,0.,0.
GRID, 2, , 1.,0.,0.
GRID, 3, , 2.,0.,0.
RBE2, 99, 1, 123456, 2, 3
        """
        p = tmp_path / "unsupported_mpc.bdf"
        p.write_text(content)

        with pytest.raises(ValueError) as excinfo:
            parse_nastran_input(str(p))
        
        error_str = str(excinfo.value)
        # Verify message mentions MPC or Rigid
        assert "MPC" in error_str or "RBE" in error_str or "Rigid" in error_str
