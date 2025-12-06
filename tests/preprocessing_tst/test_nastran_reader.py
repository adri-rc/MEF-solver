import pytest
from preprocessing.data_model import FEMModel
from tests.conftest import test_data_dir
from preprocessing.nastran_reader import (
    parse_nastran_input, 
    _parse_nastran_number, 
    _get_entry_fields
)

"""
@file test_nastran_reader.py
@brief Tests for the Nastran reading module.
@details Three levels:
    1. Parsing Logic: Unit tests for internal functions (_parse_number).
    2. Card Handling: Integration tests generating temporary files in memory.
    3. File Loading: Tests with real external files.
"""


# region --- Parsing Logic Tests --

class TestNastranParsingLogic:
    """
    @class TestNastranParsingLogic
    @brief Checks the robustness of low-level functions (text parsing).
    """

    @pytest.mark.parametrize("input_str, expected", [
        ("1.0", 1.0),
        ("  2.0  ", 2.0),       # blank spaces befora and behind
        ("2.0  ", 2.0),         # blank spaces behind
        ("  2.0", 2.0),         # blank spaces before
        ("1.E-4", 1e-4),        # Standard notation
        ("1.-4", 1e-4),         # Compressed Nastran notation
        ("1.2+3", 1200.0),      # Nastran notation with +
        ("2.D+2", 200.0),       # Double Precision Nastran Notation (D)
        ("-1.-1", -0.1),        # Initial negative sign and negative exponent
        ("", None),             # Empty field
        ("     ", None)         # Field full of blank spaces
    ])
    def test_parse_nastran_number_valid(self, input_str, expected):
        """
        @test Validation of number editing from Nastran fields.
        @brief Verify that _parse_nastran_number handles all possible Nastran formats.
        """
        result = _parse_nastran_number(input_str)
        if expected is None:
            assert result is None
        else:
            assert result == pytest.approx(expected)
    
    @pytest.mark.parametrize("invalid_case", [
        "TEXTO_INVALIDO",
        "1.2.2",   # Dos puntos decimales
        "1.P-7",   # Caracter inválido en exponente
        "1,856"    # Coma decimal (Python usa punto)
    ])
    def test_parse_nastran_number_invalid(self, invalid_case):
        """
        @test Invalid number error.
        @brief Checks that it throws a ValueError with garbage text.
        """
        with pytest.raises(ValueError) as excinfo:
            _parse_nastran_number(invalid_case)
        expected_msg = f"[ERROR] The value '{invalid_case}' can not be converted"
        assert expected_msg in str(excinfo.value)
        

    def test_get_entry_fields_free_vs_fixed(self):
        """
        @test Field separation (Free vs Fixed).
        @brief Checks that _get_entry_fields distinguishes between commas and fixed width.
        """
        
        line_free = "GRID, 1, , 1.0, 2.0" # Caso 1: Free Field (CSV)
        fields_free = _get_entry_fields(line_free)
        assert fields_free == ["GRID", "1", "", "1.0", "2.0"]

        # Caso 2: Fixed Field ($1234567)
        #            "GRID    " (8) "1       " (8) "        " (8) "1.0     " (8)
        #             $1234567$1234567$1234567$1234567$1234567$1234567$1234567
        line_fixed = "GRID           1        1.0       2.0   randoM   1 2 3  "
        fields_fixed = _get_entry_fields(line_fixed)
        assert fields_fixed[0] == "GRID"
        assert fields_fixed[1] == "1"
        assert fields_fixed[2] == ""
        assert fields_fixed[3] == "1.0"
        assert fields_fixed[4] == "2.0"
        assert fields_fixed[5] == "randoM"
        assert fields_fixed[6] == "1 2 3"

# endregion --- Parsing Logic Tests --


# region --- Dummy Nastran texts Tests --
class TestNastranIntegration:
    """
    @class TestNastranIntegration
    @brief Tests that create temporary files to test complete text reading.
    """

    def test_read_complete_model_from_temp_file(self, tmp_path):
        """
        @test Complete reading of a dummy.bdf
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
        d = tmp_path / "subdir"
        d.mkdir()
        p = d / "temp_model.bdf"
        p.write_text(bdf_content, encoding='utf-8')

        model = parse_nastran_input(str(p)) #run function to test

        # Nodes check
        assert 101 in model.nodes
        assert 102 in model.nodes
        assert model.nodes[101].x == 1.0
        assert model.nodes[101].z == 0.0
        assert model.nodes[102].y == 3.
        
        # Materials check
        assert 5 in model.materials
        assert model.materials[5].E == 2.1e11
        assert model.materials[5].G is not None # G is autocalculated

        # Elements check
        assert 20 in model.elements
        assert model.elements[20].nodes == [101, 102]

        # Constraints check
        assert len(model.constraints) == 2
        assert model.constraints[0].node_id == 101
        assert model.constraints[1].components == "123"

        # Loads check
        assert len(model.loads) == 1
        assert model.loads[0].module == 500.0
        assert model.loads[0].vector[1] == 1.

    def test_include_recursive_file(self, tmp_path):
        """
        @test Reading INCLUDE cards(Just one include)
        @brief Checks that the reader follows the referenced files.
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
        assert model.nodes[99].y == 9.9
        assert model.materials[5].G is not None

        assert 101 in model.nodes
        assert model.nodes[102].y == 3.

    def test_include_recursive_file(self, tmp_path):
        """
        @test Reading INCLUDE cards(Just one include)
        @brief Checks that the reader follows the referenced files.
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
INCLUDE '{include_file1.name}
INCLUDE '{include_file2.name}'"""
        main_file = tmp_path / "master.bdf"
        main_file.write_text(main_content)

        model = parse_nastran_input(str(main_file))

        assert 99 in model.nodes
        assert model.nodes[99].y == 9.9
        assert model.materials[5].G is not None
        assert 88 in model.nodes
        assert model.nodes[88].y == 0.
        assert 8 in model.materials
        assert model.materials[8].G is not None
        assert 101 in model.nodes
        assert model.nodes[102].y == 3.
        
# endregion --- Dummy Nastran texts Tests --

# region --- Nastran Real Input Tests --
class TestRealNastranFiles:
    """
    @class TestRealNastranFiles
    @brief Tests that load existing .dat and .bdf files in the 'files4test' folder.
    """

    def test_load_reference_bdf(self, test_data_dir):
        """
        @test Carga de 'test_model.bdf'.
        @brief Verifica que el archivo de prueba estándar se lea correctamente.
        @param test_data_dir Fixture definida en conftest.py que apunta a files4test.
        """
    
        dat_path = test_data_dir / "Ej3.dat"
        if not dat_path.exists():
            pytest.skip(f"The file {dat_path} does not exist. Create the file to run this test.")
        
        model = parse_nastran_input(str(dat_path)) # main nastran_reader function
        
        assert isinstance(model, FEMModel)
        assert len(model.nodes) > 0, "Model should not be empty."
        assert len(model.elements) > 0, "There should be elements"
        assert model.loads[0].module == 500.
        assert model.nodes[4].y == 57.735
        assert 502 in model.materials
        assert model.properties[402].A == 225.


    def test_file_not_found(self):
        """
        @test File not found error.
        @brief Checks that parse_nastran_input throws the correct error if the path is incorrect.
        """
        with pytest.raises(Exception) as excinfo:
            parse_nastran_input("ruta/inventada/no_existe.dat")

        assert "file not found" in str(excinfo.value)

# endregion --- Nastran Real Input Tests --