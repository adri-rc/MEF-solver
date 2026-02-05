"""
@file test_nastran_reader_end_to_end.py
@brief End-to-End validation tests for the Nastran parser.
@details
This script resides in 'tests/preprocessing_tst/' but consumes files 
located in 'tests/solver_tst/validation_cases/' to avoid data duplication.
"""

import pytest
import numpy as np
from pathlib import Path
from src.preprocessing.nastran_reader import parse_nastran_input
from src.preprocessing.data_model import (
    FEMModel, 
    SinglePointLoad, 
    SinglePointMoment, 
    SPCConstraint,
    RodElement,       
    IsotropicMaterial 
)

# DIRECTORY FIXTURE
@pytest.fixture
def shared_validation_dir():
    """
@brief Resolves the path to 'tests/solver_tst/validation_cases' 
           relative to this script's location.
    @details
    Assumed structure:
      root/
        tests/
          preprocessing_tst/  <-- (__file__)
          solver_tst/
            validation_cases/ <-- goal
    """
    current_dir = Path(__file__).parent
    target_dir = current_dir.parent / "solver_tst" / "validation_cases"
    
    if not target_dir.exists():
        pytest.fail(f"Shared validation directory not found at: {target_dir}")
        
    return target_dir

# ==============================================================================
# VALIDATION CASE DATABASE (PARSER)
# ==============================================================================

READER_VALIDATION_CASES = [
    # --- CASO 8: Riel 3D (Validación Básica) ---
    {
        "case_name": "Ej8_Rail_Y_Axis",
        "folder": "case_8",
        "filename": "Ej8.dat",
        "expected_counts": {
            "nodes": 2,
            "elements": 1,
            "materials": 1,
            "properties": 1,
            "loads": 1,
            "constraints": 2
        },
        "checks": {
            "nodes": {
                7: [40.0, 50.0, 30.0],
                8: [10.0, 80.0, 50.0]
            },
            "elements": {
                # 104 is CROD -> RodElement
                104: {"nodes": [7, 8], "pid": 401, "expected_class": RodElement}
            },
            "materials": {
                # 501 is MAT1 -> IsotropicMaterial
                501: {"E": 210000.0, "nu": 0.3, "expected_class": IsotropicMaterial}
            },
            "loads": [
                {"type": "FORCE", "sid": 704, "node": 8, "mag": 11.72, "vec": [-3., 3., 2.]}
            ],
            "constraints": [
                {"sid": 204, "node": 7, "comps": "123456"},
                {"sid": 204, "node": 8, "comps": "13456"}
            ]
        }
    },

    # --- CASO 10: Cercha DAA (Validación de INCLUDEs y múltiples elementos) ---
    {
        "case_name": "Ej10_DAA_Class_Truss",
        "folder": "case_10",
        "filename": "Ej_DAA_clase.dat",
        "expected_counts": {
            "nodes": 6,
            "elements": 9,
            "materials": 1,
            "properties": 7,
            "loads": 2,
            "constraints_min": 6
        },
        "checks": {
            "nodes": {
                1: [0.0, 0.0, 0.0],
                3: [4.0, 0.0, 0.0],
                6: [12.0, 0.0, 0.0]
            },
            "elements": {
                1: {"nodes": [1, 2], "pid": 11, "expected_class": RodElement},
                9: {"nodes": [4, 6], "pid": 11, "expected_class": RodElement}
            },
            "materials": {
                1: {"E": 2.0e11, "nu": 0.3, "expected_class": IsotropicMaterial}
            },
            "loads": [
                {"type": "FORCE", "sid": 1, "node": 2, "mag": 40000.0, "vec": [0., -1., 0.]},
                {"type": "FORCE", "sid": 1, "node": 5, "mag": 100000.0, "vec": [0., -1., 0.]}
            ]
        }
    },

    # --- CASO 12: Cadena 3D (Validación de MOMENT y Enforced Displacements) ---
    {
        "case_name": "Ej12_3D_Chain_Enforced",
        "folder": "case_12",
        "filename": "Ej12.dat",
        "expected_counts": {
            "nodes": 4,
            "elements": 3,
            "loads": 2,
            "constraints_min": 4
        },
        "checks": {
            "nodes": {
                2: [10., 10., 5.],
                4: [30., 10., 5.]
            },
            "elements": {
                102: {"nodes": [2, 3], "pid": 402, "expected_class": RodElement}
            },
            "materials": {
                501: {"E": 2.1e5, "expected_class": IsotropicMaterial},
                502: {"E": 7.0e4, "expected_class": IsotropicMaterial}
            },
            "loads": [
                {"type": "FORCE", "sid": 1100, "node": 2, "mag": 1000.0, "vec": [0., -1., 0.]},
                {"type": "MOMENT", "sid": 1100, "node": 3, "mag": 1000.0, "vec": [1., 0., 0.]}
            ],
            "constraints": [
                {"sid": 1100, "node": 1, "comps": "123456", "value": 0.0},
                {"sid": 1100, "node": 4, "comps": "2", "value": 2.0},
                {"sid": 1100, "node": 4, "comps": "3", "value": -5.0}
            ]
        }
    }
]

# ==============================================================================
# TEST CLASS
# ==============================================================================

class TestNastranReaderValidation:
    """
    @class TestNastranReaderValidation
    @brief Parameterized tests to verify correct Nastran file parsing.
    """

    @pytest.mark.parametrize("case_data", READER_VALIDATION_CASES, ids=lambda x: x["case_name"])
    def test_parse_full_model(self, case_data, shared_validation_dir):
        """
        @brief Loads a real file and exhaustively verifies the FEMModel content.
        @param shared_validation_dir Fixture pointing to 'tests/solver_tst/validation_cases'
        """
        # 1. Path
        folder = case_data.get("folder", "")
        filename = case_data["filename"]
        file_path = shared_validation_dir / folder / filename

        if not file_path.exists():
            pytest.fail(f"Archivo de test no encontrado en: {file_path}")

        # 2. Parsewr execution (nastran_reader)
        try:
            model = parse_nastran_input(str(file_path))
        except Exception as e:
            pytest.fail(f"Parser failed reading{filename}: {e}")

        # 3. Validate Counts (Sanity Check)
        counts = case_data.get("expected_counts", {})
        if "nodes" in counts:
            assert len(model.nodes) == counts["nodes"], \
                f"Node count mismatch: {len(model.nodes)} != {counts['nodes']}"
        if "elements" in counts:
            assert len(model.elements) == counts["elements"], \
                f"Element count mismatch: {len(model.elements)} != {counts['elements']}"
        if "loads" in counts:
            assert len(model.loads) == counts["loads"], \
                f"Load count mismatch: {len(model.loads)} != {counts['loads']}"

        # 4. Validate Specific Data
        checks = case_data.get("checks", {})

        # --- A) NODES ---
        for nid, expected_coords in checks.get("nodes", {}).items():
            assert nid in model.nodes, f"Node {nid} missing."
            node_obj = model.nodes[nid]
            np.testing.assert_allclose(node_obj.get_coords(), expected_coords, atol=1e-5, 
                err_msg=f"Incorrect coords for Node {nid}")

        # --- B) ELEMENTS (Type and property validation) ---
        for eid, data in checks.get("elements", {}).items():
            assert eid in model.elements, f"Element {eid} missing."
            elem = model.elements[eid]
            
            if "nodes" in data:
                assert elem.nodes == data["nodes"]
            if "pid" in data:
                assert elem.prop_id == data["pid"]
            
            # Strict Python class verification (RodElement vs BarElement)
            if "expected_class" in data:
                assert isinstance(elem, data["expected_class"]), \
                    f"Elem {eid} is {type(elem).__name__}, expected {data['expected_class'].__name__}"

        # --- C) MATERIALS (Type and physics validation) ---
        for mid, props in checks.get("materials", {}).items():
            assert mid in model.materials, f"Material {mid} missing."
            mat = model.materials[mid]
            
            if "E" in props:
                assert mat.E == pytest.approx(props["E"])
            if "nu" in props:
                assert mat.nu == pytest.approx(props["nu"])
                
            if "expected_class" in props:
                assert isinstance(mat, props["expected_class"]), \
                    f"Mat {mid} is {type(mat).__name__}, expected {props['expected_class'].__name__}"

        # --- D) LOADS ---
        for load_check in checks.get("loads", []):
            found = False
            # Distinguish between Force Class and Moment Class
            target_type = SinglePointLoad if load_check["type"] == "FORCE" else SinglePointMoment
            
            for load in model.loads:
                if (isinstance(load, target_type) and 
                    load.sid == load_check["sid"] and 
                    load.appnode_id == load_check["node"]):
                    
                    found = True
                    # Magnitude
                    assert load.module == pytest.approx(load_check["mag"])
                    # Vector (assuming parser keeps raw direction vector)
                    np.testing.assert_allclose(load.vector, load_check["vec"], atol=1e-4)
                    break
            
            assert found, f"Expected load not found: {load_check}"

        # --- E) CONSTRAINTS (SPCs) ---
        for constr_check in checks.get("constraints", []):
            found = False
            for spc in model.constraints:
                # Check SID, Node and Components
                if (spc.sid == constr_check["sid"] and 
                    spc.node_id == constr_check["node"] and
                    spc.components == constr_check["comps"]):
                    
                    # If SPC has enforced value, verify it
                    if "value" in constr_check:
                        assert spc.value == pytest.approx(constr_check["value"]), \
                            f"Incorrect SPC value at Node {spc.node_id}"
                    
                    found = True
                    break
            
            assert found, f"Expected constraint not found: {constr_check}"