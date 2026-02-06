import pytest
import numpy as np
from pathlib import Path

from src.preprocessing.nastran_reader_v2 import parse_nastran_input

from src.preprocessing.data_model import (
    FEMModel, 
    SinglePointLoad, 
    SinglePointMoment, 
    SPCConstraint,
    RodElement,
    BarElement,      
    IsotropicMaterial 
)

"""
@file test_nastran_reader_v2_end_to_end.py
@brief End-to-End validation tests for the new pyNastran-based parser (V2).
@details
    This script ensures that the new V2 reader produces the exact same 
    physical model (FEMModel) as the legacy reader for known validation cases.
    It consumes files located in 'tests/solver_tst/validation_cases/'.
"""

# DIRECTORY FIXTURE
@pytest.fixture
def shared_validation_dir():
    """
    @brief Resolves the path to 'tests/solver_tst/validation_cases' 
           relative to this script's location.
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
    },
    # --- CASO 14: CBAR 3D (Flexión Biaxial y Torsión) ---
    {
        "case_name": "Ej14_CBAR_3D",
        "folder": "case_14",
        "filename": "Ej14.dat",
        "expected_counts": {
            "nodes": 2, 
            "elements": 1, 
            "loads": 2, # 1 Force + 1 Moment
            "constraints": 1 # SPC1 en nodo 1
        },
        "checks": {
            "nodes": {
                1: [1.0, 0.0, 1.0],
                2: [0.0, 3.0, 2.0]
            },
            "elements": {
                101: {"nodes": [1, 2], "pid": 401, "expected_class": BarElement}
            },
            # Validamos Propiedades Físicas (Solicitado)
            "properties": {
                401: {"A": 100.0, "I1": 833.33, "I2": 833.33, "J": 1406.0}
            },
            "materials": {
                501: {"E": 210000.0, "nu": 0.3}
            },
            "loads": [
                # FORCE 701, Node 2, Mag 1000, Vec [1, 5, -2]
                {"type": "FORCE", "sid": 701, "node": 2, "mag": 1000.0, "vec": [1., 5., -2.]},
                # MOMENT 701, Node 2, Mag 500, Vec [0, 1, 0]
                {"type": "MOMENT", "sid": 701, "node": 2, "mag": 500.0, "vec": [0., 1., 0.]}
            ],
            "constraints": [
                {"sid": 201, "node": 1, "comps": "123456"}
            ]
        }
    }, 

    # --- PENDIENTE DE AÑADIR: CASO 15: Vigas en Serie con Stress Points ---
    # Descripcion: Dos vigas (CBAR) alineadas en eje Y global.
    #              Cargas combinadas en extremo (Fx, Fy, Fz, Mx, My, Mz).
    #              Validacion de Puntos de Tension (C, D, E, F) definidos en PBAR.
    
    # --- CASO 16: Estructura Marco 3D (5 Barras) ---
    {
        "case_name": "Ej16_3D_Frame_5_Bars",
        "folder": "case_16",
        "filename": "Ej16.dat",
        "expected_counts": {
            "nodes": 5, 
            "elements": 5, 
            "loads": 2, 
            "constraints": 2 # SPC1 en nodos 1 y 2
        },
        "checks": {
            "nodes": {
                1: [2.0, 0.0, 2.0],
                3: [2.0, 2.0, 3.0],
                5: [2.0, 3.0, 3.0]
            },
            "elements": {
                101: {"nodes": [1, 3], "expected_class": BarElement},
                105: {"nodes": [4, 5], "expected_class": BarElement}
            },
            "properties": {
                401: {"A": 50.0, "I1": 200.0, "J": 100.0}
            },
            "materials": {
                501: {"E": 210000.0, "nu": 0.3}
            },
            "loads": [
                # FORCE 701, Node 3, Mag 1000, Vec [2, 1, -3]
                {"type": "FORCE", "sid": 701, "node": 3, "mag": 1000.0, "vec": [2., 1., -3.]},
                # FORCE 701, Node 4, Mag 1000, Vec [2, -4, -5]
                {"type": "FORCE", "sid": 701, "node": 4, "mag": 1000.0, "vec": [2., -4., -5.]}
            ],
            "constraints": [
                {"sid": 201, "node": 1, "comps": "123456"},
                {"sid": 201, "node": 2, "comps": "123456"}
            ]
        }
    },

    # --- CASO 17: Trípode Híbrido (CBAR + CROD) ---
    {
        "case_name": "Ej17_Hybrid_Tripod",
        "folder": "case_17",
        "filename": "Ej17.dat",
        "expected_counts": {
            "nodes": 5, 
            "elements": 4, # 3 CBAR + 1 CROD
            "loads": 2, 
            "constraints": 4, # SPC1 en 1,2,3 + SPC1 en 5
            "properties": 4   # 3 PBAR + 1 PROD
        },
        "checks": {
            "nodes": {
                4: [2.0, 3.0, 4.0], # Union
                5: [2.0, 3.0, 7.0]  # Top
            },
            "elements": {
                101: {"nodes": [1, 4], "expected_class": BarElement},
                104: {"nodes": [4, 5], "expected_class": RodElement} # El elemento vertical es un CROD
            },
            # Comparativa de Propiedades (3 PBARs distintas + 1 PROD)
            "properties": {
                401: {"A": 50.0, "I1": 200.0, "J": 100.0}, # Pata Robusta
                402: {"A": 30.0, "I1": 100.0, "J": 50.0},  # Pata Media
                403: {"A": 10.0, "I1": 20.0, "J": 10.0},   # Pata Esbelta
                404: {"A": 20.0, "J": 10.0}                # Poste (PROD)
            },
            "materials": {
                501: {"E": 210000.0, "nu": 0.3}
            },
            "loads": [
                {"type": "FORCE", "sid": 701, "node": 5, "mag": 2000.0, "vec": [0., 0., 1.]},
                {"type": "MOMENT", "sid": 701, "node": 4, "mag": 1500.0, "vec": [0., 0., 1.]}
            ],
            "constraints": [
                {"sid": 201, "node": 1, "comps": "123456"}, # Suelo
                {"sid": 201, "node": 5, "comps": "1245"}    # Punta estabilizada
            ]
        }
    }
]



# ==============================================================================
# TEST CLASS
# ==============================================================================
class TestNastranReaderV2Validation:
    """
    @class TestNastranReaderV2Validation
    @brief Parameterized tests to verify correct Nastran file parsing using V2 (pyNastran).
    """

    @pytest.mark.parametrize("case_data", READER_VALIDATION_CASES, ids=lambda x: x["case_name"])
    def test_parse_full_model_v2(self, case_data, shared_validation_dir):
        """
        @brief Loads a real file using parse_nastran_input and verifies the FEMModel content.
        @param shared_validation_dir Fixture pointing to 'tests/solver_tst/validation_cases'
        """
        # 1. Path Construction
        folder = case_data.get("folder", "")
        filename = case_data["filename"]
        file_path = shared_validation_dir / folder / filename

        if not file_path.exists():
            pytest.fail(f"Test file not found at: {file_path}")

        # 2. Parser Execution (USING V2 READER)
        try:
            model = parse_nastran_input(str(file_path))
        except Exception as e:
            pytest.fail(f"Parser V2 failed reading {filename}: {e}")

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
            assert nid in model.nodes, f"Node {nid} missing in model."
            node_obj = model.nodes[nid]
            # pyNastran computes global coords automatically, we verify accuracy
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
            
            # Check correct mapping (e.g. CROD -> RodElement)
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
                # We iterate to find the specific load matching SID and Node
                if (isinstance(load, target_type) and 
                    load.sid == load_check["sid"] and 
                    load.appnode_id == load_check["node"]):
                    
                    found = True
                    # Magnitude check
                    assert load.module == pytest.approx(load_check["mag"])
                    # Vector check (pyNastran should preserve direction vector)
                    np.testing.assert_allclose(load.vector, load_check["vec"], atol=1e-4)
                    break
            
            assert found, f"Expected load not found in V2 model: {load_check}"

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
            
            assert found, f"Expected constraint not found : {constr_check}"