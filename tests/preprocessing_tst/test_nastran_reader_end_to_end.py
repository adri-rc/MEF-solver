"""End-to-end checks for mapping representative Nastran validation decks."""

from pathlib import Path

import numpy as np
import pytest

from src.preprocessing.data_model import (
    BarElement,
    IsotropicMaterial,
    RodElement,
    SinglePointLoad,
    SinglePointMoment,
)
from src.preprocessing.nastran_reader import parse_nastran_input


@pytest.fixture
def shared_validation_dir():
    current_dir = Path(__file__).parent
    target_dir = current_dir.parent / "solver_tst" / "validation_cases"

    if not target_dir.exists():
        pytest.fail(f"Shared validation directory not found at: {target_dir}")

    return target_dir


READER_VALIDATION_CASES = [
    # Case 8: basic 3D rail.
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
            "constraints": 2,
        },
        "checks": {
            "nodes": {
                7: [40.0, 50.0, 30.0],
                8: [10.0, 80.0, 50.0],
            },
            "elements": {
                104: {"nodes": [7, 8], "pid": 401, "expected_class": RodElement},
            },
            "materials": {
                501: {
                    "E": 210000.0,
                    "nu": 0.3,
                    "expected_class": IsotropicMaterial,
                },
            },
            "loads": [
                {
                    "type": "FORCE",
                    "sid": 704,
                    "node": 8,
                    "mag": 11.72,
                    "vec": [-3.0, 3.0, 2.0],
                }
            ],
            "constraints": [
                {"sid": 204, "node": 7, "comps": "123456"},
                {"sid": 204, "node": 8, "comps": "13456"},
            ],
        },
    },

    # Case 10: DAA truss with INCLUDE files and multiple elements.
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

    # Case 12: 3D chain with MOMENT and enforced displacements.
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
    # Case 14: 3D CBAR under biaxial bending and torsion.
    {
        "case_name": "Ej14_CBAR_3D",
        "folder": "case_14",
        "filename": "Ej14.dat",
        "expected_counts": {
            "nodes": 2,
            "elements": 1,
            "loads": 2,
            "constraints": 1,
        },
        "checks": {
            "nodes": {
                1: [1.0, 0.0, 1.0],
                2: [0.0, 3.0, 2.0]
            },
            "elements": {
                101: {"nodes": [1, 2], "pid": 401, "expected_class": BarElement}
            },
            "properties": {
                401: {"A": 100.0, "I1": 833.33, "I2": 833.33, "J": 1406.0}
            },
            "materials": {
                501: {"E": 210000.0, "nu": 0.3}
            },
            "loads": [
                {"type": "FORCE", "sid": 701, "node": 2, "mag": 1000.0, "vec": [1., 5., -2.]},
                {"type": "MOMENT", "sid": 701, "node": 2, "mag": 500.0, "vec": [0., 1., 0.]}
            ],
            "constraints": [
                {"sid": 201, "node": 1, "comps": "123456"}
            ]
        }
    },

    # Case 15: two serial CBAR elements with PBAR stress points.
    {
        "case_name": "Ej15_Series_Beams_StressPoints",
        "folder": "case_15",
        "filename": "Ej15.dat",
        "expected_counts": {
            "nodes": 3,
            "elements": 2,
            "properties": 2,
            "materials": 1,
            "loads": 2,
            "constraints": 1
        },
        "checks": {
            "nodes": {
                1: [0.0, 0.0, 0.0],
                2: [0.0, 2.0, 0.0],
                3: [0.0, 5.0, 0.0]
            },
            "elements": {
                101: {"nodes": [1, 2], "pid": 401, "expected_class": BarElement},
                102: {"nodes": [2, 3], "pid": 402, "expected_class": BarElement}
            },
            "properties": {
                401: {
                    "A": 100.0,
                    "I1": 1000.0,
                    "stress_points": [
                        (10.0, 0.0),
                        (0.0, 10.0),
                        (0.0, 0.0),
                        (0.0, 0.0),
                    ],
                },
                402: {
                    "A": 50.0,
                    "I1": 200.0,
                    "stress_points": [
                        (5.0, 5.0),
                        (-5.0, -5.0),
                        (0.0, 0.0),
                        (0.0, 0.0),
                    ],
                },
            },
            "loads": [
                {"type": "FORCE", "sid": 701, "node": 3, "mag": 1.0, "vec": [1., 2., 4.]},
                {"type": "MOMENT", "sid": 701, "node": 3, "mag": 1.0, "vec": [-1., -3., -2.]}
            ],
            "constraints": [
                {"sid": 201, "node": 1, "comps": "123456"}
            ]
        }
    },

    # Case 16: 3D frame with five CBAR elements.
    {
        "case_name": "Ej16_3D_Frame_5_Bars",
        "folder": "case_16",
        "filename": "Ej16.dat",
        "expected_counts": {
            "nodes": 5,
            "elements": 5,
            "loads": 2,
            "constraints": 2,
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
                {"type": "FORCE", "sid": 701, "node": 3, "mag": 1000.0, "vec": [2., 1., -3.]},
                {"type": "FORCE", "sid": 701, "node": 4, "mag": 1000.0, "vec": [2., -4., -5.]}
            ],
            "constraints": [
                {"sid": 201, "node": 1, "comps": "123456"},
                {"sid": 201, "node": 2, "comps": "123456"}
            ]
        }
    },

    # Case 17: hybrid CBAR/CROD tripod.
    {
        "case_name": "Ej17_Hybrid_Tripod",
        "folder": "case_17",
        "filename": "Ej17.dat",
        "expected_counts": {
            "nodes": 5,
            "elements": 4,
            "loads": 2,
            "constraints": 4,
            "properties": 4,
        },
        "checks": {
            "nodes": {
                4: [2.0, 3.0, 4.0],
                5: [2.0, 3.0, 7.0],
            },
            "elements": {
                101: {"nodes": [1, 4], "expected_class": BarElement},
                104: {"nodes": [4, 5], "expected_class": RodElement},
            },
            "properties": {
                401: {"A": 50.0, "I1": 200.0, "J": 100.0},
                402: {"A": 30.0, "I1": 100.0, "J": 50.0},
                403: {"A": 10.0, "I1": 20.0, "J": 10.0},
                404: {"A": 20.0, "J": 10.0},
            },
            "materials": {
                501: {"E": 210000.0, "nu": 0.3}
            },
            "loads": [
                {"type": "FORCE", "sid": 701, "node": 5, "mag": 2000.0, "vec": [0., 0., 1.]},
                {"type": "MOMENT", "sid": 701, "node": 4, "mag": 1500.0, "vec": [0., 0., 1.]}
            ],
            "constraints": [
                {"sid": 201, "node": 1, "comps": "123456"},
                {"sid": 201, "node": 5, "comps": "1245"},
            ]
        }
    }
]
class TestNastranReaderValidation:
    @pytest.mark.parametrize(
        "case_data",
        READER_VALIDATION_CASES,
        ids=lambda case: case["case_name"],
    )
    def test_parse_full_model(self, case_data, shared_validation_dir):
        folder = case_data.get("folder", "")
        filename = case_data["filename"]
        file_path = shared_validation_dir / folder / filename

        if not file_path.exists():
            pytest.fail(f"Test file not found at: {file_path}")

        try:
            model = parse_nastran_input(str(file_path))
        except Exception as exc:
            pytest.fail(f"Parser failed reading {filename}: {exc}")

        counts = case_data.get("expected_counts", {})
        for entity_name in (
            "nodes",
            "elements",
            "materials",
            "properties",
            "loads",
            "constraints",
        ):
            if entity_name not in counts:
                continue
            actual_count = len(getattr(model, entity_name))
            assert actual_count == counts[entity_name], (
                f"{entity_name} count mismatch: "
                f"{actual_count} != {counts[entity_name]}"
            )
        if "constraints_min" in counts:
            assert len(model.constraints) >= counts["constraints_min"]

        checks = case_data.get("checks", {})

        for nid, expected_coords in checks.get("nodes", {}).items():
            assert nid in model.nodes, f"Node {nid} missing in model."
            node_obj = model.nodes[nid]
            np.testing.assert_allclose(
                node_obj.get_coords(),
                expected_coords,
                atol=1e-5,
                err_msg=f"Incorrect coords for Node {nid}",
            )

        for eid, data in checks.get("elements", {}).items():
            assert eid in model.elements, f"Element {eid} missing."
            elem = model.elements[eid]

            if "nodes" in data:
                assert elem.nodes == data["nodes"]
            if "pid" in data:
                assert elem.prop_id == data["pid"]

            if "expected_class" in data:
                assert isinstance(elem, data["expected_class"]), (
                    f"Element {eid} is {type(elem).__name__}; expected "
                    f"{data['expected_class'].__name__}"
                )

        for mid, props in checks.get("materials", {}).items():
            assert mid in model.materials, f"Material {mid} missing."
            mat = model.materials[mid]

            if "E" in props:
                assert mat.E == pytest.approx(props["E"])
            if "nu" in props:
                assert mat.nu == pytest.approx(props["nu"])

            if "expected_class" in props:
                assert isinstance(mat, props["expected_class"]), (
                    f"Material {mid} is {type(mat).__name__}; expected "
                    f"{props['expected_class'].__name__}"
                )

        for pid, expected in checks.get("properties", {}).items():
            assert pid in model.properties, f"Property {pid} missing."
            prop = model.properties[pid]
            for attribute, expected_value in expected.items():
                actual_value = getattr(prop, attribute)
                if attribute == "stress_points":
                    np.testing.assert_allclose(actual_value, expected_value)
                else:
                    assert actual_value == pytest.approx(expected_value)

        for load_check in checks.get("loads", []):
            found = False
            target_type = (
                SinglePointLoad
                if load_check["type"] == "FORCE"
                else SinglePointMoment
            )

            for load in model.loads:
                if (
                    isinstance(load, target_type)
                    and load.sid == load_check["sid"]
                    and load.appnode_id == load_check["node"]
                ):
                    found = True
                    assert load.module == pytest.approx(load_check["mag"])
                    np.testing.assert_allclose(
                        load.vector, load_check["vec"], atol=1e-4
                    )
                    break

            assert found, f"Expected load not found in model: {load_check}"

        for constr_check in checks.get("constraints", []):
            found = False
            for spc in model.constraints:
                if (
                    spc.sid == constr_check["sid"]
                    and spc.node_id == constr_check["node"]
                    and spc.components == constr_check["comps"]
                ):
                    if "value" in constr_check:
                        assert spc.value == pytest.approx(constr_check["value"]), (
                            f"Incorrect SPC value at Node {spc.node_id}"
                        )

                    found = True
                    break

            assert found, f"Expected constraint not found: {constr_check}"
