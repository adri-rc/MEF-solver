"""E2E validation of assembled CQUAD4 stiffness matrices against Nastran.

The originally requested case 21 is a multi-element CTRIA3 model and has no
KAAX/DMIG matrix reference.  Case 51 is the corresponding CQUAD4 membrane
case and is therefore used together with CQUAD4 cases 52, 53, 54 and 55.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.preprocessing.data_model import CQUAD4Element
from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.elements.cquad4_nastran_element import (
    CQUAD4NastranElementSolver,
)
from src.solver.linear_static_solver import LinearStaticSolver


STIFFNESS_VALIDATION_DATA = [
    {
        "case_name": "Ej51_KGlobal_CQUAD4_Membrane_K6ROT_0.0",
        "folder": "case_51",
        "filename": "ej51_kmembrane_k6rot_0.dat",
        "k_csv_file": "ej51_kmembrane_k6rot_0_Nastran.csv",
        # MID2 is blank: Nastran leaves the complete plate/shear subspace at
        # zero.  Comparing all 24 DOFs prevents a fictitious MID1 -> MID2
        # fallback from being hidden by a membrane-only slice.
        "active_dofs": list(range(24)),
        "comparison_blocks": {
            "membrane [T1,T2]": [0, 1, 6, 7, 12, 13, 18, 19],
            "plate/shear [T3,R1,R2]": [
                2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22
            ],
        },
        "max_relative_frobenius_error": 0.01,
        "debug_export": True,
        "output_k_file": "draft_Ej51_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej52_KGlobal_CQUAD4_PureBending_K6ROT_0.0",
        "folder": "case_52",
        "filename": "ej52_k_pure_flex_mid3_0_k6rot_0.dat",
        "k_csv_file": "ej52_k_pure_flex_mid3_0_k6rot_0_Nastran.csv",
        "active_dofs": [2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22],
        "max_relative_frobenius_error": 0.01,
        "debug_export": True,
        "output_k_file": "draft_Ej52_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej53_KGlobal_CQUAD4_MembraneBending_K6ROT_0.0",
        "folder": "case_53",
        "filename": "ej53_k_memb-flex_mid3_0_k6rot_0.dat",
        "k_csv_file": "ej53_k_memb-flex_mid3_0_k6rot_0_Nastran.csv",
        "active_dofs": [
            0, 1, 2, 3, 4,
            6, 7, 8, 9, 10,
            12, 13, 14, 15, 16,
            18, 19, 20, 21, 22,
        ],
        "max_relative_frobenius_error": 0.01,
        "debug_export": True,
        "output_k_file": "draft_Ej53_KGlobal_CQUAD4_python_test.csv",
    },
    pytest.param(
        {
        "case_name": "Ej54_KGlobal_CQUAD4_Warped3D_K6ROT_0.0",
        "folder": "case_54",
        "filename": "ej54_k_memb-flex-3d_mid3_0_k6rot_0.dat",
        "k_csv_file": "ej54_k_memb-flex-3d_mid3_0_k6rot_0_Nastran.csv",
        "active_dofs": list(range(24)),
        "max_relative_frobenius_error": 0.03,
        "debug_export": True,
        "output_k_file": "draft_Ej54_KGlobal_CQUAD4_python_test.csv",
        },
        marks=pytest.mark.xfail(
            strict=True,
            reason=(
                "CQUAD4 is intentionally limited to coplanar nodes; this "
                "historical warped DMIG deck must be rejected."
            ),
        ),
    ),
    {
        "case_name": "Ej55_KGlobal_CQUAD4_Inclined3D_K6ROT_100.0",
        "folder": "case_55",
        "filename": "ej55_k_memb-flex-3d_mid3_0_k6rot_100.dat",
        "k_csv_file": "ej55_k_memb-flex-3d_mid3_0_k6rot_100_Nastran.csv",
        "active_dofs": list(range(24)),
        "max_relative_frobenius_error": 0.03,
        "validate_drilling_subspace": True,
        "debug_export": True,
        "output_k_file": "draft_Ej55_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej65_KGlobal_CQUAD4_TransverseShear_K6ROT_0.0",
        "folder": "case_65",
        "filename": "ej65_cquad4_transverse_shear_stiffness_k6rot_0.dat",
        "k_csv_file": (
            "ej65_cquad4_transverse_shear_stiffness_k6rot_0_Nastran.csv"
        ),
        "active_dofs": list(range(24)),
        "max_relative_frobenius_error": 1.0e-8,
        "term_rtol": 1.0e-8,
        "term_atol": 1.0e-10,
        "comparison_blocks": {
            "membrane [T1,T2]": [0, 1, 6, 7, 12, 13, 18, 19],
            "plate/shear [T3,R1,R2]": [
                2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22
            ],
        },
        "debug_export": True,
        "output_k_file": "draft_Ej65_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej67_KGlobal_CQUAD4_FullShear_K6ROT_100.0",
        "folder": "case_67",
        "filename": "ej67_cquad4_full_stiffness_shear_k6rot_100.dat",
        "k_csv_file": (
            "ej67_cquad4_full_stiffness_shear_k6rot_100_Nastran.csv"
        ),
        "active_dofs": list(range(24)),
        "max_relative_frobenius_error": 1.0e-5,
        "term_rtol": 1.0e-5,
        "term_atol": 1.0e-1,
        "comparison_blocks": {
            "membrane [T1,T2]": [0, 1, 6, 7, 12, 13, 18, 19],
            "plate/shear [T3,R1,R2]": [
                2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22
            ],
            "drilling [R3]": [5, 11, 17, 23],
        },
        "comparison_couplings": {
            "drilling-membrane [R3]x[T1,T2]": {
                "rows": [5, 11, 17, 23],
                "columns": [0, 1, 6, 7, 12, 13, 18, 19],
            }
        },
        "debug_export": True,
        "output_k_file": "draft_Ej67_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej69_KGlobal_CQUAD4_MID3Blank_K6ROT_0.0",
        "folder": "case_69",
        "filename": "ej69_cquad4_transverse_shear_baseline_mid3_0_k6rot_0.dat",
        "k_csv_file": "ej69_cquad4_transverse_shear_baseline_mid3_0_k6rot_0_Nastran.csv",
        "active_dofs": list(range(24)),
        "max_relative_frobenius_error": 1.0e-8,
        "term_rtol": 1.0e-8,
        "term_atol": 1.0e-10,
        "debug_export": True,
        "output_k_file": "draft_Ej69_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej70_KGlobal_CQUAD4_G3Half_K6ROT_0.0",
        "folder": "case_70",
        "filename": "ej70_cquad4_transverse_shear_g3_half_k6rot_0.dat",
        "k_csv_file": "ej70_cquad4_transverse_shear_g3_half_k6rot_0_Nastran.csv",
        "active_dofs": list(range(24)),
        "max_relative_frobenius_error": 1.0e-8,
        "term_rtol": 1.0e-8,
        "term_atol": 1.0e-10,
        "debug_export": True,
        "output_k_file": "draft_Ej70_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej71_KGlobal_CQUAD4_G3Double_K6ROT_0.0",
        "folder": "case_71",
        "filename": "ej71_cquad4_transverse_shear_g3_double_k6rot_0.dat",
        "k_csv_file": "ej71_cquad4_transverse_shear_g3_double_k6rot_0_Nastran.csv",
        "active_dofs": list(range(24)),
        "max_relative_frobenius_error": 1.0e-8,
        "term_rtol": 1.0e-8,
        "term_atol": 1.0e-10,
        "debug_export": True,
        "output_k_file": "draft_Ej71_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej72_KGlobal_CQUAD4_MID3Blank_K6ROT_100.0",
        "folder": "case_72",
        "filename": "ej72_cquad4_full_stiffness_mid3_0_k6rot_100.dat",
        "k_csv_file": "ej72_cquad4_full_stiffness_mid3_0_k6rot_100_Nastran.csv",
        "active_dofs": list(range(24)),
        "max_relative_frobenius_error": 1.0e-5,
        "term_rtol": 1.0e-5,
        "term_atol": 1.0e-1,
        "debug_export": True,
        "output_k_file": "draft_Ej72_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej73_KGlobal_CQUAD4_TSTHalf_K6ROT_0.0",
        "folder": "case_73",
        "filename": "ej73_cquad4_transverse_shear_tst_half_k6rot_0.dat",
        "k_csv_file": "ej73_cquad4_transverse_shear_tst_half_k6rot_0_Nastran.csv",
        "active_dofs": list(range(24)),
        "max_relative_frobenius_error": 1.0e-8,
        "term_rtol": 1.0e-8,
        "term_atol": 1.0e-10,
        "debug_export": True,
        "output_k_file": "draft_Ej73_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej74_KGlobal_CQUAD4_SkewMID3_K6ROT_0.0",
        "folder": "case_74",
        "filename": "ej74_cquad4_skew_transverse_shear_mid3_active_k6rot_0.dat",
        "k_csv_file": "ej74_cquad4_skew_transverse_shear_mid3_active_k6rot_0_Nastran.csv",
        "active_dofs": list(range(24)),
        "max_relative_frobenius_error": 0.01,
        "term_rtol": 0.01,
        "term_atol": 1.0e-3,
        "debug_export": True,
        "output_k_file": "draft_Ej74_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej75_KGlobal_CQUAD4_SkewMID3Blank_K6ROT_0.0",
        "folder": "case_75",
        "filename": "ej75_cquad4_skew_transverse_shear_mid3_0_k6rot_0.dat",
        "k_csv_file": "ej75_cquad4_skew_transverse_shear_mid3_0_k6rot_0_Nastran.csv",
        "active_dofs": list(range(24)),
        "max_relative_frobenius_error": 0.01,
        "term_rtol": 0.01,
        "term_atol": 1.0e-3,
        "debug_export": True,
        "output_k_file": "draft_Ej75_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej81_KGlobal_CQUAD4_SkewAffineFull_K6ROT_100.0",
        "folder": "case_81",
        "filename": "ej81_cquad4_skew_parallelogram_full_dmig_k6rot_100.dat",
        "k_csv_file": (
            "ej81_cquad4_skew_parallelogram_full_dmig_k6rot_100_Nastran.csv"
        ),
        "active_dofs": list(range(24)),
        # This full-property affine parallelogram validates the oblique
        # augmented-shear metric.  The complete matrix now differs by only
        # 1.03e-7 and its isolated plate/shear block by 3.15e-6.
        "max_relative_frobenius_error": 2.0e-7,
        "term_rtol": 1.0e-5,
        "term_atol": 2.0e-3,
        "block_coordinate_system": "element",
        "comparison_blocks": {
            "membrane [T1,T2]": [0, 1, 6, 7, 12, 13, 18, 19],
            "plate/shear [T3,R1,R2]": [
                2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22
            ],
            "drilling [R3]": [5, 11, 17, 23],
        },
        "comparison_couplings": {
            "drilling-membrane [R3]x[T1,T2]": {
                "rows": [5, 11, 17, 23],
                "columns": [0, 1, 6, 7, 12, 13, 18, 19],
            },
        },
        "comparison_coupling_limits": {
            "drilling-membrane [R3]x[T1,T2]": {
                "max_relative_frobenius_error": 1.0e-7,
                "term_rtol": 1.0e-7,
                "term_atol": 1.0e-8,
            },
        },
        "comparison_block_limits": {
            "membrane [T1,T2]": {
                "max_relative_frobenius_error": 1.0e-8,
                "term_rtol": 1.0e-8,
                "term_atol": 1.0e-5,
            },
            "plate/shear [T3,R1,R2]": {
                "max_relative_frobenius_error": 4.0e-6,
                "term_rtol": 1.0e-5,
                "term_atol": 2.0e-3,
            },
            "drilling [R3]": {
                "max_relative_frobenius_error": 1.0e-7,
                "term_rtol": 1.0e-7,
                "term_atol": 1.0e-8,
            },
        },
        "debug_export": True,
        "output_k_file": "draft_Ej81_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej83_KGlobal_CQUAD4_SkewAffineMID3Isolated_K6ROT_100.0",
        "folder": "case_83",
        "filename": (
            "ej83_cquad4_skew_affine_mid3_isolation_dmig_k6rot_100.dat"
        ),
        "k_csv_file": (
            "ej83_cquad4_skew_affine_mid3_isolation_dmig_k6rot_100_Nastran.csv"
        ),
        "active_dofs": list(range(24)),
        # MID1 is deliberately soft while MID2/MID3 stay stiff.  This makes
        # the transverse-shear flexibility dominate the DMIG comparison.
        "max_relative_frobenius_error": 4.0e-6,
        "term_rtol": 1.0e-5,
        "term_atol": 2.0e-3,
        "block_coordinate_system": "element",
        "comparison_blocks": {
            "membrane [T1,T2]": [0, 1, 6, 7, 12, 13, 18, 19],
            "plate/shear [T3,R1,R2]": [
                2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22
            ],
            "drilling [R3]": [5, 11, 17, 23],
        },
        "comparison_block_limits": {
            "membrane [T1,T2]": {
                "max_relative_frobenius_error": 5.0e-7,
                "term_rtol": 1.0e-5,
                "term_atol": 1.0e-6,
            },
            "plate/shear [T3,R1,R2]": {
                "max_relative_frobenius_error": 4.0e-6,
                "term_rtol": 1.0e-5,
                "term_atol": 2.0e-3,
            },
            "drilling [R3]": {
                # Its norm is intentionally tiny because it derives from
                # soft MID1; the absolute term check is the useful metric.
                "max_relative_frobenius_error": 1.0e-3,
                "term_rtol": 1.0e-3,
                "term_atol": 1.0e-8,
            },
        },
        "debug_export": True,
        "output_k_file": "draft_Ej83_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej89_KGlobal_CQUAD4_SkewAffineMID3Blank_K6ROT_100.0",
        "folder": "case_89",
        "filename": "ej89_cquad4_skew_affine_mid3_blank_dmig_k6rot_100.dat",
        "k_csv_file": (
            "ej89_cquad4_skew_affine_mid3_blank_dmig_k6rot_100_Nastran.csv"
        ),
        "active_dofs": list(range(24)),
        # MID3 blank disables physical transverse-shear flexibility but not
        # MacNeal's residual bending/shear correction.  It is therefore the
        # control paired with case 83.
        "max_relative_frobenius_error": 2.0e-6,
        "term_rtol": 1.0e-5,
        "term_atol": 2.0e-3,
        "block_coordinate_system": "element",
        "comparison_blocks": {
            "membrane [T1,T2]": [0, 1, 6, 7, 12, 13, 18, 19],
            "plate/shear [T3,R1,R2]": [
                2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22
            ],
            "drilling [R3]": [5, 11, 17, 23],
        },
        "comparison_block_limits": {
            "membrane [T1,T2]": {
                "max_relative_frobenius_error": 2.0e-6,
                "term_rtol": 1.0e-5,
                "term_atol": 1.0e-6,
            },
            "plate/shear [T3,R1,R2]": {
                "max_relative_frobenius_error": 2.0e-6,
                "term_rtol": 1.0e-5,
                "term_atol": 2.0e-3,
            },
            "drilling [R3]": {
                "max_relative_frobenius_error": 4.0e-3,
                "term_rtol": 4.0e-3,
                "term_atol": 3.0e-8,
            },
        },
        "debug_export": True,
        "output_k_file": "draft_Ej89_KGlobal_CQUAD4_python_test.csv",
    },
    {
        "case_name": "Ej90_KGlobal_CQUAD4_HighSkewAffine_K6ROT_100.0",
        "folder": "case_90",
        "filename": "ej90_cquad4_high_skew_affine_dmig_k6rot_100.dat",
        "k_csv_file": (
            "ej90_cquad4_high_skew_affine_dmig_k6rot_100_Nastran.csv"
        ),
        "active_dofs": list(range(24)),
        # High-skew calibration control: its 50.6-degree midline angle exposed
        # that the former one-point continuation was too soft. Together with
        # cases 81/89 it fixes the bounded rational torsion metric; loaded
        # case 60 checks an intermediate angle independently.
        "max_relative_frobenius_error": 1.0e-6,
        "term_rtol": 1.0e-5,
        "term_atol": 2.0e-3,
        "block_coordinate_system": "element",
        "comparison_blocks": {
            "membrane [T1,T2]": [0, 1, 6, 7, 12, 13, 18, 19],
            "plate/shear [T3,R1,R2]": [
                2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22
            ],
            "drilling [R3]": [5, 11, 17, 23],
        },
        "comparison_block_limits": {
            "membrane [T1,T2]": {
                "max_relative_frobenius_error": 1.0e-6,
                "term_rtol": 1.0e-5,
                "term_atol": 2.0e-3,
            },
            "plate/shear [T3,R1,R2]": {
                "max_relative_frobenius_error": 4.0e-6,
                "term_rtol": 1.0e-5,
                "term_atol": 2.0e-3,
            },
            "drilling [R3]": {
                "max_relative_frobenius_error": 1.0e-6,
                "term_rtol": 1.0e-5,
                "term_atol": 1.0e-6,
            },
        },
        "debug_export": True,
        "output_k_file": "draft_Ej90_KGlobal_CQUAD4_python_test.csv",
    },
]


class TestCQUAD4StiffnessMatrix:
    """Compare the solver's assembled CQUAD4 K matrix with Nastran KAAX."""

    @staticmethod
    def _export_stiffness_matrix(K: np.ndarray, output_path: Path) -> None:
        """Write the complete calculated matrix using the Nastran CSV layout."""
        dof_labels = []
        for node_index in range(1, K.shape[0] // 6 + 1):
            dof_labels.extend(
                [
                    f"u{node_index}",
                    f"v{node_index}",
                    f"w{node_index}",
                    f"theta_x{node_index}",
                    f"theta_y{node_index}",
                    f"theta_z{node_index}",
                ]
            )

        dataframe = pd.DataFrame(K, columns=dof_labels)
        dataframe.insert(0, "DoF", dof_labels)
        dataframe.insert(0, "Index", np.arange(K.shape[0]))
        dataframe.to_csv(output_path, index=False, float_format="%.9e")
        print(f"\n[DEBUG] Matriz CSV del solver exportada a: {output_path}")

    @pytest.mark.parametrize(
        "case_data",
        STIFFNESS_VALIDATION_DATA,
        ids=lambda case: case["case_name"],
    )
    def test_global_stiffness_matrix(self, case_data, validation_data_dir):
        case_path = (
            validation_data_dir.parent.parent
            / "validation_cases"
            / case_data["folder"]
        )
        dat_path = case_path / case_data["filename"]
        csv_path = case_path / case_data["k_csv_file"]

        if not dat_path.is_file():
            pytest.fail(f"Falta el archivo de entrada Nastran: {dat_path}")
        if not csv_path.is_file():
            pytest.fail(f"Falta el CSV de rigidez de Nastran: {csv_path}")

        reference = pd.read_csv(csv_path).iloc[:, 2:].to_numpy(dtype=float)
        model = parse_nastran_input(str(dat_path))
        solver = LinearStaticSolver(model)
        solver.build_stiffness_matrix()
        calculated = solver.K_global.toarray()

        if case_data.get("debug_export", False):
            output_k_file = case_data.get("output_k_file")
            if not output_k_file:
                pytest.fail("debug_export requiere configurar output_k_file")
            self._export_stiffness_matrix(
                calculated, case_path / output_k_file
            )

        assert calculated.shape == reference.shape, (
            "La dimensión de K_global no coincide: "
            f"solver={calculated.shape}, Nastran={reference.shape}"
        )
        np.testing.assert_allclose(
            calculated,
            calculated.T,
            rtol=1.0e-12,
            atol=1.0e-8,
            err_msg="K_global calculada no es simétrica",
        )

        active_dofs = case_data["active_dofs"]
        calculated_active = calculated[np.ix_(active_dofs, active_dofs)]
        reference_active = reference[np.ix_(active_dofs, active_dofs)]
        reference_norm = np.linalg.norm(reference_active)
        assert reference_norm > 0.0, "La matriz de referencia activa es nula"

        relative_error = (
            np.linalg.norm(calculated_active - reference_active) / reference_norm
        )
        tolerance = case_data["max_relative_frobenius_error"]

        block_calculated = calculated
        block_reference = reference
        if case_data.get("block_coordinate_system") == "element":
            cquad_elements = [
                element
                for element in model.elements.values()
                if isinstance(element, CQUAD4Element)
            ]
            assert len(cquad_elements) == 1, (
                "La comparación en ejes del elemento requiere un solo CQUAD4"
            )
            element = cquad_elements[0]
            prop = model.properties[element.prop_id]
            element_solver = CQUAD4NastranElementSolver(
                *(model.get_node(node_id) for node_id in element.nodes),
                prop,
                mat_membrane=model.materials[prop.mat_id],
                mat_bending=model.materials.get(prop.mat_id_bending),
                mat_shear=model.materials.get(prop.mat_id_shear),
                k6rot=solver.k6rot,
                theta=element.theta,
            )
            transform = element_solver.T
            block_calculated = transform @ calculated @ transform.T
            block_reference = transform @ reference @ transform.T

        block_diagnostics = []
        for block_name, block_dofs in case_data.get(
            "comparison_blocks", {}
        ).items():
            calculated_block = block_calculated[np.ix_(block_dofs, block_dofs)]
            reference_block = block_reference[np.ix_(block_dofs, block_dofs)]
            block_reference_norm = np.linalg.norm(reference_block)
            if block_reference_norm == 0.0:
                block_error = np.linalg.norm(calculated_block)
                error_label = "error absoluto"
            else:
                block_error = (
                    np.linalg.norm(calculated_block - reference_block)
                    / block_reference_norm
                )
                error_label = "error relativo de Frobenius"
            block_diagnostics.append(
                f"  - {block_name}: {error_label}={block_error:.6e}"
            )
            block_limit = case_data.get(
                "comparison_block_limits", {}
            ).get(block_name)
            if block_limit is not None:
                assert block_error <= block_limit[
                    "max_relative_frobenius_error"
                ], (
                    f"El bloque {block_name} excede su tolerancia en "
                    f"{case_data['case_name']}: error={block_error:.6e}, "
                    "límite="
                    f"{block_limit['max_relative_frobenius_error']:.6e}"
                )
                if "term_rtol" in block_limit:
                    np.testing.assert_allclose(
                        calculated_block,
                        reference_block,
                        rtol=block_limit["term_rtol"],
                        atol=block_limit["term_atol"],
                        err_msg=(
                            f"El bloque {block_name} no coincide término a "
                            f"término en {case_data['case_name']}"
                        ),
                    )

        for coupling_name, coupling_dofs in case_data.get(
            "comparison_couplings", {}
        ).items():
            rows = coupling_dofs["rows"]
            columns = coupling_dofs["columns"]
            calculated_coupling = block_calculated[np.ix_(rows, columns)]
            reference_coupling = block_reference[np.ix_(rows, columns)]
            coupling_reference_norm = np.linalg.norm(reference_coupling)
            if coupling_reference_norm == 0.0:
                coupling_error = np.linalg.norm(calculated_coupling)
                error_label = "error absoluto"
            else:
                coupling_error = (
                    np.linalg.norm(calculated_coupling - reference_coupling)
                    / coupling_reference_norm
                )
                error_label = "error relativo de Frobenius"
            block_diagnostics.append(
                f"  - {coupling_name}: {error_label}={coupling_error:.6e}"
            )
            coupling_limit = case_data.get(
                "comparison_coupling_limits", {}
            ).get(coupling_name)
            if coupling_limit is not None:
                assert coupling_error <= coupling_limit[
                    "max_relative_frobenius_error"
                ], (
                    f"El acoplamiento {coupling_name} excede su tolerancia "
                    f"en {case_data['case_name']}: error={coupling_error:.6e}, "
                    "límite="
                    f"{coupling_limit['max_relative_frobenius_error']:.6e}"
                )
                if "term_rtol" in coupling_limit:
                    np.testing.assert_allclose(
                        calculated_coupling,
                        reference_coupling,
                        rtol=coupling_limit["term_rtol"],
                        atol=coupling_limit["term_atol"],
                        err_msg=(
                            f"El acoplamiento {coupling_name} no coincide "
                            f"término a término en {case_data['case_name']}"
                        ),
                    )

        absolute_differences = np.abs(calculated_active - reference_active)
        worst_flat_index = int(np.argmax(absolute_differences))
        worst_local_row, worst_local_column = np.unravel_index(
            worst_flat_index, absolute_differences.shape
        )
        worst_row = active_dofs[worst_local_row]
        worst_column = active_dofs[worst_local_column]
        diagnostics = ""
        if block_diagnostics:
            diagnostics += "\nDiagnostico por bloques:\n" + "\n".join(
                block_diagnostics
            )
        diagnostics += (
            "\nMayor diferencia termino a termino: "
            f"K[{worst_row},{worst_column}] solver="
            f"{calculated[worst_row, worst_column]:.9e}, Nastran="
            f"{reference[worst_row, worst_column]:.9e}, diferencia="
            f"{absolute_differences[worst_local_row, worst_local_column]:.9e}"
        )
        assert relative_error < tolerance, (
            f"K_global no coincide para {case_data['case_name']}: "
            f"error relativo de Frobenius={relative_error:.6e}, "
            f"tolerancia={tolerance:.6e}{diagnostics}"
        )

        if "term_rtol" in case_data:
            np.testing.assert_allclose(
                calculated_active,
                reference_active,
                rtol=case_data["term_rtol"],
                atol=case_data["term_atol"],
                err_msg=(
                    "La matriz activa no coincide termino a termino en "
                    f"{case_data['case_name']}"
                ),
            )

        if case_data.get("validate_drilling_subspace", False):
            element = next(iter(model.elements.values()))
            coordinates = np.asarray(
                [model.get_node(nid).get_coords() for nid in element.nodes]
            )
            normal = np.cross(
                coordinates[2] - coordinates[0],
                coordinates[3] - coordinates[1],
            )
            normal /= np.linalg.norm(normal)

            projection = np.zeros((calculated.shape[0], len(element.nodes)))
            for column, node_id in enumerate(element.nodes):
                node_dofs = solver.dof_mapper.get_dof_indices(node_id)
                projection[node_dofs[3:6], column] = normal

            # Project each nodal normal rotation against every global DOF.
            # This isolates both theta_z-theta_z and theta_z-macro-rotation
            # terms, preventing the much larger membrane block from hiding
            # an incorrect K6ROT distribution in the global Frobenius norm.
            calculated_drilling = projection.T @ calculated
            reference_drilling = projection.T @ reference
            np.testing.assert_allclose(
                calculated_drilling,
                reference_drilling,
                rtol=1.0e-5,
                atol=1.0e-4,
                err_msg=(
                    "El subespacio drilling K6ROT no coincide termino a "
                    f"termino en {case_data['case_name']}"
                ),
            )
