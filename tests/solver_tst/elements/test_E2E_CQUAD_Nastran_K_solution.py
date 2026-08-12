"""E2E validation of assembled CQUAD4 stiffness matrices against Nastran.

The originally requested case 21 is a multi-element CTRIA3 model and has no
KAAX/DMIG matrix reference.  Case 51 is the corresponding CQUAD4 membrane
case and is therefore used together with CQUAD4 cases 52, 53, 54 and 55.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.linear_static_solver import LinearStaticSolver


STIFFNESS_VALIDATION_DATA = [
    {
        "case_name": "Ej51_KGlobal_CQUAD4_Membrane_K6ROT_0.0",
        "folder": "case_51",
        "filename": "ej51_kmembrane_k6rot_0.dat",
        "k_csv_file": "ej51_kmembrane_k6rot_0_Nastran.csv",
        "active_dofs": [0, 1, 6, 7, 12, 13, 18, 19],
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
        assert relative_error < tolerance, (
            f"K_global no coincide para {case_data['case_name']}: "
            f"error relativo de Frobenius={relative_error:.6e}, "
            f"tolerancia={tolerance:.6e}"
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
