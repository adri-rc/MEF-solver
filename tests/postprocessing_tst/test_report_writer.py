from types import SimpleNamespace

import numpy as np

from src.postprocessing.report_writer import write_report


class SequentialDofMapper:
    def get_dof_indices(self, node_id):
        start = (node_id - 1) * 6
        return list(range(start, start + 6))


def _shell_results(element_type):
    return {
        "type": element_type,
        "membrane_force_x": 1.0,
        "membrane_force_y": 2.0,
        "membrane_force_xy": 3.0,
        "bending_moment_x": 4.0,
        "bending_moment_y": 5.0,
        "bending_moment_xy": 6.0,
        "shear_force_x": 7.0,
        "shear_force_y": 8.0,
        "stress_top_fiber_x": 11.0,
        "stress_top_fiber_y": 12.0,
        "tau_top_fiber_xy": 13.0,
        "major_principal_top": 14.0,
        "minor_principal_top": 15.0,
        "von_mises_top": 16.0,
        "stress_bot_fiber_x": 21.0,
        "stress_bot_fiber_y": 22.0,
        "tau_bot_fiber_xy": 23.0,
        "major_principal_bot": 24.0,
        "minor_principal_bot": 25.0,
        "von_mises_bot": 26.0,
    }


def test_report_separates_element_families_and_uses_physical_headers(tmp_path):
    model = SimpleNamespace(
        nodes={1: object(), 2: object()},
        materials={},
        properties={},
        elements={},
        loads=[],
        constraints=[],
    )
    solver = SimpleNamespace(dof_mapper=SequentialDofMapper())
    cbar_point_a = {
        "P0": {
            "y": 0.1, "z": 0.2, "sigma_axial": 30.0,
            "sigma_bending": 31.0, "sigma_total": 61.0, "tau_xy": 5.0,
        }
    }
    cbar_point_b = {
        "P0": {
            "y": 0.1, "z": 0.2, "sigma_axial": 30.0,
            "sigma_bending": -31.0, "sigma_total": -1.0, "tau_xy": 5.0,
        }
    }
    results = {
        "displacements": np.arange(12, dtype=float) * 1.0e-3,
        "reactions": {1: np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])},
        "stresses": {
            10: {
                "type": "CROD", "axial_force": 100.0, "torque": 25.0,
                "axial_stress": 50.0, "torsional_stress": 12.5,
                "axial_strain": 0.01, "twist_rate": 0.02,
            },
            20: {
                "type": "CBAR",
                "forces_n1": {"Fx": 1.0, "Fy": 2.0, "Fz": 3.0, "Mx": 4.0, "My": 5.0, "Mz": 6.0},
                "forces_n2": {"Fx": -1.0, "Fy": -2.0, "Fz": -3.0, "Mx": -4.0, "My": -5.0, "Mz": -6.0},
                "stress_points_end_a": cbar_point_a,
                "stress_points_end_b": cbar_point_b,
            },
            30: _shell_results("CTRIA3"),
            40: _shell_results("CQUAD4"),
        },
    }

    report_path = tmp_path / "nested" / "mixed_model.rpt"
    write_report(report_path, "mixed_model", model, results, solver)
    report = report_path.read_text(encoding="utf-8")

    assert "PRIMARY 1" not in report
    assert "FORCES IN CROD ELEMENTS" in report
    assert "END FORCES AND MOMENTS IN CBAR ELEMENTS" in report
    assert "FORCE AND MOMENT RESULTANTS IN CTRIA3 ELEMENTS" in report
    assert "TOP AND BOTTOM FIBER STRESSES IN CTRIA3 ELEMENTS" in report
    assert "FORCE AND MOMENT RESULTANTS IN CQUAD4 ELEMENTS" in report
    assert "TOP AND BOTTOM FIBER STRESSES IN CQUAD4 ELEMENTS" in report
    assert "MEMBRANE NX" in report
    assert "BENDING MX" in report
    assert "SHEAR QX" in report
    assert "NORMAL X" in report
    assert "MAJOR PRINCIPAL" in report
    assert report.count("TOP (+T/2)") == 2
    assert report.count("BOTTOM (-T/2)") == 2
    assert report.index("FORCE AND MOMENT RESULTANTS IN CTRIA3") < report.index(
        "TOP AND BOTTOM FIBER STRESSES IN CTRIA3"
    )

    lines = report.splitlines()
    shell_header_index = next(
        index
        for index, line in enumerate(lines)
        if "MEMBRANE NX" in line and "MEMBRANE NY" in line
    )
    shell_header = lines[shell_header_index]
    shell_data = lines[shell_header_index + 2]
    assert shell_header.index("ELEMENT ID") + len("ELEMENT ID") == shell_data.index("30") + len("30")
    assert shell_header.index("MEMBRANE NX") + len("MEMBRANE NX") == shell_data.index("1.000000E+00") + len("1.000000E+00")

    fiber_header_index = next(
        index
        for index, line in enumerate(lines)
        if "FIBER" in line and "NORMAL X" in line
    )
    fiber_header = lines[fiber_header_index]
    fiber_data = lines[fiber_header_index + 2]
    assert fiber_header.index("FIBER") == fiber_data.index("TOP (+T/2)")


def test_report_handles_models_without_reactions_or_element_results(tmp_path):
    model = SimpleNamespace(
        nodes={1: object()}, materials={}, properties={}, elements={},
        loads=[], constraints=[],
    )
    solver = SimpleNamespace(dof_mapper=SequentialDofMapper())
    results = {
        "displacements": np.zeros(6),
        "reactions": {},
        "stresses": {},
    }

    report_path = tmp_path / "empty_results.rpt"
    write_report(report_path, "empty_results", model, results, solver)
    report = report_path.read_text(encoding="utf-8")

    assert "(No non-zero SPC reactions recorded)" in report
    assert "(No element results recorded)" in report
