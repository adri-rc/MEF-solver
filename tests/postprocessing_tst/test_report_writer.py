from types import SimpleNamespace

import numpy as np

from src.preprocessing.data_model import (
    AnalysisConfiguration,
    AnalysisType,
    NonlinearStepControl,
)
from src.postprocessing.report_writer import write_report
from src.solver.nonlinear_analysis import (
    NonlinearIncrementResult,
    NonlinearIterationResult,
)


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
        source_entry_counts={
            "CROD": 1,
            "FORCE": 2,
            "GRID": 2,
            "MAT1": 1,
            "MOMENT": 1,
            "SPC": 1,
            "SPC1": 2,
        },
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
                "axial_force": 101.0,
                "torque": 25.0,
                "axial_stress": 30.0,
                "torsional_stress": 5.0,
                "axial_strain": 1.5e-3,
                "twist_rate": 2.5e-3,
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

    assert report.startswith("LINEAR STATIC ANALYSIS REPORT")
    assert "0A. ANALYSIS CONFIGURATION" in report
    assert "LINEAR STATIC" in report
    summary = report[report.index("0. MODEL SUMMARY"):report.index("0A. ANALYSIS CONFIGURATION")]
    summary_counts = {
        fields[0]: int(fields[1])
        for line in summary.splitlines()
        if len(fields := line.split()) == 2 and fields[1].isdigit()
    }
    assert summary_counts["SPC"] == 1
    assert summary_counts["SPC1"] == 2
    assert summary_counts["FORCE"] == 2
    assert summary_counts["MOMENT"] == 1
    assert "T1 (X)" not in report
    assert "R1 (RX)" not in report
    assert "F1 (X)" not in report
    assert "M1 (MX)" not in report
    assert "T1" in report
    assert "R1" in report
    assert "F1" in report
    assert "M1" in report
    assert "PRIMARY 1" not in report
    assert "FORCES IN CROD ELEMENTS" in report
    assert "SECTION RESPONSE IN CBAR ELEMENTS" in report
    assert "END FORCES AND MOMENTS IN CBAR ELEMENTS" in report
    assert "NONLINEAR STATE IN CBAR ELEMENTS" not in report
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


def test_nonlinear_cbar_report_includes_section_response_and_converged_state(tmp_path):
    model = SimpleNamespace(
        nodes={1: object(), 2: object()},
        materials={},
        properties={},
        elements={},
        loads=[],
        constraints=[],
        analysis=AnalysisConfiguration(
            solution_sequence=400,
            analysis_type=AnalysisType.NONLINEAR_STATIC,
            analysis_keyword="NLSTAT",
            large_displacement_parameter=1,
            large_displacements=True,
            follower_force_stiffness=False,
        ),
    )
    solver = SimpleNamespace(dof_mapper=SequentialDofMapper())
    results = {
        "displacements": np.zeros(12),
        "reactions": {},
        "stresses": {
            21: {
                "type": "CBAR",
                "axial_force": 1.25e5,
                "torque": 2.50e4,
                "axial_stress": 250.0,
                "torsional_stress": 12.5,
                "axial_strain": 1.25e-3,
                "twist_rate": 2.50e-4,
                "current_length": 1.0025e3,
                "strain_energy": 3.75e3,
                "max_element_stress": 310.0,
                "min_element_stress": -85.0,
                "forces_n1": {
                    "Fx": -1.25e5, "Fy": 1.0e4, "Fz": -2.0e4,
                    "Mx": -2.50e4, "My": 3.0e4, "Mz": -4.0e4,
                },
                "forces_n2": {
                    "Fx": 1.25e5, "Fy": -1.0e4, "Fz": 2.0e4,
                    "Mx": 2.50e4, "My": -3.0e4, "Mz": 4.0e4,
                },
                "stress_points_end_a": {
                    "P0": {
                        "y": 10.0, "z": 20.0, "sigma_axial": 250.0,
                        "sigma_bending": 60.0, "sigma_total": 310.0,
                        "tau_xy": 12.5,
                    }
                },
                "stress_points_end_b": {
                    "P0": {
                        "y": 10.0, "z": 20.0, "sigma_axial": 250.0,
                        "sigma_bending": -335.0, "sigma_total": -85.0,
                        "tau_xy": 12.5,
                    }
                },
            }
        },
    }

    report_path = tmp_path / "nonlinear_cbar.rpt"
    write_report(report_path, "nonlinear_cbar", model, results, solver)
    report = report_path.read_text(encoding="utf-8")

    assert report.startswith("NONLINEAR STATIC ANALYSIS REPORT")
    assert "SECTION RESPONSE IN CBAR ELEMENTS" in report
    assert "AXIAL FORCE" in report
    assert "MAX TORSIONAL STRESS" in report
    assert "AXIAL STRAIN" in report
    assert "TWIST RATE" in report
    assert "1.250000E+05" in report
    assert "2.500000E-04" in report

    assert "NONLINEAR STATE IN CBAR ELEMENTS" in report
    assert "CURRENT LENGTH" in report
    assert "ELEMENT STRAIN ENERGY" in report
    assert "MAX NORMAL STRESS" in report
    assert "MIN NORMAL STRESS" in report
    assert "1.002500E+03" in report
    assert "3.750000E+03" in report
    assert "3.100000E+02" in report
    assert "-8.500000E+01" in report

    assert "END FORCES AND MOMENTS IN CBAR ELEMENTS" in report
    assert "STRESSES AT RECOVERY POINTS IN CBAR ELEMENTS" in report
    assert "P0" not in report
    assert any(
        line.split()[2] == "C"
        for line in report.splitlines()
        if line.strip().startswith("21 A")
    )
    assert report.index("SECTION RESPONSE IN CBAR ELEMENTS") < report.index(
        "NONLINEAR STATE IN CBAR ELEMENTS"
    )
    assert report.index("NONLINEAR STATE IN CBAR ELEMENTS") < report.index(
        "END FORCES AND MOMENTS IN CBAR ELEMENTS"
    )


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


def test_nonlinear_report_includes_configuration_histories_and_final_state(tmp_path):
    control = NonlinearStepControl(
        sid=127,
        source="NLSTEP",
        total_time=2.0,
        control_definition="SEVERELY",
        increment_method="ADAPT",
        number_of_increments=None,
        initial_step_fraction=0.05,
        minimum_step_fraction=0.01,
        maximum_step_fraction=0.20,
        desired_iterations=4,
        step_scale_factor=1.25,
        output_intervals=10,
        maximum_increments=100,
        maximum_iterations=25,
        minimum_iterations=1,
        maximum_bisections=8,
        convergence="PV",
        displacement_tolerance=1.0e-7,
        residual_tolerance=2.0e-7,
        work_tolerance=3.0e-7,
        stiffness_update_method="PFNT",
    )
    model = SimpleNamespace(
        nodes={1: object()},
        materials={},
        properties={},
        elements={},
        loads=[],
        constraints=[],
        analysis=AnalysisConfiguration(
            solution_sequence=400,
            analysis_type=AnalysisType.NONLINEAR_STATIC,
            analysis_keyword="NLSTAT",
            large_displacement_parameter=1,
            large_displacements=True,
            follower_force_stiffness=False,
            nonlinear_step=control,
        ),
    )
    solver = SimpleNamespace(dof_mapper=SequentialDofMapper())
    increments = [
        NonlinearIncrementResult(
            increment=1,
            load_factor=0.05,
            step_size=0.05,
            iterations=2,
            residual_norm=4.0e-9,
            relative_residual_norm=5.0e-10,
            cutbacks=1,
            displacements=np.zeros(6),
            strain_energy=12.5,
        )
    ]
    iterations = [
        NonlinearIterationResult(
            increment=1,
            load_factor=0.05,
            iteration=1,
            residual_norm=6.0e-4,
            relative_residual_norm=7.0e-5,
            correction_norm=8.0e-3,
        ),
        NonlinearIterationResult(
            increment=1,
            load_factor=0.05,
            iteration=2,
            residual_norm=4.0e-9,
            relative_residual_norm=5.0e-10,
            correction_norm=9.0e-6,
        ),
    ]
    results = {
        "displacements": np.array([1.0, 2.0, 3.0, 0.0, 0.0, 0.0]),
        "reactions": {1: np.array([10.0, 0.0, 0.0, 0.0, 0.0, 0.0])},
        "stresses": {
            10: {
                "type": "CROD",
                "axial_force": 100.0,
                "axial_stress": 50.0,
                "axial_strain": 0.01,
            }
        },
        "increments": increments,
        "iterations": iterations,
        "load_factor": 0.05,
        "strain_energy": 12.5,
        "cutbacks": 1,
    }

    report_path = tmp_path / "nonlinear.rpt"
    write_report(report_path, "nonlinear", model, results, solver)
    report = report_path.read_text(encoding="utf-8")

    assert report.startswith("NONLINEAR STATIC ANALYSIS REPORT")
    assert "0A. ANALYSIS CONFIGURATION" in report
    assert "SOL 400" in report
    assert "NLSTAT" in report
    assert "LARGE DISPLACEMENTS" in report
    assert "YES" in report
    assert "0B. NONLINEAR STEP CONTROL" in report
    assert "NLSTEP" in report
    assert "ADAPT" in report
    assert "INITIAL STEP FRACTION" in report
    assert "RESIDUAL TOLERANCE" in report
    assert "0C. NONLINEAR SOLUTION SUMMARY" in report
    assert "4.1. CONVERGED INCREMENT HISTORY" in report
    assert "STEP SIZE" in report
    assert "STRAIN ENERGY" in report
    assert "4.2. NEWTON ITERATION HISTORY" in report
    assert "CORRECTION NORM" in report
    assert "6.000000E-04" in report
    assert "8.000000E-03" in report
    assert "FINAL CONVERGED STATE RESULTS" in report
    assert report.index("FINAL CONVERGED STATE RESULTS") < report.index(
        "NODAL DISPLACEMENTS"
    )
    assert report.index("ELEMENT RESULTS") < report.index(
        "CONVERGED INCREMENT HISTORY"
    )
    assert report.index("CONVERGED INCREMENT HISTORY") < report.index(
        "NEWTON ITERATION HISTORY"
    )
