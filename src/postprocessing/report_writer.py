"""Formatted text reports for linear-static analyses.

The report layout lives here so the application entry point only orchestrates
parse, solve and post-processing.  Element results are intentionally split by
element family: a single set of generic columns cannot describe rods, bars and
shell resultants without losing their physical meaning.
"""

from collections import defaultdict
from numbers import Integral, Real
from pathlib import Path
from typing import Iterable, Mapping, TextIO

import numpy as np

from src.preprocessing.data_model import (
    BarElement,
    BarProperty,
    CQUAD4Element,
    CTRIA3Element,
    RodElement,
    RodProperty,
    ShellProperty,
    SinglePointLoad,
    SinglePointMoment,
    SPCConstraint,
)


_REPORT_WIDTH = 148
_RESULT_TYPE_ORDER = ("CROD", "CBAR", "CTRIA3", "CQUAD4")


def _format_cell(value, width: int) -> str:
    if isinstance(value, Integral):
        return f"{int(value):>{width}d}"
    if isinstance(value, Real):
        return f"{float(value):>{width}.6E}"
    return f"{str(value):<{width}}"


def _write_table(
    stream: TextIO,
    headers: tuple[str, ...],
    widths: tuple[int, ...],
    rows: Iterable[tuple],
    left_aligned_columns: frozenset[int] = frozenset(),
) -> None:
    """Write a fixed-width table using the same alignment for headers and data."""
    if len(headers) != len(widths):
        raise ValueError("Each report column requires a matching width.")

    line_width = sum(widths) + len(widths) - 1
    header_cells = (
        f"{header:<{width}}" if index in left_aligned_columns else f"{header:>{width}}"
        for index, (header, width) in enumerate(zip(headers, widths))
    )
    stream.write(" ".join(header_cells) + "\n")
    stream.write("-" * line_width + "\n")
    for row in rows:
        if len(row) != len(headers):
            raise ValueError("Report row length does not match its table header.")
        stream.write(" ".join(_format_cell(value, width) for value, width in zip(row, widths)) + "\n")
    stream.write("\n")


def _write_report_heading(stream: TextIO, job_name: str) -> None:
    stream.write("LINEAR STATIC ANALYSIS REPORT\n")
    stream.write(f"JOB: {job_name}\n")
    stream.write("All values use the consistent unit system defined by the input model.\n")
    stream.write("=" * _REPORT_WIDTH + "\n\n")


def _write_model_summary(stream: TextIO, model) -> None:
    property_counts = {
        "PROD": sum(isinstance(prop, RodProperty) for prop in model.properties.values()),
        "PBAR": sum(isinstance(prop, BarProperty) for prop in model.properties.values()),
        "PSHELL": sum(isinstance(prop, ShellProperty) for prop in model.properties.values()),
    }
    element_counts = {
        "CROD": sum(isinstance(element, RodElement) for element in model.elements.values()),
        "CBAR": sum(isinstance(element, BarElement) for element in model.elements.values()),
        "CTRIA3": sum(isinstance(element, CTRIA3Element) for element in model.elements.values()),
        "CQUAD4": sum(isinstance(element, CQUAD4Element) for element in model.elements.values()),
    }
    counts = {
        "GRID": len(model.nodes),
        "MAT1": len(model.materials),
        **property_counts,
        **element_counts,
        "FORCE": sum(isinstance(load, SinglePointLoad) for load in model.loads),
        "MOMENT": sum(isinstance(load, SinglePointMoment) for load in model.loads),
        "SPC": sum(isinstance(constraint, SPCConstraint) for constraint in model.constraints),
    }

    stream.write("0. MODEL SUMMARY\n")
    rows = ((name, count) for name, count in counts.items() if count)
    _write_table(
        stream,
        ("ENTRY TYPE", "NUMBER OF ENTRIES"),
        (22, 20),
        rows,
        left_aligned_columns=frozenset({0}),
    )


def _write_displacements(stream: TextIO, model, displacements, dof_mapper) -> None:
    stream.write("1. NODAL DISPLACEMENTS\n")
    headers = ("NODE ID", "T1 (X)", "T2 (Y)", "T3 (Z)", "R1 (RX)", "R2 (RY)", "R3 (RZ)")
    widths = (10, 17, 17, 17, 17, 17, 17)

    def rows():
        for node_id in sorted(model.nodes):
            indices = dof_mapper.get_dof_indices(node_id)
            values = displacements[indices]
            yield (node_id, *values)

    _write_table(stream, headers, widths, rows())


def _write_reactions(stream: TextIO, reactions: Mapping[int, np.ndarray]) -> None:
    stream.write("2. SPC REACTION FORCES AND MOMENTS (NON-ZERO NODES)\n")
    headers = ("NODE ID", "F1 (X)", "F2 (Y)", "F3 (Z)", "M1 (MX)", "M2 (MY)", "M3 (MZ)")
    widths = (10, 17, 17, 17, 17, 17, 17)
    rows = [
        (node_id, *reactions[node_id])
        for node_id in sorted(reactions)
        if np.linalg.norm(reactions[node_id]) > 1.0e-9
    ]
    if rows:
        _write_table(stream, headers, widths, rows)
    else:
        stream.write("(No non-zero SPC reactions recorded)\n\n")


def _write_crod_results(stream: TextIO, section: int, items: list[tuple[int, Mapping]]) -> int:
    stream.write(f"3.{section}. FORCES IN CROD ELEMENTS\n")
    _write_table(
        stream,
        ("ELEMENT ID", "AXIAL FORCE", "TORQUE"),
        (12, 18, 18),
        ((eid, data.get("axial_force", 0.0), data.get("torque", 0.0)) for eid, data in items),
    )
    section += 1

    stream.write(f"3.{section}. STRESSES AND STRAINS IN CROD ELEMENTS\n")
    _write_table(
        stream,
        ("ELEMENT ID", "AXIAL STRESS", "TORSIONAL STRESS", "AXIAL STRAIN", "TWIST RATE"),
        (12, 18, 20, 18, 18),
        (
            (
                eid,
                data.get("axial_stress", 0.0),
                data.get("torsional_stress", 0.0),
                data.get("axial_strain", 0.0),
                data.get("twist_rate", 0.0),
            )
            for eid, data in items
        ),
    )
    return section + 1


def _write_cbar_results(stream: TextIO, section: int, items: list[tuple[int, Mapping]]) -> int:
    stream.write(f"3.{section}. END FORCES AND MOMENTS IN CBAR ELEMENTS\n")
    force_rows = []
    for eid, data in items:
        for end, key in (("A", "forces_n1"), ("B", "forces_n2")):
            forces = data.get(key, {})
            force_rows.append(
                (
                    eid, end,
                    forces.get("Fx", 0.0), forces.get("Fy", 0.0), forces.get("Fz", 0.0),
                    forces.get("Mx", 0.0), forces.get("My", 0.0), forces.get("Mz", 0.0),
                )
            )
    _write_table(
        stream,
        ("ELEMENT ID", "END", "AXIAL FX", "SHEAR FY", "SHEAR FZ", "TORQUE MX", "BENDING MY", "BENDING MZ"),
        (12, 6, 17, 17, 17, 17, 17, 17),
        force_rows,
        left_aligned_columns=frozenset({1}),
    )
    section += 1

    stream.write(f"3.{section}. STRESSES AT RECOVERY POINTS IN CBAR ELEMENTS\n")
    stress_rows = []
    for eid, data in items:
        for end, key in (("A", "stress_points_end_a"), ("B", "stress_points_end_b")):
            for point_id, point in sorted(data.get(key, {}).items()):
                stress_rows.append(
                    (
                        eid, end, point_id, point.get("y", 0.0), point.get("z", 0.0),
                        point.get("sigma_axial", 0.0), point.get("sigma_bending", 0.0),
                        point.get("sigma_total", 0.0), point.get("tau_xy", 0.0),
                    )
                )
    _write_table(
        stream,
        ("ELEMENT ID", "END", "POINT", "LOCAL Y", "LOCAL Z", "AXIAL", "BENDING", "TOTAL", "TORSIONAL"),
        (12, 6, 8, 16, 16, 17, 17, 17, 17),
        stress_rows,
        left_aligned_columns=frozenset({1, 2}),
    )
    return section + 1


def _write_shell_results(
    stream: TextIO,
    section: int,
    element_type: str,
    items: list[tuple[int, Mapping]],
) -> int:
    stream.write(f"3.{section}. FORCE AND MOMENT RESULTANTS IN {element_type} ELEMENTS\n")
    stream.write("     N and Q are force resultants per unit length; M are moment resultants per unit length.\n")
    _write_table(
        stream,
        ("ELEMENT ID", "MEMBRANE NX", "MEMBRANE NY", "MEMBRANE NXY", "BENDING MX", "BENDING MY", "TWIST MXY", "SHEAR QX", "SHEAR QY"),
        (12, 17, 17, 18, 17, 17, 17, 17, 17),
        (
            (
                eid,
                data.get("membrane_force_x", 0.0), data.get("membrane_force_y", 0.0),
                data.get("membrane_force_xy", 0.0), data.get("bending_moment_x", 0.0),
                data.get("bending_moment_y", 0.0), data.get("bending_moment_xy", 0.0),
                data.get("shear_force_x", 0.0), data.get("shear_force_y", 0.0),
            )
            for eid, data in items
        ),
    )
    section += 1

    stream.write(f"3.{section}. TOP AND BOTTOM FIBER STRESSES IN {element_type} ELEMENTS\n")
    stream.write("     TOP is local z=+t/2; BOTTOM is local z=-t/2.\n")
    stress_rows = []
    for eid, data in items:
        stress_rows.extend(
            (
                (
                    eid, "TOP (+T/2)", data.get("stress_top_fiber_x", 0.0),
                    data.get("stress_top_fiber_y", 0.0), data.get("tau_top_fiber_xy", 0.0),
                    data.get("major_principal_top", 0.0), data.get("minor_principal_top", 0.0),
                    data.get("von_mises_top", 0.0),
                ),
                (
                    eid, "BOTTOM (-T/2)", data.get("stress_bot_fiber_x", 0.0),
                    data.get("stress_bot_fiber_y", 0.0), data.get("tau_bot_fiber_xy", 0.0),
                    data.get("major_principal_bot", 0.0), data.get("minor_principal_bot", 0.0),
                    data.get("von_mises_bot", 0.0),
                ),
            )
        )
    _write_table(
        stream,
        ("ELEMENT ID", "FIBER", "NORMAL X", "NORMAL Y", "SHEAR XY", "MAJOR PRINCIPAL", "MINOR PRINCIPAL", "VON MISES"),
        (12, 15, 17, 17, 17, 20, 20, 17),
        stress_rows,
        left_aligned_columns=frozenset({1}),
    )
    return section + 1


def _write_element_results(stream: TextIO, stresses: Mapping[int, Mapping]) -> None:
    stream.write("3. ELEMENT RESULTS IN LOCAL ELEMENT AXES\n")
    stream.write("=" * _REPORT_WIDTH + "\n\n")
    if not stresses:
        stream.write("(No element results recorded)\n\n")
        return

    by_type: dict[str, list[tuple[int, Mapping]]] = defaultdict(list)
    for eid, data in sorted(stresses.items()):
        by_type[str(data.get("type", "UNKNOWN")).upper()].append((eid, data))

    section = 1
    for element_type in _RESULT_TYPE_ORDER:
        items = by_type.pop(element_type, [])
        if not items:
            continue
        if element_type == "CROD":
            section = _write_crod_results(stream, section, items)
        elif element_type == "CBAR":
            section = _write_cbar_results(stream, section, items)
        else:
            section = _write_shell_results(stream, section, element_type, items)

    if by_type:
        stream.write(f"3.{section}. UNRECOGNIZED ELEMENT RESULT TYPES\n")
        rows = (
            (eid, element_type)
            for element_type, items in sorted(by_type.items())
            for eid, _ in items
        )
        _write_table(
            stream,
            ("ELEMENT ID", "RESULT TYPE"),
            (12, 20),
            rows,
            left_aligned_columns=frozenset({1}),
        )


def write_report(path: Path | str, job_name: str, model, results: Mapping, solver) -> None:
    """Write the complete analysis report to ``path``.

    ``solver.dof_mapper`` is used only to map the displacement vector back to
    node IDs.  All other content comes from the model and result dictionaries.
    """
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    with report_path.open("w", encoding="utf-8", newline="\n") as stream:
        _write_report_heading(stream, job_name)
        _write_model_summary(stream, model)
        _write_displacements(stream, model, results["displacements"], solver.dof_mapper)
        _write_reactions(stream, results["reactions"])
        _write_element_results(stream, results["stresses"])
