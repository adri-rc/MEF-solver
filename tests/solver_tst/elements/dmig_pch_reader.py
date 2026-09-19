"""Minimal reader for symmetric DMIG matrices in Nastran punch files."""

from pathlib import Path

import numpy as np


def _parse_nastran_float(value: str) -> float:
    """Parse standard and compact Nastran floating-point notation."""
    value = value.strip().upper().replace("D", "E")
    if "E" not in value:
        exponent_pos = max(value.rfind("-", 1), value.rfind("+", 1))
        if exponent_pos != -1:
            value = f"{value[:exponent_pos]}E{value[exponent_pos:]}"
    return float(value)


def read_dmig_pch(
    file_path: str | Path,
    matrix_name: str = "KAAX",
    node_count: int = 3,
) -> np.ndarray:
    """Read a symmetric six-DOF-per-node DMIG matrix from a punch file.

    The validation decks number the element nodes consecutively from one, so
    each ``(GRID, component)`` pair maps directly into the returned matrix.
    """
    stiffness = np.zeros((node_count * 6, node_count * 6))
    in_matrix = False
    column_node = None
    column_dof = None

    with Path(file_path).open(encoding="latin-1") as punch_file:
        for line in punch_file:
            if line.startswith("DMIG "):
                fields = line.split()
                in_matrix = len(fields) > 1 and fields[1] == matrix_name
            elif in_matrix and line.startswith("DMIG*"):
                column_node = int(line[32:48].strip())
                column_dof = int(line[48:64].strip())
            elif in_matrix and line.startswith("*"):
                row_node = int(line[8:24].strip())
                row_dof = int(line[24:40].strip())
                value = _parse_nastran_float(line[40:56])

                column_index = (column_node - 1) * 6 + column_dof - 1
                row_index = (row_node - 1) * 6 + row_dof - 1
                stiffness[row_index, column_index] = value
                stiffness[column_index, row_index] = value

    return stiffness
