"""Compare both CQUAD4 formulations with the stored Nastran KAAX DMIGs.

The script reads the actual validation decks, builds one element with either
the Nastran reconstruction or classical Reissner-Mindlin formulation, and
reports relative Frobenius and maximum absolute errors on the active block.

Run from the repository root with::

    python tests/solver_tst/elements/run_cquad4_nastran_comparison.py

Use ``--export-dir path`` to also write each calculated 24x24 matrix as CSV.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.preprocessing.nastran_reader_v2 import parse_nastran_input
from src.solver.elements.cquad4_nastran_element import CQUAD4NastranElementSolver
from src.solver.elements.cquad4_reissner_mindlin_element import (
    CQUAD4ReissnerMindlinElementSolver,
)


VALIDATION_ROOT = Path(__file__).parents[1] / "validation_cases"
CASES = {
    "51": (
        "case_51",
        "ej51_kmembrane_k6rot_0.dat",
        "ej51_kmembrane_k6rot_0_Nastran.csv",
        [0, 1, 6, 7, 12, 13, 18, 19],
    ),
    "52": (
        "case_52",
        "ej52_k_pure_flex_mid3_0_k6rot_0.dat",
        "ej52_k_pure_flex_mid3_0_k6rot_0_Nastran.csv",
        [2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22],
    ),
    "53": (
        "case_53",
        "ej53_k_memb-flex_mid3_0_k6rot_0.dat",
        "ej53_k_memb-flex_mid3_0_k6rot_0_Nastran.csv",
        [
            0, 1, 2, 3, 4,          # Nudo 1: u, v, w, th_x, th_y (th_z inactivo)
            6, 7, 8, 9, 10,         # Nudo 2: u, v, w, th_x, th_y
            12, 13, 14, 15, 16,     # Nudo 3: u, v, w, th_x, th_y
            18, 19, 20, 21, 22      # Nudo 4: u, v, w, th_x, th_y
        ],
    ),
    "54": (
        "case_54",
        "ej54_k_memb-flex-3d_mid3_0_k6rot_0.dat",
        "ej54_k_memb-flex-3d_mid3_0_k6rot_0_Nastran.csv",
        # Inclined element: its local null drilling direction is not aligned
        # with a single global rotational component, hence all DOFs are kept.
        list(range(24)),
    ),
}


def calculate(case: str, formulation: str) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """Return calculated/reference matrices and active DOFs for one case."""
    folder, deck, reference, active_dofs = CASES[case]
    case_path = VALIDATION_ROOT / folder
    model = parse_nastran_input(str(case_path / deck))
    element = next(iter(model.elements.values()))
    nodes = [model.get_node(node_id) for node_id in element.nodes]
    prop = model.properties[element.prop_id]

    solver_class = {
        "nastran": CQUAD4NastranElementSolver,
        "reissner-mindlin": CQUAD4ReissnerMindlinElementSolver,
    }[formulation]
    solver = solver_class(
        *nodes,
        prop,
        mat_membrane=model.materials[prop.mat_id],
        mat_bending=model.materials.get(prop.mat_id_bending),
        mat_shear=model.materials.get(prop.mat_id_shear),
        k6rot=0.0,
        theta=element.theta,
    )
    calculated = solver.get_stiffness_matrix()
    nastran = pd.read_csv(case_path / reference).iloc[:, 2:].to_numpy(float)
    return calculated, nastran, active_dofs


def main() -> None:
    """Run selected comparisons and optionally export calculated matrices."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES))
    parser.add_argument("--export-dir", type=Path)
    args = parser.parse_args()

    print("case  formulation  relative Frobenius error  maximum absolute error")
    for case in args.cases:
        for formulation in ("nastran", "reissner-mindlin"):
            calculated, reference, active = calculate(case, formulation)
            index = np.ix_(active, active)
            difference = calculated[index] - reference[index]
            relative_error = np.linalg.norm(difference) / np.linalg.norm(reference[index])
            print(
                f"{case:>4}  {formulation:<11} {relative_error:>23.9%}"
                f"  {np.max(np.abs(difference)):>22.9g}"
            )

            if args.export_dir is not None:
                args.export_dir.mkdir(parents=True, exist_ok=True)
                np.savetxt(
                    args.export_dir / f"case_{case}_{formulation}.csv",
                    calculated,
                    delimiter=",",
                )


if __name__ == "__main__":
    main()
