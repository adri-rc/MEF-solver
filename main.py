"""Command-line entry point for static finite element analyses."""

import argparse
from pathlib import Path
import sys
from time import time

FILE_PATH = Path(__file__).resolve()
PROJECT_ROOT = FILE_PATH.parent

try:
    from src.preprocessing.nastran_reader import parse_nastran_input
    from src.preprocessing.data_model import AnalysisType
    from src.solver.linear_static_solver import LinearStaticSolver
    from src.solver.nonlinear_analysis import NonlinearStaticSolver
    from src.postprocessing.report_writer import write_report
except ImportError as exc:
    print("\n[CRITICAL ERROR] The project modules could not be imported.")
    print(f"Details: {exc}")
    sys.exit(1)


INPUT_DIR = PROJECT_ROOT / "input"
OUTPUT_DIR = PROJECT_ROOT / "output"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="MEF-solver static structural analysis"
    )
    parser.add_argument("--file", type=str, help="Path to a .dat or .bdf file")
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all .dat files in the input folder",
    )
    parser.add_argument(
        "--subcase",
        type=int,
        help="SUBCASE ID to solve when a deck contains several",
    )

    args = parser.parse_args()

    if not INPUT_DIR.exists():
        INPUT_DIR.mkdir()
        print(f"[INFO] Created '{INPUT_DIR}' folder. Place your .dat files here.")
    OUTPUT_DIR.mkdir(exist_ok=True)

    if args.file:
        file_path = Path(args.file)
        run_analysis_job(file_path, subcase_id=args.subcase)
    elif args.all:
        print(f"[INFO] Batch processing all .dat files in {INPUT_DIR}...")
        files = sorted(INPUT_DIR.glob("*.dat"))
        if not files:
            print("[WARNING] No .dat files found in input directory.")
        for file_path in files:
            run_analysis_job(file_path, subcase_id=args.subcase)
    else:
        print("[INFO] No arguments were specified. Searching in /input...")

        files = sorted(INPUT_DIR.glob("*.dat"))
        if files:
            run_analysis_job(files[0], subcase_id=args.subcase)
        else:
            print("[ERROR] No .dat found in input/")


def run_analysis_job(file_path: Path, subcase_id: int | None = None):
    """Parse, solve, and report one input deck."""
    job_name = file_path.stem
    print("=" * 80)
    print(f"STARTING JOB: {job_name}")
    print(f"File: {file_path}")
    print("=" * 80)

    t0 = time()

    try:
        print(f"[{time() - t0:.4f}s] Parsing input file...")
        model = parse_nastran_input(str(file_path), subcase_id=subcase_id)

        print(f"[{time() - t0:.4f}s] Initializing solver...")
        if model.analysis.analysis_type is AnalysisType.NONLINEAR_STATIC:
            solver = NonlinearStaticSolver(model)
            analysis_label = "Nonlinear Static"
        else:
            solver = LinearStaticSolver(model)
            analysis_label = "Linear Static"
        print(f"[{time() - t0:.4f}s] Assembling system matrices...")
        solver.assemble()

        print(f"[{time() - t0:.4f}s] Solving {analysis_label} system...")
        results = solver.solve()

        print(f"[{time() - t0:.4f}s] Generating report...")
        report_path = OUTPUT_DIR / f"{job_name}.rpt"
        write_report(report_path, job_name, model, results, solver)
        print(f"[{time() - t0:.4f}s] JOB COMPLETED SUCCESSFULLY.")
        print(f"Result file: {report_path}")
        return results

    except Exception as exc:
        print(f"\n[ERROR] Job Failed: {exc}")
        raise


if __name__ == "__main__":
    main()
