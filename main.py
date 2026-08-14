"""
@file main.py
@brief Entry point for the FEM Solver application.
@details
Usage:
    python main.py --file input/Ejemplo.dat   -> Runs a specific file
    python main.py --all                      -> Runs all .dat files in /input directory
"""
import sys
import argparse
from pathlib import Path
from time import time

FILE_PATH = Path(__file__).resolve()  
PROJECT_ROOT = FILE_PATH.parent      

try:
    from src.preprocessing.nastran_reader_v2 import parse_nastran_input
    from src.solver.linear_static_solver import LinearStaticSolver
    from src.postprocessing.report_writer import write_report
except ImportError as e:
    print("\n[CRITICAL ERROR] The project modules could not be imported.")
    print(f"Details: {e}")
    sys.exit(1)


# --- PATHS CONFIGURATION ---
INPUT_DIR = PROJECT_ROOT / "input"
OUTPUT_DIR = PROJECT_ROOT / "output"


def main():
    # 1. Setup Arguments
    parser = argparse.ArgumentParser(description="PyFEM Solver - Linear Static Analysis")
    parser.add_argument('--file', type=str, help="Path to a specific .dat file")
    parser.add_argument('--all', action='store_true', help="Run all .dat files in the input folder")
    parser.add_argument('--subcase', type=int, help="SUBCASE ID to solve when a deck contains several")
    
    args = parser.parse_args()

    if not INPUT_DIR.exists():
        INPUT_DIR.mkdir()
        print(f"[INFO] Created '{INPUT_DIR}' folder. Place your .dat files here.")
        return
    
    if not OUTPUT_DIR.exists():
        OUTPUT_DIR.mkdir()

    if args.file:
        file_path = Path(args.file)
        run_analysis_job(file_path, subcase_id=args.subcase)
    elif args.all:
        print(f"[INFO] Batch processing all .dat files in {INPUT_DIR}...")
        files = list(INPUT_DIR.glob("*.dat"))
        if not files:
            print("[WARNING] No .dat files found in input directory.")
        for f in files:
            run_analysis_job(f, subcase_id=args.subcase)
    else:
        print("[INFO] No arguments were specified. Searching in /input...")

        files = list(INPUT_DIR.glob("*.dat"))
        
        if len(files) > 0:
            file_to_run = files[0]
            run_analysis_job(file_to_run, subcase_id=args.subcase)
        else:
            print("[ERROR] No .dat found in input/")

def run_analysis_job(file_path: Path, subcase_id: int | None = None):
    """
    @brief Orchestrates the analysis for a single file: Parse -> Solve -> Report.
    """
    job_name = file_path.stem
    print("="*80)
    print(f"STARTING JOB: {job_name}")
    print(f"File: {file_path}")
    print("="*80)

    t0 = time()

    try:
        # --- STEP 1: PRE-PROCESSING ---
        print(f"[{time()-t0:.4f}s] Parsing input file...")
        model = parse_nastran_input(str(file_path), subcase_id=subcase_id)
        
        # --- STEP 2: SOLVER INITIALIZATION ---
        print(f"[{time()-t0:.4f}s] Initializing Solver...")
        solver = LinearStaticSolver(model)
        
        # --- STEP 3: ASSEMBLY ---
        print(f"[{time()-t0:.4f}s] Assembling System Matrices...")
        solver.assemble()
        
        # --- STEP 4: SOLUTION ---
        print(f"[{time()-t0:.4f}s] Solving Linear System...")
        results = solver.solve()
        
        # --- STEP 5: POST-PROCESSING / REPORTING ---
        print(f"[{time()-t0:.4f}s] Generating Report...")
        report_path = OUTPUT_DIR / f"{job_name}.rpt"
        write_report(report_path, job_name, model, results, solver)
        
        print(f"[{time()-t0:.4f}s] JOB COMPLETED SUCCESSFULLY.")
        print(f"Result file: {report_path}")

    except Exception as e:
        print(f"\n[ERROR] Job Failed: {e}")
        raise

if __name__ == "__main__":
    main()



