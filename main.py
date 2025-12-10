"""
@file main.py
@brief Entry point for the FEM Solver application.
@details
Usage:
    python main.py --file input/Ejemplo.dat   -> Runs a specific file
    python main.py --all                      -> Runs all .dat files in /input directory
"""
import sys
import os
import argparse
import numpy as np
from pathlib import Path
from time import time


FILE_PATH = Path(__file__).resolve()  # Ruta absoluta de este archivo main.py
PROJECT_ROOT = FILE_PATH.parent       # Carpeta raíz del proyecto

'''# Añadimos al path solo si no está ya
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))'''

# --- IMPORTS DEL PROYECTO (Deben ir DESPUÉS del sys.path.insert) ---
try:
    from src.preprocessing.nastran_reader import parse_nastran_input
    from src.solver.linear_static_solver import LinearStaticSolver
    from src.preprocessing.data_model import (
        SinglePointLoad, SinglePointMoment, SPCConstraint, RodElement
    )
except ImportError as e:
    print("\n[ERROR CRÍTICO] No se pudieron importar los módulos del proyecto.")
    print(f"Detalle del error: {e}")
    sys.exit(1)



# --- CONFIGURATION ---
INPUT_DIR = PROJECT_ROOT / "input"
OUTPUT_DIR = PROJECT_ROOT / "output"

def main():
    # 1. Setup Arguments
    parser = argparse.ArgumentParser(description="PyFEM Solver - Linear Static Analysis")
    parser.add_argument('--file', type=str, help="Path to a specific .dat file")
    parser.add_argument('--all', action='store_true', help="Run all .dat files in the input folder")
    
    args = parser.parse_args()

    # 2. Ensure directories exist
    if not INPUT_DIR.exists():
        INPUT_DIR.mkdir()
        print(f"[INFO] Created '{INPUT_DIR}' folder. Place your .dat files here.")
        return
    
    if not OUTPUT_DIR.exists():
        OUTPUT_DIR.mkdir()

    # 3. Job Dispatcher
    if args.file:
        file_path = Path(args.file)
        run_analysis_job(file_path)
    elif args.all:
        print(f"[INFO] Batch processing all .dat files in {INPUT_DIR}...")
        files = list(INPUT_DIR.glob("*.dat"))
        if not files:
            print("[WARNING] No .dat files found in input directory.")
        for f in files:
            run_analysis_job(f)
    else:
        # --- MODIFICACIÓN: Comportamiento por defecto ---
        print("[INFO] No se indicaron argumentos. Buscando en carpeta 'input'...")
        
        # Buscar todos los .dat
        files = list(INPUT_DIR.glob("*.dat"))
        
        if len(files) > 0:
            # Si hay archivos, cogemos el primero (o el que quieras)
            archivo_a_correr = files[0] 
            print(f"[AUTO] Se ha encontrado: '{archivo_a_correr.name}'. Ejecutando...")
            run_analysis_job(archivo_a_correr)
        else:
            print("[ERROR] No hay archivos .dat en la carpeta input/ para ejecutar.")

def run_analysis_job(file_path: Path):
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
        model = parse_nastran_input(str(file_path))
        print(f"   -> Nodes: {len(model.nodes)}")
        print(f"   -> Elements: {len(model.elements)}")
        print(f"   -> Loads: {len(model.loads)}")
        
        # --- STEP 2: SOLVER INITIALIZATION ---
        print(f"[{time()-t0:.4f}s] Initializing Solver...")
        solver = LinearStaticSolver(model)
        
        # --- STEP 3: ASSEMBLY ---
        print(f"[{time()-t0:.4f}s] Assembling System Matrices...")
        solver.assemble()
        
        # --- STEP 4: SOLUTION ---
        print(f"[{time()-t0:.4f}s] Solving Linear System...")
        results = solver.solve() # Returns dict {displacements, reactions, stresses}
        
        # --- STEP 5: POST-PROCESSING / REPORTING ---
        print(f"[{time()-t0:.4f}s] Generating Report...")
        report_path = OUTPUT_DIR / f"{job_name}.rpt"
        write_report(report_path, job_name, model, results, solver)
        
        print(f"[{time()-t0:.4f}s] JOB COMPLETED SUCCESSFULLY.")
        print(f"Result file: {report_path}")

    except Exception as e:
        print(f"\n[ERROR] Job Failed: {e}")
        # import traceback; traceback.print_exc() # Uncomment for deep debugging

def write_report(path: Path, job_name: str, model, results, solver):
    """
    @brief Writes a formatted text file with the analysis results.
    @details Uses the DofMapper inside 'solver' to translate indices back to Node IDs.
    """
    displacements = results['displacements']
    reactions = results['reactions']
    stresses = results['stresses']
    
    # --- CONTEO PARA EL RESUMEN ---
    count_grid = len(model.nodes)
    count_mat1 = len(model.materials)
    count_prod = len(model.properties)

    # Contar elementos por tipo
    count_crod = sum(1 for l in model.loads if isinstance(l, SinglePointLoad))
    
    # Contar Cargas por tipo
    count_force = sum(1 for l in model.loads if isinstance(l, SinglePointLoad))
    count_moment = sum(1 for l in model.loads if isinstance(l, SinglePointMoment))
    
    # Contar Constraints (SPC1)
    count_spc = sum(1 for c in model.constraints if isinstance(c, SPCConstraint))
    
    with open(path, 'w') as f:
        f.write(f"ANALYSIS REPORT: {job_name}\n")
        f.write("="*60 + "\n\n")

       # --- 0. MODEL SUMMARY ---
        f.write("                       M O D E L   S U M M A R Y\n")
        f.write("\n")
        f.write(f"      {'ENTRY NAME':<15} {'NUMBER OF ENTRIES':>20}\n")
        f.write(f"      {'-'*15:<15} {'-'*20:>20}\n")
        
        if count_grid > 0:   f.write(f"      {'GRID':<15} {count_grid:>20}\n")
        if count_crod > 0:   f.write(f"      {'CROD':<15} {count_crod:>20}\n")
        if count_mat1 > 0:   f.write(f"      {'MAT1':<15} {count_mat1:>20}\n")
        if count_prod > 0:   f.write(f"      {'PROD':<15} {count_prod:>20}\n")
        if count_force > 0:  f.write(f"      {'FORCE':<15} {count_force:>20}\n")
        if count_moment > 0: f.write(f"      {'MOMENT':<15} {count_moment:>20}\n")
        if count_spc > 0:    f.write(f"      {'SPC1':<15} {count_spc:>20}\n")
        
        f.write("\n" + "="*80 + "\n\n")
   
        # --- 1. DISPLACEMENTS ---
        f.write("1. NODAL DISPLACEMENTS\n")
        f.write("-" * 85 + "\n")
        f.write(f"{'NODE ID':<10} {'T1 (X)':<12} {'T2 (Y)':<12} {'T3 (Z)':<12} {'R1 (Rx)':<12} {'R2 (Ry)':<12} {'R3 (Rz)':<12}\n")
        f.write("-" * 85 + "\n")

        # Iterate over sorted node IDs for clean output
        for nid in sorted(model.nodes.keys()):
            try:
                # Use Mapper to get the 6 indices for this node
                indices = solver.mapper.get_dof_indices(nid)
                u = displacements[indices] # Extract 6 components
                
                f.write(f"{nid:<10} {u[0]:.4e}   {u[1]:.4e}   {u[2]:.4e}   {u[3]:.4e}   {u[4]:.4e}   {u[5]:.4e}\n")
            except Exception:
                f.write(f"{nid:<10} {'--- ERROR RETRIEVING DATA ---'}\n")
        f.write("\n")

        # --- 2. REACTION FORCES (Only constrained nodes) ---
        f.write("2. SPC REACTION FORCES (Non-zero only)\n")
        f.write("-" * 85 + "\n")
        f.write(f"{'NODE ID':<10} {'F1 (X)':<12} {'F2 (Y)':<12} {'F3 (Z)':<12} {'M1 (Mx)':<12} {'M2 (My)':<12} {'M3 (Mz)':<12}\n")
        f.write("-" * 85 + "\n")

        if not reactions:
            f.write("   (No reactions recorded)\n")
        else:
            for nid in sorted(reactions.keys()):
                r = reactions[nid]
              # Filter to show only if magnitude is significant (> 1e-9)
                if np.linalg.norm(r) > 1e-9:
                    f.write(f"{nid:<10} {r[0]:.4e}   {r[1]:.4e}   {r[2]:.4e}   {r[3]:.4e}   {r[4]:.4e}   {r[5]:.4e}\n")
        f.write("\n")

        # --- 3. ELEMENT FORCES & STRESSES ---
        f.write("3. ELEMENT FORCES AND STRESSES (CROD)\n")
        f.write("-" * 60 + "\n")
        f.write(f"{'ELEM':<8} {'AXIAL FORCE':<15} {'AXIAL STRESS':<15} {'TORQUE':<15} {'TORSIONAL STRESS':<18}\n")
        f.write("-" * 60 + "\n")

        for eid in sorted(stresses.keys()):
            data = stresses[eid]
            # Handle potential missing keys gracefully
            ax_f = data.get('axial_force', 0.0)
            ax_s = data.get('axial_stress', 0.0)
            torq = data.get('torque', 0.0)
            tau  = data.get('torsional_stress', 0.0)
            
            # Formateamos con la nueva columna
            f.write(f"{eid:<8} {ax_f:<15.4e} {ax_s:<15.4e} {torq:<15.4e} {tau:<18.4e}\n")
            

if __name__ == "__main__":
    main()



