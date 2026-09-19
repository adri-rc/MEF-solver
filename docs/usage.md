# Using MEF-solver

## Verified environment

The repository's local virtual environment was inspected on 19 September
2026 with these versions:

| Component | Version |
|---|---:|
| Python | 3.12.14 |
| NumPy | 1.26.4 |
| SciPy | 1.13.1 |
| pandas | 2.2.3 |
| pyNastran | 1.4.1 |
| pytest | 8.3.5 |

`requirements.txt` declares minimum dependency versions. The project does not
yet provide a lock file or a `pyproject.toml`, so the table records the tested
local environment rather than a reproducible dependency contract.

On a new Windows clone, create an environment and install the declared
dependencies with:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Validation decks and Nastran reference outputs are stored through Git LFS.
After cloning, materialize them and verify that a `.dat`, `.bdf`, or `.f06`
file contains its real content rather than an LFS pointer:

```powershell
git lfs pull
```

From the repository root, either activate the environment:

```powershell
.\.venv\Scripts\Activate.ps1
python main.py --file tests\solver_tst\validation_cases\case_1\Ej1_v1.dat
```

or call its interpreter directly:

```powershell
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_1\Ej1_v1.dat
```

## Execution modes

```text
python main.py --file PATH_TO_MODEL.dat
python main.py --file PATH_TO_MODEL.dat --subcase 2
python main.py --all
python main.py
```

- `--file` analyzes the specified file. The reader accepts either a complete
  Nastran deck or a Bulk Data-only punch file.
- `--subcase ID` selects one Case Control subcase explicitly.
- `--all` processes every `*.dat` file directly under `input/`.
- With no arguments, the lexicographically first `*.dat` found directly under
  `input/` is processed.

`parse_nastran_input()` accepts `.bdf` as well as `.dat`, but automatic input
discovery only searches for `.dat`. pyNastran resolves `INCLUDE` paths relative
to the input deck.

For a linear model, `main.py` performs:

1. read and map the deck;
2. create `LinearStaticSolver`;
3. assemble `K` and `F`;
4. solve for displacements and recover reactions and element results;
5. write `output/<model-name>.rpt`.

A `SOL 400` deck with `ANALYSIS=NLSTAT` is routed to
`NonlinearStaticSolver`. That route is separate from `LinearStaticSolver` and
is limited to `MAT1` with `PROD/CROD` and/or `PBAR/PBARN1/CBAR`, conservative
nodal loads, and full Newton-Raphson. Its result dictionary additionally
contains accepted increments, Newton iterations, load factor, elastic strain
energy, and cutback count.

`run_analysis_job()` prints error context and re-raises the exception. The CLI
therefore exits with a nonzero status after an analysis error. In `--all`
mode, the first failed model stops the batch; continuing and collecting
multiple failures is not implemented.

## Supported input subset

| Group | Supported cards/parameters | Notes |
|---|---|---|
| Nodes | `GRID` | `CP` positions are converted to basic Cartesian coordinates; `CD` defines the nodal degree-of-freedom frame. |
| Materials | `MAT1` | At least two of `E`, `G`, and `NU` must be supplied; the third is calculated. |
| Properties | `PROD`, `PBAR`, `PBARN1`, `PSHELL` | `PBARN1` retains `SECT`, `BEH2`, and `INT2` for `SOL 400`. |
| Elements | `CROD`, `CBAR`, `CTRIA3`, `CQUAD4` | CTRIA3 and CQUAD4 read `12I/T^3` and `TS/T`; the active CQUAD4 is the coplanar Nastran reconstruction. |
| Loads | `FORCE`, `MOMENT`, `LOAD` | Vectors are transformed from their `CID`; nested load combinations are reduced. |
| Constraints | `SPC`, `SPC1`, `SPCADD`, `GRID.PS` | Nonzero prescribed values are supported. |
| Parameters | `PARAM,K6ROT`, `PARAM,LGDISP`, `PARAM,FOLLOWK` | `K6ROT` acts on shells; `LGDISP=-1/1/2` controls nonlinear CROD/CBAR kinematics; follower-load stiffness is outside the supported solver scope. |
| Analysis control | `SOL 101`; `SOL 400`, `ANALYSIS=NLSTAT`, `NLSTEP`, `NLPARM` | `SOL 400` is solved only for the CROD/CBAR subset above. |

The reader fails explicitly when it encounters unsupported elements, loads,
or constraints. In particular, distributed loads, MPCs, RBEs, and other rigid
elements are not ignored. Unsupported kinematic sets such as `SUPORT`,
`ASET`/`BSET`/`CSET`/`QSET`/`USET`/`OMIT`, `SPCOFF`, and superelement DOF
releases are likewise rejected instead of being dropped silently.

Only linear `SOL 101` and nonlinear-static `SOL 400` with
`ANALYSIS=NLSTAT` are accepted as complete-deck analysis controls. Other
solution sequences, `SOL 400` without `NLSTAT`, and `NLSTAT` paired with a
different solution sequence fail before assembly. Bulk Data-only punch files
remain valid and have no solution sequence.

### Subcase selection

When a complete deck has one `SUBCASE`, it is selected automatically. When it
has more than one, pass `--subcase ID`; omission is an error so that unrelated
load cases cannot be mixed accidentally. The selected `LOAD` and `SPCADD`
trees are reduced before mapping. Permanent constraints in `GRID.PS` are also
applied.

In a Bulk Data-only file without Case Control, supported base load and
constraint cards are mapped directly.

### PSHELL material semantics

For CTRIA3, `MID1`, `MID2`, and `MID3` are independent:

- blank `MID1` disables the membrane block;
- blank `MID2` disables bending;
- `MID3` cannot be supplied without `MID2`;
- blank `MID3` does not inherit another material and disables the explicit
  physical transverse-shear term in the implemented CTRIA3 formulation;
- `12I/T^3` scales bending and the effective fiber section modulus;
- `TS/T` replaces the default `5/6` transverse-shear factor.

The active CQUAD4 formulation requires `MID1`; a bending-only CQUAD4 is outside
the supported subset and is rejected explicitly. Blank `MID2` disables bending
and blank `MID3` disables physical transverse-shear flexibility without
removing the element's residual bending flexibility. Neither field inherits
another MID. Only coplanar geometry is supported.

### Nonlinear CROD/CBAR input

The reader normalizes `NLSTEP` continuations `GENERAL`, `FIXED` or `ADAPT`,
and `MECH`. Legacy `NLPARM` is also accepted. The result is stored in
`model.analysis` together with `LGDISP` and `FOLLOWK`. The solver performs
fixed or adaptive stepping, bisecting cutbacks, and a tangent update at every
Newton iteration.

Unsupported nonlinear features include `ARCLN`, thermal or coupled controls,
material nonlinearity, contact, follower loads, partial stiffness updates,
and more than one stepping method on one control entry.

Nonzero prescribed displacements use the `D` field of `SPC` and are scaled by
the current load factor until their full value is reached. `SPCD` is not
implemented. `K6ROT` applies only to shell elements and does not stabilize
the inactive rotational degrees of freedom of a CROD.

The auditable MSC Nastran CROD campaign is stored in cases 92--127. Thirty-five
cases reach full load and are compared with explicit references; case 120
documents the expected force-control limit point. Cases 123--127 extend the
campaign to structures of 10--30 CRODs. See
[`NONLINEAR_CROD_CASES.md`](../tests/solver_tst/validation_cases/NONLINEAR_CROD_CASES.md)
and [`nonlinear_crod.md`](nonlinear_crod.md).

For CBAR, the validated property combination is `PBAR` plus
`PBARN1,SECT=S,BEH2=BAR`, with `INT2=LC` or `INT2=LS`. `I12`, `K1/K2`,
offsets, releases, and nonzero SPC values are supported. `LC` selects
linear/cubic Euler--Bernoulli interpolation and ignores `K1/K2`; `LS` selects
the one-point linear-shear element and uses those factors. Cases 128--166 and
their references are described in [`nonlinear_cbar.md`](nonlinear_cbar.md).
Mixed CROD/CBAR integration cases are 167--176.

## Units

MEF-solver does not convert units. Use one coherent system throughout the
model, for example:

- SI: m, N, Pa;
- a common Nastran engineering system: mm, N, MPa (`N/mm^2`).

Areas, inertias, moments, loads, and material constants must all match the
chosen system.

## Minimal CROD example

The complete deck
[`Ej1_v1.dat`](../tests/solver_tst/validation_cases/case_1/Ej1_v1.dat)
includes
[`model_ej1.bdf`](../tests/solver_tst/validation_cases/case_1/model_ej1.bdf).
It contains one CROD, one PROD, one MAT1, an axial force, and enough
constraints to remove every degree of freedom without stiffness.

Constraining inactive degrees of freedom is essential. The solver reserves
six degrees of freedom for every node even when an element contributes only
to some of them. An isolated CROD, for example, contributes axial and
torsional stiffness only along its local axis.

## Results

`LinearStaticSolver.solve()` returns:

```python
{
    "displacements": numpy.ndarray,
    "reactions": {node_id: numpy.ndarray},
    "stresses": {element_id: dict},
}
```

`NonlinearStaticSolver.solve()` retains those three keys and adds:

```python
{
    "increments": list[NonlinearIncrementResult],
    "iterations": list[NonlinearIterationResult],
    "load_factor": float,
    "strain_energy": float,
    "cutbacks": int,
}
```

- `displacements` stores six consecutive positions per node in `DofMapper`
  order;
- `reactions` maps constrained node IDs to six-component vectors;
- `stresses` maps each element ID to family-specific physical quantities.

The report is generated by
[`report_writer.py`](../src/postprocessing/report_writer.py) at
`output/<model-name>.rpt`. Its header identifies `LINEAR STATIC` or
`NONLINEAR STATIC`, then reports the model summary, solution sequence,
analysis type, and applicable global parameters.

For `SOL 400`, the report also includes:

- normalized `NLSTEP`/`NLPARM` parameters, step limits, Newton limits,
  convergence flags, and tolerances;
- final load factor, accepted increments, Newton iterations, cutbacks, and
  elastic strain energy;
- after the final state, the accepted-increment history and every Newton
  correction with load factor, residuals, and correction norm.

`FINAL CONVERGED STATE RESULTS` is organized by family:

- CROD: axial force, torque, stress, strain, and twist per unit length;
- CBAR: axial/torsional response, end A/B forces and moments, C--F recovery
  stresses, and, in `SOL 400`, current length, element energy, and extreme
  normal stresses;
- CTRIA3 and CQUAD4: separate tables for `Nx`, `Ny`, `Nxy`, `Mx`, `My`,
  `Mxy`, `Qx`, and `Qy`, and for top/bottom fiber stresses;
- displacements and reactions: independent nodal tables.

The report uses physical component names instead of generic numbered columns.
For shells, `TOP` is local `z=+t/2` and `BOTTOM` is local `z=-t/2`. For CBAR,
internal keys `P0`--`P3` are displayed as Nastran recovery points C--F. The
`NONLINEAR STATE IN CBAR ELEMENTS` table is emitted only when converged
nonlinear geometry and energy are present.

See [`result_report.md`](result_report.md) for column definitions, iteration
interpretation, and F06 correspondence.

## Running the tests

Run the complete collection from the repository root:

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

The currently verified result and the meaning of strict expected failures are
recorded in [`validation.md`](validation.md).
