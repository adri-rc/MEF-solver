# Technical handover

## Purpose and current status

MEF-solver is an academic finite element solver used for numerical comparison
with MSC Nastran. It implements a deliberately limited subset of Nastran and
must not be presented as a general replacement for it. The supported workflow
is:

```text
DAT/BDF -> pyNastran -> FEMModel -> linear or nonlinear solver -> RPT report
```

Development for the master's thesis is complete. The repository serves as an
academic reference and is not actively maintained. User support, issue triage,
and review of external contributions are not provided. Future-work items are
optional research suggestions, not commitments to further development.

The complete suite was run on 20 September 2026 after the publication asset
cleanup and reported `714 passed, 6 xfailed, 2 warnings in 422.63s`. The
nonlinear cleanup removed four obsolete tests together with the smaller legacy
bridge case. The final asset pass removed 30 duplicate, superseded, or invalid
validation artifacts, restored Case 18 to its canonical directory, and
canonicalized the Case 37 F06 name under Git LFS. Deck comments were translated
without changing any active Nastran line. Both warnings were
`PytestCacheWarning` messages caused by denied write access to `.pytest_cache`,
not numerical failures.

The active final nonlinear integration model is exclusively
[`case_nonlinear_test_final_updated`](../tests/solver_tst/validation_cases/case_nonlinear_test_final_updated/README.md).
The `_updated` suffix is retained because the delivered thesis refers to that
name; it does not imply that another public validation case is required.

## Supported scope

### Linear static analysis

The `SOL 101` path supports small deformations, isotropic `MAT1`, six degrees
of freedom per node, and the following elements:

- `PROD/CROD`;
- `PBAR/CBAR`;
- `PSHELL/CTRIA3`;
- `PSHELL/CQUAD4` for coplanar geometry.

The supported input subset also includes `FORCE`, `MOMENT`, nested `LOAD`
combinations, `SPC`, `SPC1`, `SPCADD`, `GRID.PS`, prescribed nonzero values,
explicit subcase selection, local `CP/CD/CID` coordinate systems, and
`PARAM,K6ROT` for shell drilling stiffness.

### Geometrically nonlinear static analysis

The `SOL 400` path requires `ANALYSIS=NLSTAT` and supports geometric
nonlinearity for isotropic elastic `MAT1` models containing `PROD/CROD` and/or
`PBAR/PBARN1/CBAR`. It uses conservative nodal loads, full Newton iterations,
and `FIXED` or `ADAPT` increments. The active element implementations are
[`crod_nonlinear.py`](../src/solver/elements/crod_nonlinear.py) and
[`cbar_nonlinear.py`](../src/solver/elements/cbar_nonlinear.py).

This route does not implement material plasticity, contact, buckling
assessment, arc-length control, follower loads, unilateral cable/slackening
behaviour, distributed loads, solids, or shells. Exact formulation and
validation boundaries are documented in [nonlinear CROD](nonlinear_crod.md)
and [nonlinear CBAR](nonlinear_cbar.md).

## Architecture and numerical invariants

The main execution boundary is intentionally small:

| Layer | Canonical implementation | Responsibility |
|---|---|---|
| CLI orchestration | [`main.py`](../main.py) | Read a model, select the solver, write a report |
| Input adapter | [`nastran_reader.py`](../src/preprocessing/nastran_reader.py) | Validate and map pyNastran entities into `FEMModel` |
| Data model | [`data_model.py`](../src/preprocessing/data_model.py) | Format-independent nodes, materials, properties, elements, loads, constraints, and analysis controls |
| Linear solution | [`linear_static_solver.py`](../src/solver/linear_static_solver.py) | Sparse assembly, prescribed/free partition, solution, reactions, and recovery |
| Nonlinear solution | [`nonlinear_analysis.py`](../src/solver/nonlinear_analysis.py) | Increment control, residuals, tangents, cutbacks, Newton convergence, and recovery |
| Presentation | [`report_writer.py`](../src/postprocessing/report_writer.py) | F06-like text report only; it does not own numerical recovery |

Preserve these invariants when modifying the project:

- complete decks must select supported `SOL 101` or `SOL 400`; bulk-only
  punch files retain the documented linear default;
- unsupported elements, loads, constraints, rigid elements, MPCs, and DOF-set
  cards fail explicitly instead of being silently omitted;
- `GRID.CP` positions and load vectors in their `CID` are mapped into basic
  coordinates during preprocessing; assembly then accounts for `GRID.CD`;
- there is no general `AUTOSPC`: a real mechanism or an unconstrained active
  degree of freedom must fail;
- a blank PSHELL `MID2` disables bending and a blank `MID3` disables physical
  transverse-shear flexibility; neither inherits another material field;
- supported CQUAD4 requires `MID1` and coplanar geometry;
- CQUAD4 stiffness must not be tuned solely to improve centroidal `Qx/Qy`
  recovery without a new auditable DMIG reference;
- all model units must be coherent because the solver performs no conversion.

See [architecture.md](architecture.md) for the full coordinate-transformation,
assembly, solution, recovery, and reporting flow, and [formulation.md](formulation.md)
for the governing equations.

## Test map

The suite is organized by responsibility:

- `tests/preprocessing_tst/test_nastran_reader.py` covers card mapping and
  fail-fast rejection of unsupported structural entities;
- `tests/preprocessing_tst/test_nastran_reader_end_to_end.py` and
  `test_case_control_and_coordinates.py` cover complete decks, subcases,
  nested load/constraint sets, and coordinate transformations;
- `tests/preprocessing_tst/test_nonlinear_analysis_control.py` covers the
  `SOL 101`/`SOL 400`/punch distinction and normalized nonlinear controls;
- `tests/solver_tst/elements/` covers element matrices, invariants, invalid
  geometry, material-field semantics, and nonlinear element states;
- `tests/solver_tst/test_linear_static_solver.py` and
  `test_solver_regressions.py` cover global linear assembly and failure modes;
- `tests/solver_tst/test_nonlinear_static_solver.py` covers Newton assembly,
  increments, cutbacks, convergence, and state history;
- `tests/solver_tst/test_nonlinear_crod_*`,
  `test_nonlinear_cbar_*`, and `test_nonlinear_mixed_validation_end_to_end.py`
  compare the nonlinear campaigns with transcribed MSC references;
- `tests/solver_tst/test_nonlinear_final_bridge.py` checks the authoritative
  202-node bridge topology, constraints, loads, and converged solution;
- `tests/postprocessing_tst/test_report_writer.py` owns the report-format
  contract independently of the numerical solvers.

Reference values are literals extracted from preserved MSC output. Tests do
not parse F06/PCH files dynamically, so a missing or changed reference artifact
cannot silently alter an assertion. The complete validation rationale and
error measures are in [validation.md](validation.md).

## Auditable validation cases

### Linear

- [`case_linear_test_final`](../tests/solver_tst/validation_cases/case_linear_test_final/README.md):
  planar frame/web coupon with 384 shell elements, 76 CBAR, and 4 CROD; full
  response errors against MSC remain below `1.5e-5` in relative L2 norm.
- [`case_linear_test_final_T`](../tests/solver_tst/validation_cases/case_linear_test_final_T/README.md):
  three-dimensional T frame with 417 GRID and 506 elements; global response,
  shell membrane resultants, and fibre stresses correlate with MSC.
- [`case_linear_test_final_T_no_CTRIA_MID3`](../tests/solver_tst/validation_cases/case_linear_test_final_T_no_CTRIA_MID3/README.md):
  near-singular diagnostic that isolates CTRIA3 transverse-shear recovery; it
  is not a physically valid alternative model.

The CTRIA3 `MID1/MID2/MID3`, `TS/T`, and `12I/T^3` campaigns and the CQUAD4
DMIG/E2E matrix cases remain the element-level reference basis.

### Nonlinear

- Cases 92--127 cover CROD small-displacement limits, finite extension,
  stepping, coordinates, topology, tangent probes, torsion, and prescribed
  values. Case 120 documents expected force-control failure near a limit point.
- Cases 128--166 cover nonlinear CBAR axial, bending, torsion, orientation,
  offsets, releases, section recovery, tangent checks, and `INT2=LC/LS`.
- Cases 167--176 contain mixed CROD/CBAR structures with increasing size and
  combined coordinate, load, property, and prescribed-displacement features.
- The only final bridge case is
  [`case_nonlinear_test_final_updated`](../tests/solver_tst/validation_cases/case_nonlinear_test_final_updated/README.md):
  202 GRID, 144 CBAR, 428 CROD, 44 forces, 3 moments, explicit foundations,
  imposed settlements, and 20 fixed increments. Relative L2 errors are about
  `10^-5` for nodal/global and axial quantities; CBAR fibres retain a separate
  `1.673e-2` L2 envelope and a localized 11.57% significant maximum.

Campaign indexes are preserved in
[`NONLINEAR_CROD_CASES.md`](../tests/solver_tst/validation_cases/NONLINEAR_CROD_CASES.md),
[`NONLINEAR_CBAR_CASES.md`](../tests/solver_tst/validation_cases/NONLINEAR_CBAR_CASES.md),
and
[`NONLINEAR_FINAL_CLOSURE_CASES.md`](../tests/solver_tst/validation_cases/NONLINEAR_FINAL_CLOSURE_CASES.md).

## Expected failures and explicit limits

Six strict `xfail` cases define CQUAD4 scope rather than known regressions:

| Cases | Reason | Required behaviour |
|---|---|---|
| 54, 59, 62, 63, 64 | Warped quadrilateral | Reject geometry outside the coplanar formulation |
| 78 | Default `SNORM` smoothing | Remain unsupported; the `SNORM,0` variant is validated separately |

Centroidal CTRIA3 moments/shears in the final T frame and CQUAD4 centroidal
`Qx/Qy` are documented recovery limitations. Do not hide either limitation by
loosening unrelated global tolerances or recalibrating already validated
stiffness.

## Environment and reproducible commands

The current local environment uses Python 3.12.14, NumPy 1.26.4, SciPy 1.13.1,
pandas 2.2.3, pyNastran 1.4.1, and pytest 8.3.5. `requirements.txt` currently
contains compatible lower bounds rather than an exact lock.

From PowerShell:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
git lfs pull
.\.venv\Scripts\python.exe -m pytest -q
```

Useful focused runs are:

```powershell
.\.venv\Scripts\python.exe -m pytest -q tests\preprocessing_tst
.\.venv\Scripts\python.exe -m pytest -q tests\solver_tst\elements
.\.venv\Scripts\python.exe -m pytest -q tests\solver_tst\test_nonlinear_static_solver.py
.\.venv\Scripts\python.exe -m pytest -q tests\solver_tst\test_nonlinear_final_bridge.py
```

Run one model through the CLI with:

```powershell
.\.venv\Scripts\python.exe main.py --file input\model.dat
```

Select `--subcase <ID>` when a complete deck has more than one eligible
subcase. Reports are written under ignored `output/` paths. DAT/BDF, F06, and
image evidence are managed through Git LFS; PCH references remain plain text.
Confirm LFS-backed files are materialized before diagnosing a parser or
numerical failure.

## Publication preparation and optional future work

The project's original source code is licensed under
[BSD-3-Clause](../LICENSE). The [README](../README.md) records the thesis
purpose, absence of warranties, limitation of liability, and scope of the
MSC Nastran comparisons. The included
[master's thesis (PDF)](thesis/MEF-solver_MastersThesis_AdrianRuizCalderon.pdf)
is available as an academic reference. Code, tests, and the retained
technical documentation define the current implementation.

Before making the repository public, run a final secret, absolute-path, and
Git LFS audit, verify the documentation links, and record the final test-suite
result. These publication checks do not require extending the solver's scope.

Possible follow-up work for independent study includes continuous integration
with Git LFS checkout, `pyproject.toml` and/or an exact dependency lock,
mesh-refinement convergence studies, and an explicit CLI batch-error policy.
These ideas and the extensions in [roadmap.md](roadmap.md) have no committed
schedule and do not imply ongoing maintenance or support.

Do not commit generated `output/`, caches, `draft_*` material, or local thesis
sources. Thesis material under `docs/TFM/` or `docs/draft_TFM/` is outside the
repository source of truth and must remain untouched unless explicitly
requested.
