# MEF-solver

> **Academic disclaimer**
>
> This project was developed solely to complete a master's thesis (TFM).
> The software, models, documentation, and results may contain errors,
> omissions, and limitations. They are provided "as is", without warranties
> of correctness, accuracy, completeness, or fitness for a particular purpose.
> To the maximum extent permitted by applicable law, the author and
> contributors disclaim liability for errors, use of this project, or any
> resulting loss or damage. Any engineering use requires independent
> verification by a qualified professional. See [LICENSE](LICENSE).
>
> MSC Nastran results are used as external references for assessing the
> implementation developed for the thesis. This project does not assess,
> certify, endorse, or warrant the correctness of MSC Nastran or its results,
> and does not claim to reproduce or replace that software. The project is
> independent and is not endorsed or certified by MSC Nastran's developers
> or vendors.

MEF-solver is an academic finite element solver for static structural analysis.
Its supported analyses, elements, and input formats are deliberately limited
to the documented scope.

For questions, implementation details, and mathematical formulations, consult
the [master's thesis (PDF)](docs/thesis/MEF-solver_MastersThesis_AdrianRuizCalderon.pdf).

## Project status

This repository contains the completed work developed for a master's thesis
and is provided as an academic reference. It is not actively maintained.
No user support, handling of bug reports, or review of external contributions
is provided. Readers may fork and adapt the original source code under the
[BSD 3-Clause License](LICENSE).

The [roadmap](docs/roadmap.md) records possible directions for future academic
work; it does not constitute a development plan or a commitment to implement
changes.

## Capabilities

The linear static path supports small deformations, isotropic `MAT1`
materials, and six degrees of freedom per node for:

- `CROD` with `PROD`;
- `CBAR` with `PBAR`;
- `CTRIA3` with a subset of `PSHELL`;
- coplanar `CQUAD4` with a subset of `PSHELL`.

The geometric nonlinear `SOL 400` path supports `MAT1` models made of
`CROD` and/or `CBAR`, conservative nodal loads, full Newton-Raphson, and
`FIXED` or `ADAPT` load increments. See
[nonlinear CROD](docs/nonlinear_crod.md) and
[nonlinear CBAR](docs/nonlinear_cbar.md) for its exact scope.

The reader accepts linear `SOL 101` and nonlinear-static `SOL 400` with
`ANALYSIS=NLSTAT`. Other solution sequences and inconsistent `SOL`/`ANALYSIS`
pairs are rejected explicitly.

Supported input also includes:

- `FORCE`, `MOMENT`, and nested `LOAD` combinations;
- `SPC`, `SPC1`, `SPCADD`, `GRID.PS`, and nonzero prescribed values;
- explicit `SUBCASE` selection;
- `GRID.CP` positions, `GRID.CD` degrees of freedom, and loads defined in
  a local `CID`;
- `PARAM,K6ROT` for CTRIA3 and CQUAD4.

## Quick start

Python 3.12 is the currently verified interpreter. From PowerShell:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
git lfs pull
.\.venv\Scripts\python.exe main.py --file input\model.dat
```

Select a subcase explicitly when a deck contains more than one:

```powershell
.\.venv\Scripts\python.exe main.py --file input\model.dat --subcase 2
```

Generated reports are written to `output/<model>.rpt`. The output directory
is intentionally ignored by Git.

## Verification

Run the complete suite with:

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

The verified baseline on 20 September 2026 is `714 passed, 6 xfailed`. The
six strict expected failures describe explicit CQUAD4 scope limits: five
warped geometries and one faceted model that requires unsupported `SNORM`
smoothing. They are not failures of the supported coplanar element.

Validation assets include element matrices, end-to-end decks, and preserved
MSC Nastran 2026.1 references. The main linear integration models are:

- [airframe section](tests/solver_tst/validation_cases/case_linear_test_final/README.md),
  with 384 shell elements, 76 CBAR, and 4 CROD;
- [three-dimensional T frame](tests/solver_tst/validation_cases/case_linear_test_final_T/README.md),
  with 417 GRID and 506 elements;
- [CTRIA3 shear diagnostic](tests/solver_tst/validation_cases/case_linear_test_final_T_no_CTRIA_MID3/README.md),
  which isolates a known recovery limitation and is not a physically valid
  alternative model.

The T-frame correlates global response, shell membrane resultants, and fibre
stresses. Its centroidal CTRIA3 `M/Q` output remains a documented recovery
limitation and is not hidden behind a relaxed tolerance. CQUAD4 stiffness must
not be recalibrated solely to improve centroidal `Qx/Qy` without a new DMIG
reference.

## Important modelling semantics

- No general `AUTOSPC` is applied. A structural mechanism must fail rather
  than be stabilized silently.
- PSHELL material fields are independent: blank `MID2` disables bending and
  blank `MID3` disables physical transverse-shear flexibility. They do not
  inherit another MID.
- The supported CQUAD4 subset requires `MID1`; bending-only CQUAD4 is rejected.
- CTRIA3 applies `12I/T^3` to bending and `TS/T` to transverse shear.
- CQUAD4 rejects warped geometry and supports independent facet normals, not
  `SNORM` smoothing.
- The program performs no unit conversion; every deck must use one coherent
  unit system.

## License

The project's original source code is licensed under the
[BSD 3-Clause License](LICENSE). Existing third-party notices remain
applicable.

## Documentation

- [Master's thesis (PDF)](docs/thesis/MEF-solver_MastersThesis_AdrianRuizCalderon.pdf)
- [Usage](docs/usage.md)
- [Architecture](docs/architecture.md)
- [Mathematical formulation](docs/formulation.md)
- [Validation](docs/validation.md)
- [Result report reference](docs/result_report.md)
- [CQUAD4 implementation](docs/cquad4.md)
- [Project roadmap](docs/roadmap.md)
- [Technical handover](docs/HANDOVER.md)

Repository contributors and AI assistants must also follow
[AGENTS.md](AGENTS.md).
