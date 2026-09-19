# Roadmap

## Verified baseline

The DAT/BDF to data model to sparse assembly to solution to report workflow is
operational for CROD, CBAR, CTRIA3, and coplanar CQUAD4. The complete suite
reported `718 passed, 6 xfailed, 2 warnings` on 19 September 2026. The six
strict expected failures are documented CQUAD4 scope limits.

The linear implementation currently includes:

- independent PSHELL `MID1`, `MID2`, and `MID3` semantics for CTRIA3 and
  CQUAD4;
- CTRIA3 application of `12I/T^3` to bending and `TS/T` to transverse
  shear;
- validated coplanar CQUAD4 stiffness through cases 90 and 91;
- report generation isolated in `src/postprocessing/report_writer.py`;
- two large linear integration models and one targeted CTRIA3 shear
  diagnostic.

The first integration model, `case_linear_test_final`, contains 384 shell
elements, 76 CBAR, and 4 CROD. Its unmodified deck completed in MSC Nastran
2026.1, and full-field relative L2 errors remain below `1.5e-5`.

`case_linear_test_final_T` extends the model to a three-dimensional T frame
with 417 GRID and 506 elements. Global response, shell membrane resultants,
and fibre stresses correlate with MSC Nastran. Centroidal CTRIA3 moments and
shears do not: their L2 errors are `0.682` and `17.30`. This is an explicit
recovery limitation and must not be hidden by changing stiffness or relaxing a
tolerance.

`case_linear_test_final_T_no_CTRIA_MID3` changes only the CTRIA3 transverse
shear material to `E=1e-20`. It reduces the full-field CTRIA3 moment error to
`7.106e-5`, but produces a nearly singular physical model with
`cond(K)=5.660208e28`. The case is a diagnostic, not an alternative validated
structure.

The nonlinear implementation distinguishes `SOL 101` from `SOL 400`,
normalizes `NLSTEP`/`NLPARM`, `LGDISP`, `FOLLOWK`, and relevant
`PBARN1` fields, and solves geometric nonlinearity for `MAT1` CROD/CBAR
models. The audited MSC campaigns cover CROD cases 92--127, CBAR cases
128--166, mixed cases 167--176, and two final bridge models.

## Completed: geometric nonlinear CROD

1. Preserved auditable MSC results for cases 92--127.
2. Transcribed intermediate states, tangent probes, and final results.
3. Implemented `CRODNonlinear` and `NonlinearStaticSolver` without changing
   the linear formulation.
4. Integrated CLI routing, recovery, and nonlinear report history.
5. Documented the force-controlled limit point in case 120 and the absence of
   arc-length control.

## Completed: geometric nonlinear CBAR

1. Implemented a corotational `CBARNonlinear` without modifying the linear
   path.
2. Integrated `I12`, `K1/K2`, offsets, releases, `G0` orientation, loads,
   and prescribed displacements.
3. Identified logarithmic axial strain with cases 130/131 and checked the
   tangent with cases 137--139.
4. Added 39 explicit CBAR end-to-end references; all pass.
5. Validated `INT2=LC` as Euler--Bernoulli and `INT2=LS` as one-point
   linear shear interpolation with cases 154--160.
6. Added converged geometry, end actions, and C--F stresses to the report
   without changing the linear tables.

Some near-neutral fibres and severe combined-rotation segments in cases
165/166 remain diagnostics. They do not invalidate the closed global response
and centroidal axial-stress scope.

## Priority 1: public-release consolidation

The preprocessing and linear-analysis cleanup is the first publication pass.
It covers active-module naming, removal of disconnected prototypes, English
comments/docstrings, and English documentation checked against code and tests.
The nonlinear implementation will receive the same focused cleanup in a later
pass.

Remaining consolidation work:

1. add continuous integration for the full pytest suite;
2. add `pyproject.toml` and/or a reproducible dependency lock;
3. add mesh-refinement convergence tests in addition to DMIG and E2E checks;
4. decide whether batch mode stops on the first error or accumulates failures;
5. complete a public-release audit of nonlinear code and its dedicated
   documentation;
6. choose and add a repository license before public release.

The CTRIA3 centroidal `M/Q` discrepancy in the T frame is already isolated as
a formulation/recovery limit. It is not a consolidation task and does not
justify recalibrating validated stiffness.

## Priority 2: performance

1. assemble from vectorized COO triplets;
2. cache transformations and element solver objects between assembly and
   recovery;
3. measure parsing, element stiffness, assembly, solution, and recovery on
   increasing mesh sizes;
4. investigate reordering and iterative solvers when model size warrants it.

## Priority 3: new capabilities

Every physical extension must include a minimal model, an auditable external
reference, and an end-to-end test before it is advertised as supported.
Candidates include:

- shell pressure, gravity, and other distributed loads;
- MPC and rigid elements;
- advanced beam properties;
- orthotropic and laminated materials;
- higher-order and solid elements;
- export to visualization formats.

## Optional CQUAD4 work

- Investigate MSC Nastran's centroidal `Qx/Qy` reporting convention without
  changing validated stiffness.
- Add recovery at additional stations or corners.
- Implement `SNORM` only if a supported faceted model requires it.
- Keep rejecting warped quadrilaterals. A future warped formulation would be a
  new element, not a relaxed tolerance in the current coplanar one.
