# Validation

## Verified state

The complete suite was run on 20 September 2026 with:

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

Result:

```text
714 passed, 6 xfailed, 2 warnings
```

Both warnings are `PytestCacheWarning` messages caused by denied write access
to `.pytest_cache`; they do not affect test execution.

## Validation layers

1. Data model, parser, case-control selection, and coordinate systems.
2. Shape functions, Jacobians, transformations, and quadrature.
3. Element invariants: symmetry, rigid modes, orientation, decoupling, and
   invalid geometry.
4. Element/global stiffness matrices against Nastran DMIG references.
5. End-to-end displacements, applied loads, reactions, and element results
   against values transcribed from preserved Nastran output.

Report formatting has a separate contract in
`tests/postprocessing_tst/test_report_writer.py`. Case 80 is the principal
mixed linear regression because one model combines CBAR, CTRIA3, and CQUAD4
and checks assembly and recovery across families.

Preprocessing also has explicit fail-fast coverage. Complete decks must select
`SOL 101` or `SOL 400`; punch files retain the documented linear default; and
`SOL 400` must select `ANALYSIS=NLSTAT`. Unsupported constraint and DOF-set
cards are rejected instead of being silently ignored, including `MPC`, rigid
elements, `SUPORT`, `ASET`/`BSET`/`CSET`/`QSET`/`USET`, `OMIT`, `SPCOFF`, and
superelement `RELEASE` variants.

## CTRIA3 PSHELL material semantics

[case_ctria3_pshell_material_semantics](../tests/solver_tst/validation_cases/case_ctria3_pshell_material_semantics/README.md)
contains seven one-element decks executed with MSC Nastran Student Edition
2026.1. Distinct `MID1`, `MID2`, and `MID3` moduli and `K6ROT=0` make
the `KAAX` DMIG blocks independently observable.

The references establish that:

- `MID1` alone produces membrane without plate stiffness;
- `MID2` alone produces plate stiffness without membrane;
- `MID2+MID3` without `MID1` is valid;
- `MID3` alone and `MID1+MID3` produce
  `USER FATAL MESSAGE 9994 (BULAEM)`, because `MID3` requires `MID2`.

The five valid models preserve F06/PCH output. The two rejected models preserve
their F06. Automated checks cover those artifacts, inactive zero blocks,
membrane DMIG values, blank `MID1` parsing, and the absence of `MID2/MID3`
fallbacks. They do not claim complete equality of the Reissner--Mindlin plate
block.

## CTRIA3 PSHELL factors

[case_ctria3_pshell_factors_k](../tests/solver_tst/validation_cases/case_ctria3_pshell_factors_k/README.md)
contains four DMIG models with `AUTOSPC=NO` and `K6ROT=0`. Two vary
`TS/T` between `0.25` and `1.25`; two vary `12I/T^3` between `0.40`
and `1.60`.

The bending matrices match term by term with relative Frobenius errors below
`1.9e-10`. Shear-control errors are `9.96e-5` and `4.98e-4`; the largest
relative errors among selected terms are `2.46e-4` and `1.23e-3`. These
bounded differences are consistent with the documented CTRIA3 plate
limitation and are not hidden by one broad absolute tolerance.

[case_ctria3_pshell_factors_e2e](../tests/solver_tst/validation_cases/case_ctria3_pshell_factors_e2e/README.md)
contains three four-CTRIA3 cantilevers. All finish in MSC Nastran 2026.1 with
`AUTOSPC=NO` and no `FATAL`. Selected-value maxima are:

| Quantity | Maximum difference |
|---|---:|
| Displacements | 1.10% |
| Reactions | 0.61% |
| Membrane results | 2.7e-6% |
| Moments | 4.02% |
| Fibre stresses | 4.41% |
| Significant shears | 0.87% |

Near-zero shear components differ by less than `8.3e-3 N/mm` and use a
separate absolute threshold. These cases confirm that `TS/T` replaces the
old fixed `5/6` and that `12I/T^3` separately scales bending and the fibre
section modulus. They do not establish general equivalence of CTRIA3
centroidal `M/Q` recovery.

## Final linear integration model

`case_linear_test_final` is a reproducible `1200 x 600 mm` planar
frame/web coupon with two circular cut-outs. It contains:

- 256 CQUAD4;
- 128 CTRIA3 around the cut-outs;
- 76 CBAR cap/stiffener members;
- 4 CROD load-introduction members;
- two forces, one moment, six fixed stations, and one prescribed
  `0.8 mm` displacement.

`test_linear_final_validation_case.py` checks the input inventory, finite
solution, prescribed value, global force/moment equilibrium, and explicit MSC
references for every element family. The unmodified deck completed normally
in MSC Nastran Student Edition 2026.1 on 21 August 2026.

Full OP2-array relative L2 errors were:

| Field | Relative L2 error |
|---|---:|
| Displacements | `1.052e-6` |
| Reactions | `1.450e-5` |
| CROD/CBAR results | below `5.4e-6` |
| CTRIA3/CQUAD4 results | below `1.1e-5` |

The preserved F06 is auditable, while pytest compares transcribed values and
does not read F06/OP2 dynamically. This is numerical cross-validation of one
deck, not experimental correlation or general validation of Nastran. See the
[case README](../tests/solver_tst/validation_cases/case_linear_test_final/README.md).

## Three-dimensional T-frame model

`case_linear_test_final_T` adds a perpendicular
`1200 x 300 mm` CQUAD4 flange to the perforated web. Six diagonals split
into twelve CBAR members form three symmetric ribs. The complete inventory is
417 GRID, 288 CQUAD4, 128 CTRIA3, 88 CBAR, and 2 CROD. It includes four
discrete flange connections, two `6000 N` axial forces, a
`250000 Nmm` moment, and two prescribed `T2=0.2 mm` values with
`AUTOSPC=NO`.

`test_linear_final_t_validation_case.py` checks geometry, inventory,
boundary conditions, equilibrium, and literal F06 references for 15 GRID,
six reactions, both CROD, and 15 elements from every other family.

Full OP2-array relative L2 errors were:

| Field | Relative L2 error |
|---|---:|
| Translations | `8.206e-3` |
| Rotations | `2.547e-2` |
| Force reactions | `1.705e-4` |
| Moment reactions | `3.553e-4` |
| CBAR normal stresses | `3.930e-3` |
| CQUAD4/CTRIA3 membrane | `3.268e-5 / 7.218e-6` |
| CQUAD4/CTRIA3 fibre stresses | `6.582e-3 / 1.020e-2` |

Centroidal CTRIA3 moments and shears are explicitly not accepted as
correlated: their L2 errors are `0.682` and `17.30`. This known recovery
limit is neither covered by a wide tolerance nor used to retune stiffness.
The MSC Nastran 2026.1 F06 reaches `END OF JOB` with `AUTOSPC=NO`. See the
[case README](../tests/solver_tst/validation_cases/case_linear_test_final_T/README.md).

### Complementary diagnostic with negligible MID3

`case_linear_test_final_T_no_CTRIA_MID3` changes only `PSHELL 11.MID3` to
a new `MAT1 99` with `E=1e-20`. MSC reaches `END OF JOB`, but the
stiffness condition number is `5.660208e28`. This is a diagnostic probe, not
a physical alternative model.

The full-field CTRIA3 moment error falls from `0.682` to `7.106e-5`, and
both implementations recover engineering-zero shears: norms below
`4.1e-11` locally and `8.7e-24` in MSC. Other element families, reactions,
and DOFs unrelated to `T3` retain their correlation. However, GRID 203
transverse translation becomes `7.3586e10 mm` locally versus `-0.3421 mm`
in MSC, confirming a near-mechanism.

Together, the physical and diagnostic models isolate the `M/Q` discrepancy
to transverse-shear treatment/recovery without validating the degenerate
structure. See the
[diagnostic README](../tests/solver_tst/validation_cases/case_linear_test_final_T_no_CTRIA_MID3/README.md).

## CQUAD4 stiffness matrices

`test_E2E_CQUAD_Nastran_K_solution.py` applies complementary criteria:

- relative Frobenius error for the whole matrix;
- term-by-term `allclose`;
- element-axis blocks for membrane, plate/shear, drilling, and coupling.

The cases cover membrane, bending, their combination, `K6ROT=0/100`, blank
and active `MID3`, `G3` and `TS/T` changes, rectangles, trapezoids,
skew affine quadrilaterals, and arbitrarily oriented planes.

Case 90 exposed an overly soft former skew law. Keeping the edge operator,
residual flexibility, and `MID3` contribution while changing only the
rational torsional metric gives:

- full-matrix Frobenius error below `1e-6`;
- plate/shear block error below `4e-6`;
- term checks with `rtol=1e-5`, `atol=2e-3`.

Loaded case 91 confirms the adjustment on two highly skew panels:

- displacement-vector error: `2.07e-5`;
- reaction-vector error: `6.30e-5`;
- `N`, `M`, and fibre stresses: below 0.1% per component;
- largest `Qx/Qy` difference: 47.1%, or 0.0216 absolute.

The `Q` difference belongs to centroidal recovery, not stiffness; there is no
associated discrepancy in displacement, reaction, or moment.

Cases 57 and 58 are severely ill-conditioned
(`2.90e25` and `8.52e26`) and retain documented sensitivity. Case 59 is
warped and is expected to be rejected.

## Meaning of the six expected failures

| Cases | Reason | Interpretation |
|---|---|---|
| 54, 59, 62, 63, 64 | Warped CQUAD4 | Outside the coplanar domain; rejection is expected |
| 78 | Default `SNORM` smoothing | Not implemented; use the `SNORM,0` variant |

Cases 82 and 88 validate faceted assemblies with independent normals. Former
CBAR expected failure 142 was closed by the `LC/LS` work in cases 154--160.

The coplanar CQUAD4 implementation is considered closed because rigid modes
and invariants are covered, active DMIG matrices match globally and by block,
multi-element models reproduce global response, and unsupported geometry is
rejected explicitly. More accurate replication of centroidal `Qx/Qy`
printing remains optional recovery work.

## Final nonlinear integration model

`case_nonlinear_test_final_updated` is the final three-span bridge documented
in the submitted thesis. It contains 202 GRID and 572 elements: 144 CBAR and
428 CROD. Both solvers complete 20 increments without cutbacks. Relative L2
errors are `3.280e-5` in translations, `2.416e-5` in rotations,
`2.898e-5` in force reactions, `2.125e-5` in CROD stress, and
`1.379e-5` in centroidal CBAR axial stress. Full CBAR fibre results have
`1.673e-2` L2 error and a localized 11.57% significant maximum. Separate
envelopes prevent this recovery difference from relaxing global or axial
checks.

## Nonlinear CROD campaign

Cases 92--127 cover the small-displacement limit, finite tension/compression,
fixed and adaptive stepping, `CD/CID`, series/parallel layouts, `LOAD`,
2D trusses, 3D structures, tangent probes, torsion, prescribed values, and
models of 10--30 rods.

Thirty-five cases reach full load in MSC and locally. Case 120 preserves the
expected non-convergence near `lambda=0.98387` at a force-controlled limit
point. Cases 102, 105, and 125 cross snap-through by jumping to a later branch;
MSC accumulated work in those paths is not asserted as instantaneous elastic
energy. Explicit final values are stored in tests rather than read from F06 at
runtime.

## Nonlinear CBAR and mixed campaign

Cases 128--166 comprise 39 CBAR decks covering the linear limit, axial
response, biaxial bending, moments, torsion, combined loading, a central
tangent probe, `G0`, load `CID`, `I12`, `K1/K2`, prescribed values,
multiple properties/materials, refinement, releases, offsets, and C--F
recovery.

Cases 154--160 establish `INT2=LC` as Euler--Bernoulli and `INT2=LS` as
one-point linear shear interpolation. Their maximum vector errors are
`8.04e-5` in displacements, `6.03e-5` in reactions, and `1.52e-4` in
energy.

Cases 167--176 mix 2--30 CROD/CBAR elements. All 16 closure decks (161--176)
pass preprocessing, have full-rank initial stiffness where applicable, reach
full MSC/local load, and contain no MSC `FATAL`. The independent mixed E2E
module covers those ten decks and the final bridge. It uses a 3% relative
envelope for the numbered campaign; observed preselection vector errors remain
below 0.04% in displacement, 0.03% in reaction, 2.2% in CBAR fibres, and 1.3%
in CROD results. Near-neutral fibres and the most rotated portions of cases
165/166 remain documented diagnostics rather than `xfail` cases.
