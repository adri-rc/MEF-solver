# Geometrically nonlinear CROD analysis

## Status and scope

This document describes the input contract, the corotational element in
`src/solver/elements/crod_nonlinear.py`, and the incremental solver for
geometrically nonlinear static analysis with CROD. `LinearStaticSolver`
remains physically unchanged; `main.py` routes SOL 400 to
`src/solver/nonlinear_analysis.py` and SOL 101 to the linear path.

The CROD scope implemented is:

- isotropic elastic `MAT1`, `PROD`, and `CROD`;
- large displacements and rotations with finite axial strain;
- conservative nodal `FORCE`/`MOMENT` loads and `LOAD` combinations;
- normalized `SPC`, `SPC1`, `SPCADD`, and `GRID.PS` constraints, including
  nonzero `SPC` values applied incrementally;
- `FIXED`/`ADAPT` stepping, full Newton-Raphson, residual evaluation,
  cutbacks, and a consistent tangent;
- no plasticity, contact, arc-length control, or follower loads.

CBAR is also supported by the nonlinear solver through its own corotational
formulation and validation campaign. It is documented separately in
[`nonlinear_cbar.md`](nonlinear_cbar.md); CROD must not be treated as if it
included bending or beam-section kinematics.

The SOL 400 guide states that its 1D elements support large
displacements and rotations. Large material strain is a separate option and
is not enabled here. Results from cases 94--96 with `LGDISP=1` give
`epsilon=(L-L0)/L0`, `sigma=E*epsilon`, and `N=A*sigma` using the reference
area; they do not show an isochoric update of `A`. This is the law reproduced
by the solver, which is limited to geometric nonlinearity.

## Main MSC Nastran entries

### Analysis selection

```text
SOL 400
CEND
SUBCASE 1
  ANALYSIS = NLSTAT
  SPC = 10
  LOAD = 20
  NLSTEP = 30
```

- `SOL 400` selects MSC Nastran's implicit nonlinear solver.
- `ANALYSIS=NLSTAT` identifies a nonlinear static step. The reader also
  normalizes the `NLSTATIC` and `NLSTATICS` aliases.
- `SPC` and `LOAD` select the boundary conditions and the **total** load value
  for the step. `NLSTEP` defines how that value is reached.
- The `NLSTEP=30` Case Control integer must match the `SID` of an `NLSTEP`
  Bulk Data entry.

### Geometric nonlinearity

```text
PARAM,LGDISP,1
PARAM,FOLLOWK,NO
```

`PARAM,LGDISP,1` activates large displacements. This requires updating the
geometry, internal forces, and tangent stiffness; it is not equivalent to
repeatedly solving the linear elastic matrix without a geometric term.

Preprocessing preserves the complete value: `-1` disables geometric
nonlinearity, `0` retains differential stiffness only, `1` activates large
displacements and follower loads, and `2` activates large displacements
without follower loads. The campaign uses `-1` in control cases 92/121 and
`1` in the remaining cases except control case 118, which uses `2`.

MSC also enables follower-load treatment when large displacements are used.
The initial campaign sets `PARAM,FOLLOWK,NO` to exclude follower-load
stiffness and isolate the structural tangent. `FORCE` entries are defined in
the basic system or in a fixed `CID`. The project implementation remains
limited to conservative nodal loads; any extension to follower loads requires
independent references.

### Prescribed displacements

The preprocessor and `NonlinearStaticSolver` support prescribed displacements
through a nonzero value in an `SPC` entry. For example:

```text
SPC1,124,23456,5
SPC,124,5,1,-4.0
```

GRID 5 is constrained in all six DOFs and reaches `T1=-4.0` when the load
factor reaches one. During the solution, `T1(lambda)=-4.0*lambda` is
prescribed consistently with the other actions in the step. Cases 124 and
126 exercise nodal forces together with this type of motion. The `SPCD` load
entry is not part of the supported reader subset; it must not be confused
with a nonzero `SPC` value.

### Why `K6ROT` does not appear in the CROD cases

`PARAM,K6ROT` is numerical stiffness for the drilling rotation, namely the
rotation about the normal of shell elements such as CTRIA3 and CQUAD4. It
does not add physical stiffness or stabilization to a CROD. In a CROD truss
with `J=0`, nodal rotational DOFs do not contribute to the energy and are
constrained explicitly with `SPC`; leaving them free would create a real
mechanism in the discrete model. If `PROD.J>0`, only the relative rotation
projected onto the rod axis produces torsion, as verified by cases 116/117.

Adding `K6ROT` to a CROD-only deck would therefore have no effect and could
obscure the reason for the `456` constraints. The five decks in cases
123--127 also include an explicit comment stating this point.

### Increment and iteration control

Example with fixed increments:

```text
NLSTEP,30,1.0,MILDLY
,GENERAL,30,1,10
,FIXED,20
,MECH,PV,1.E-8,1.E-8,1.E-8,PFNT
```

Adaptive example:

```text
NLSTEP,30,1.0,SEVERELY
,GENERAL,50,1,20
,ADAPT,.025,1.E-6,.10,4,1.20,0,500
,MECH,PV,1.E-9,1.E-9,1.E-9,PFNT
```

The subset read during preprocessing is:

- header: `ID`, total time, and `CTRLDEF` (`QLINEAR`, `MILDLY`, or
  `SEVERELY`);
- `GENERAL`: `MAXITER`, `MINITER`, and `MAXBIS`;
- `FIXED`: the number `NINC` of equal increments;
- `ADAPT`: initial/minimum/maximum fractions, desired iterations, growth
  factor, output intervals, and maximum number of increments;
- `MECH`: `U/P/W/V/N/A` criteria, tolerances, and stiffness-update policy.

`PFNT` requests full Newton-Raphson and is the option used throughout the
campaign. For `ADAPT`, `INTOUT=0` requests output at every converged increment
and `-1` only at the end; diagnostic cases use `0` to populate the
intermediate OP2 files requested by `NLOPRM`.
`ITER,KSTEP=1` is also normalized, but is not used as the primary reference.
The `ARCLN` arc-length methods, thermal control, and coupling are rejected
because they are outside the defined scope.

The implemented solver requires `PFNT` and requires `MECH` to include the
force-residual criterion `P`. Additional flags are retained in the model, but
increment acceptance is decided by this residual; input without `P` is
rejected instead of silently ignoring its criterion.

The reader also accepts legacy `NLPARM` input and converts it to the same data
model. New decks use `NLSTEP`, which is the preferred SOL 400 control.

### Output and diagnostics

```text
NLOPRM OUTCTRL = (SOLUTION,INTERM) NLDBG = ADVDBG
DISPLACEMENT(PRINT,PUNCH) = ALL
SPCFORCES(PRINT,PUNCH) = ALL
GPFORCE(PRINT,PUNCH) = ALL
STRESS(PRINT,PUNCH) = ALL
FORCES(PRINT,PUNCH) = ALL
```

`NLOPRM` must appear only once per control case; its options are placed next
to one another in the same statement. `NLOPRM` and the requests above do not
change the physics. They preserve converged states, the iteration/residual
table, and grid-point balances. `NLOUT` is normally reserved for selecting
advanced nonlinear results and is not required for the conventional CROD
tables in this campaign.

#### Locating final results in the F06

A SOL 400 F06 repeats its tables for every converged state. Before
transcribing an end-to-end reference, locate the final
`LOAD STEP = 1.00000E+00` block and confirm the message
`JOB CONVERGES FOR THE CURRENT STEP`. Within that block:

- `DISPLACEMENT VECTOR` contains `T1 T2 T3 R1 R2 R3` for each `POINT ID`;
- `FORCES OF SINGLE-POINT CONSTRAINT` contains the reactions in the same six
  columns;
- `NONLINEAR STRESSES IN ROD ELEMENTS (CROD)` contains `AXIAL STRESS`,
  `TOTAL STRAIN`, and `LIN. TORSIONAL STRESS` for each `ELEMENT ID`;
- `GLOBAL WORKS & ENERGIES` contains `TOT STRAIN ENERGY`;
- the `FORCE RESIDUAL CONVERGENCE RATIO` lines and convergence messages allow
  the iterations, residuals, and number of increments to be audited.

For example, in the final block for case 112, GRID 3 reports
`T2=-1.538383E+01`; GRID 1 reports reactions `T1=1.181807E+06` and
`T2=1.000000E+06`; and CROD 1 reports `AXIAL STRESS=-1.548118E+04` and
`TOTAL STRAIN=-7.371989E-02`. These values become the `u`, `r`, and `s` keys,
respectively, in `test_nonlinear_crod_validation_end_to_end.py`. Axial force is
obtained as `N=A*sigma` when it is not reported directly, and torque as
`T=(J/C)*tau_torsion`. The test module contains the complete example and the
mapping for all its keys.

There is no standard output equivalent to a DMIG for **every** Newton
tangent. The initial stiffness audit therefore relies on:

1. the 92/93 pair for the small-load limit;
2. `F_ext - F_int` equilibrium through `GPFORCE`, reactions, and CROD forces;
3. convergence factors and the worst residual from `.f06`/`.sts`;
4. sensitivity to increment discretization in 94/95 and 104/109/110.

If these layers are insufficient, a dedicated DMAP extraction must be
prepared for the MSC Nastran version in use; a linear DMIG must not be assumed
to represent the converged tangent.

The MSC Nastran 2026.1 audit also confirmed that 114/115 can be used as a
numerical probe: the central difference in displacements under `+/-1000`
horizontal loads approximates the tangent compliance around the deformed
state of 112. Cases 116/117 separate pure torsion from combined tension and
torsion; 119 defines a stable reference before the jump observed in 102; 120
deliberately exercises adaptive control through that jump; and the 121/122
pair verifies the linear limit in a structure whose orientation changes.

## Primary references

- [MSC Nastran 2021.1 Nonlinear User's Guide SOL 400](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2021.1_Nonlinear_SOL_400_User_Guide/raw/resource/enus/MSC_Nastran_2021.1_Nonlinear_SOL_400_User_Guide.pdf): analysis types, `NLSTEP`, convergence, and `LGDISP`.
- [MSC Nastran 2021.3 SOL 400 Getting Started Guide](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2021.3_SOL_400_Getting_Started_Guide/raw/resource/enus/MSC_Nastran_2021.3_SOL_400_Getting_Started_Guide.pdf): 1D elements, diagnostics, and intermediate `NLOPRM` output.
- [MSC Nastran 2020 Reference Guide](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2020_Reference_Guide/raw/resource/enus/MSC_Nastran_2020_Reference_Guide.pdf): geometrically nonlinear `SOL 400` examples with `PARAM,LGDISP`.
- [MSC Nastran 2023.4 Reference Guide](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2023.4_Reference_Guide/raw/resource/enus/MSC_Nastran_2023.4_Reference_Guide.pdf): scope of `K6ROT` as numerical stiffness for shell drilling rotation.

## Preprocessing contract

`FEMModel.analysis` contains an `AnalysisConfiguration` with:

- solution sequence and analysis type;
- active `ANALYSIS` keyword;
- `LGDISP` value, derived large-displacement activation, and follower-load
  stiffness;
- `NonlinearStepControl` normalized from `NLSTEP` or `NLPARM`.

pyNastran 1.4.1 does not implement the `NLSTEP` Bulk Data entry. It preserves
it in `reject_lines`; `nastran_reader.py` interprets its free-field or
small-field rows there and applies strict validation. An unknown continuation
is not ignored: it stops preprocessing.

## Implemented formulation

At every load-factor increment `lambda`, the solver seeks:

```text
r(u, lambda) = lambda F_ext - F_int(u) = 0
K_T(u) Delta_u = r(u, lambda)
u_(i+1) = u_i + Delta_u
```

For a current chord `d`, length `l`, and direction `n`, it uses:

```text
epsilon = (l-L0)/L0
N = E*A*epsilon
f_int = [-N*n, +N*n]
k = E*A/L0*(n outer n) + N/l*(I - n outer n)
```

The first term is material and the second is geometric. The element test
compares the analytical matrix with the central derivative of `f_int`; the
global test compares the deformed compliance `K_T^-1[ux,ux]` with the Nastran
114/115 probe. Torsion uses `GJ/L0` about the reference axis, reproduces `R1`,
torque, and `tau=T*C/J` from 116/117, and remains uncoupled from axial
extension.

## Validation status

The 35 cases that reach full load (92--119 and 121--127) are solved end to end
and compared with explicit Nastran values in
`test_nonlinear_crod_validation_end_to_end.py`. They cover displacements,
reactions, axial stress/strain, torsion, CD/CID systems, `LOAD`, 2D/3D
assemblies, prescribed displacements, and fixed/adaptive
controls. New entries use descriptive keys (`case_name`, `folder`, `filename`,
`description`, tolerances, and expected results) documented in the test module
itself.

The history for 112 matches at 25, 50, 75, and 100 percent load. The
horizontal tangent compliance obtained at the final state is `6.127186e-6`,
equal to the central difference from 114/115. Adaptive controls 95, 104, and
111 reproduce 9, 14, and 15 increments, respectively.

The project's `.rpt` report preserves this type of traceability during a
normal run. It identifies SOL 400 in the header, prints the normalized
configuration, and separates the converged-increment history from the Newton
correction history. Displacements, reactions, and CROD states under
`FINAL CONVERGED STATE RESULTS` belong exclusively to the final converged
equilibrium; the two histories follow. The [`result_report.md`](result_report.md)
guide explains every column and its location in the F06.

Structural cases 123--127 reproduce `25/30/30/40/14` increments; the complete
adaptive load-factor sequence for 127 also matches. The complete audit of
their F06 files obtained maximum relative errors of `3.77e-7` in
displacements, `3.76e-7` in reactions, and `3.72e-7` in CROD results. In 125,
MSC jumps from `T3=-13.5215` to `-69.40416` between factors `0.533333` and
`0.566667`; its accumulated `TotStrEn` is therefore not used as final elastic
strain energy, although the entire final mechanical state matches.

Case 120 does not reach full load under force control: both MSC and the solver
stop near `lambda=0.98387` while approaching the limit point. The exact number
of cutbacks is not matched because MSC's internal adaptation heuristic is not
published; the phenomenon and the final stable state do match. Arc-length
control is not implemented.

The energy returned by the project is the elastic strain energy of the current
state. Along continuous branches it agrees with Nastran apart from incremental
integration error. In 102/105, Nastran's global `TotStrEn` field accumulates
work across the jump and does not represent the same quantity; it is therefore
not used as a final reference in those two cases.
