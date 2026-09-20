# Geometrically nonlinear CBAR analysis

## Status and scope

`NonlinearStaticSolver` supports `MAT1/PBAR/PBARN1/CBAR` models with large
displacements and rotations. The implementation lives in
`src/solver/elements/cbar_nonlinear.py` and does not modify the linear
`BarElementSolver` path. The external CBAR campaign comprises MSC Nastran
SOL 400 cases 128--166.

The main functional reference is the official
[MSC Nastran SOL 400 2022.4 documentation](https://documentation-be.hexagon.com/bundle/MSC_Nastran_2022.4_SOL_400_Getting_Started_Guide/raw/resource/enus/MSC_Nastran_2022.4_SOL_400_Getting_Started_Guide.pdf?save_local=true).
For CBAR/PBAR, SOL 400 supports large displacements and rotations but not large
section strains. Forces and stresses are expressed in axes that corotate with
the element.

The implemented scope is:

- isotropic elastic `MAT1` and a constant `PBAR` section;
- `PBARN1,SECT=S,BEH2=BAR`, with validated `INT2=LC` and `INT2=LS`;
- axial, torsional, and bending response in both planes, including `I12`;
- linear/cubic Euler--Bernoulli interpolation for `LC` and one-point linear
  shear interpolation for `LS`;
- transverse-shear factors `K1/K2` in `LS`, including the unit value that
  SOL 400 applies when the field is zero;
- `WA/WB` offsets, `PA/PB` releases, and vector or `G0` orientation;
- `FORCE`, `MOMENT`, `SPC`, `SPC1`, `SPCADD`, and prescribed values;
- conservative loads (`FOLLOWK=NO`) and full Newton (`PFNT`).

Plasticity, large section strains, warping, distributed or follower loads,
PBEAM, CBARAO, and arc-length control are outside the implemented scope.

## Corotational kinematics

Each element stores its reference length and triad. Total nodal rotations are
interpreted as rotation vectors and converted to matrices through the
exponential map. Offsets rotate with their nodes:

```text
xa = Xa + ua + Ra*WA
xb = Xb + ub + Rb*WB
l  = ||xb-xa||
ex = (xb-xa)/l
```

The current triad is constructed in two steps. First, the shortest rotation
that maps the reference axis to the current chord is applied. The average
physical twist of both sections about the chord is then added. This
*swing--twist* decomposition prevents bending rotations from contaminating
the section twist.

Nodal orientations relative to this triad provide the local end rotations.
The basic deformation vector is:

```text
d = [L0*log(l/L0), theta_xb-theta_xa,
     theta_ya, theta_za, theta_yb, theta_zb]
```

The first term is the logarithmic strain integrated over the reference
length. Axial cases 130/131 identify it directly: for a stress `sigma`,
Nastran returns `l/L0=exp(sigma/E)`. By contrast, `LGDISP=-1` preserves the
exact reference linear formulation.

## PBARN1 interpolation

The official SOL 400 documentation distinguishes two two-node CBAR elements.
`INT2=LC` uses linear interpolation along the axis and cubic interpolation for
transverse displacement, with three integration points. This is the
Euler--Bernoulli branch: the `PBAR K1/K2` factors do not change its stiffness.
Cases 154 and 160, which differ only in those factors, produce exactly the
same MSC result.

`INT2=LS` uses linear transverse displacement and section rotation with one
integration point. In the local XY plane, for the DOFs
`[vA, theta_zA, vB, theta_zB]`, the implementation uses:

```text
B_curvature = [0, -1/L, 0, 1/L]
B_shear     = [-1/L, -1/2, 1/L, -1/2]

K_LS = E*I*L*(B_curvature outer B_curvature)
     + K*A*G*L*(B_shear outer B_shear)
```

The second plane applies the CBAR-specific sign change for `theta_y`. The
`I12` coupling acts between the two constant curvatures. Releases and offsets
are then applied with the same conventions as in `LC`.

For a cantilever with a transverse force at its tip, the two interpolations
give, in their linear limit:

```text
delta_LC = F*L^3/(3*E*I)
delta_LS = F*L^3/(4*E*I) + F*L/(K*A*G)
theta    = F*L^2/(2*E*I)
```

The `1/4` versus `1/3` difference explains why replacing `LS` with the closed
Timoshenko matrix of the linear CBAR did not reproduce case 142. Pairs 155/156
and 157/158 isolate both planes. Case 159 demonstrates that SOL 400 uses
`K=1` when `K1/K2` are zero in `LS`; those fields are ignored in `LC`.

## Section law, internal force, and tangent

The `LC` branch reuses the Euler matrix from `BarElementSolver` with `K1/K2`
explicitly set to zero. The `LS` branch constructs the matrix above inside
`CBARNonlinear`. The six basic measures are embedded in a 12-DOF local vector,
and the actions are computed as:

```text
q_local = K_local d_local
```

Forces and moments are then rotated into basic coordinates. Offset moments add
`r x f`. This spatial vector is assembled and balanced against the nodal
loads.

The tangent is the central derivative of the complete internal force with
respect to the 12 total displacement and rotation-vector coordinates. It
therefore includes changes in the chord, triad, offsets, and geometric
actions. Symmetry is not imposed: with spatial moments and total rotation
vectors, these incremental coordinates do not form a conjugate pair that
guarantees a symmetric matrix. `spsolve` solves the nonsymmetric system
directly.

This choice prioritizes an auditable physical force and a consistent
derivative. Its cost is high: one tangent evaluation requires 24 CBAR force
evaluations. Any future analytical linearization must first reproduce the
finite-difference and Nastran-probe tests.

## Recovery

`CBARNonlinearState` contains:

- the 12x12 internal-force tangent and internal force;
- current length and triad;
- basic/local deformations and local end actions;
- axial strain, stress, and force, plus twist and torque;
- normal stresses at recovery points C--F at both ends;
- stress extrema and state strain energy.

Axial stress uses `E*log(l/L0)`. Bending reuses the linear PBAR recovery in the
corotated triad. Centroidal axial stresses are compared in the end-to-end
tests. The exact MSC convention for fibre stresses under combined 3D rotation
is not fully resolved; a centroidal stress must not be inferred by averaging
C--F when the recovery points are not symmetric.

### Presentation in the `.rpt` report

The final converged state of each CBAR is presented in four blocks:

1. `SECTION RESPONSE`: axial force and stress, torque, maximum torsional stress
   at the configured points, logarithmic axial strain, and twist per unit
   length;
2. `NONLINEAR STATE`: current offset-end length, element strain energy, and
   maximum/minimum total normal stress across both ends;
3. `END FORCES AND MOMENTS`: local `Fx/Fy/Fz/Mx/My/Mz` actions at A and B;
4. `STRESSES AT RECOVERY POINTS`: axial, bending, total, and torsional
   contributions at C--F for both ends.

Dictionary keys `P0`--`P3` are converted to C--F only during presentation;
the numerical result is not changed. The `NONLINEAR STATE` table is enabled
by the joint presence of `current_length` and `strain_energy`, so a linear
CBAR retains the other three tables without showing a fictitious nonlinear
state. The detailed F06 correspondence is documented in
[`result_report.md`](result_report.md).

Decks 161--166 isolate this convention progressively and request C--F results
at final and intermediate states. Their F06/PCH files end without `FATAL`,
reach load factor 1.0, and are part of the CBAR end-to-end suite. The
centroidal axial stress fitted from C--F agrees approximately within 0.01% to
1.7%; fibres away from the neutral axis are also compared. Some nearly
neutral fibres and the moment distribution in the second segment of case 165
and third segment of case 166 remain diagnostics: they are not hidden behind
broad tolerances or `xfail`.

Mixed models 167--176 verify that CBAR recovery and assembly coexist with the
corotational axial CROD response. Their ten references, together with the
updated final bridge, form an independent 11-case end-to-end suite and reach
full load. See
[`NONLINEAR_FINAL_CLOSURE_CASES.md`](../tests/solver_tst/validation_cases/NONLINEAR_FINAL_CLOSURE_CASES.md).

## Validation

The tests are organized into five complementary groups:

- `test_cbar_nonlinear.py`: ten test functions, collecting eleven checks of
  the linear limit, `LC/LS` matrices, the unit shear default, rigid-rotation
  objectivity, logarithmic strain, torsion, releases, offsets, and the tangent
  derivative;
- `test_nonlinear_cbar_intermediate_states.py`: paths for cases 132/137,
  tangent compliance from 138/139, and states 0.2/0.5/0.8 for the `LC/LS`
  diagnostics 155--159;
- `test_nonlinear_cbar_validation_end_to_end.py`: an explicit dictionary of
  the 39 cases, with displacements, reactions, axial forces, torque, selected
  increment checks, and C--F fibres for 161--166;
- `test_nonlinear_mixed_validation_end_to_end.py`: eleven integrated
  references--cases 167--176 and the updated final bridge--with displacements,
  reactions, and CBAR/CROD results;
- `test_nonlinear_final_bridge.py`: structural checks for the updated
  three-span bridge, including forces/moments, imposed settlements,
  equilibrium, and a signal relative to `LGDISP=-1`.

The `dT2/dFy` compliance in case 137 differs from the MSC probe by 0.1% to
0.8% at load factors 0.2--1.0. All 39 CBAR and 11 mixed end-to-end references
reach full load. Tolerances reflect the complexity: basic cases remain below
0.1%, larger assemblies are compared between 0.3% and 4%, and released frame
147 uses 7% because of a secondary transverse reaction.

Diagnostic cases 154--160 isolate `LC/LS`, each shear plane, and the presence
or absence of `K1/K2`. All converge without `FATAL` and are part of the final
end-to-end suite. Against MSC, the maximum vector errors across these seven
cases are `8.04e-5` in displacements, `6.03e-5` in reactions, and `1.52e-4` in
energy. Case 142 passes with `rtol=5e-4` and is no longer an `xfail`.

The DAT files and interpretation of their outputs are documented in
[`NONLINEAR_CBAR_CASES.md`](../tests/solver_tst/validation_cases/NONLINEAR_CBAR_CASES.md)
and [`validation.md`](validation.md). The sole final integration model is the
updated three-span bridge documented in its
[`README`](../tests/solver_tst/validation_cases/case_nonlinear_test_final_updated/README.md).
It contains 202 GRID, 144 CBAR, 428 CROD, and explicit foundations. Its full
correlation retains L2 errors of order `10^-5` in global and axial response;
continuous-girder C--F recovery is quantified separately (`1.673%` L2 and
`11.57%` maximum significant error), without applying that envelope to other
fields.
