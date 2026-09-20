# Planar CQUAD4: final formulation and validation

## Status

The active implementation is `CQUAD4NastranElementSolver`. It reconstructs
the Nastran QUAD4 stiffness for **coplanar**, non-degenerate quadrilaterals
with a positive Jacobian. The plane may have any three-dimensional
orientation, and different elements in one structure may lie in different
planes.

Warped CQUAD4 elements are deliberately outside the project scope. The base
class checks the distance of all four nodes from the local plane and raises
`ValueError` if their range exceeds `max(1e-10, 1e-8*L)`, where `L` is the
greatest side length. The geometry is neither projected nor corrected for
warping.

The complete-suite baseline observed on 19 September 2026 was `714 passed,
6 xfailed, 2 warnings`. All six `xfail` cases concern this formulation: warped
cases 54, 59, and 62--64, and case 78 with default `SNORM` smoothing. The
equivalent faceted cases with `SNORM,0` pass. The two warnings were local
`PytestCacheWarning` messages, not numerical failures.

## Classes

- `CQUAD4ElementBase`: planar geometry, shape functions, constitutive
  matrices, transformation, `K6ROT`, and shared recovery.
- `CQUAD4NastranElementSolver`: the formulation used by the global solver.
- `CQUAD4ReissnerMindlinElementSolver`: an independent classical formulation
  retained for academic comparison; it is not the active element.

## Geometry and local system

The origin is the average of the four nodes. The normal is obtained from the
cross product of diagonals `G1-G3` and `G2-G4`; local x is the bisector of the
unit diagonals `G1-G3` and `G4-G2`; and local y completes the right-handed
triad.

The bilinear shape functions on `[-1,1]^2` are:

```text
N1 = 1/4 (1-xi)(1-eta)       N2 = 1/4 (1+xi)(1-eta)
N3 = 1/4 (1+xi)(1+eta)       N4 = 1/4 (1-xi)(1+eta)
```

The Jacobian is evaluated at every point and vertex required by the
formulation. A non-positive determinant is rejected.

## PSHELL

- `MID1`: membrane material; it is required by the active CQUAD4.
- `MID2`: bending material; a blank value disables bending and does **not**
  inherit `MID1`.
- `MID3`: material for physical transverse-shear flexibility; a blank value
  does not inherit `MID2`. QUAD4's residual bending flexibility remains.
- `12I/T^3`: multiplies `t^3/12`, with a default value of 1.
- `TS/T`: multiplies the shear thickness, with a default value of `5/6`.
- `THETA`: orients the material; with isotropic `MAT1` it does not change the
  constitutive matrix, and output is recovered in element axes.

A bending-only CQUAD4 with blank `MID1` is outside the supported subset and is
rejected explicitly.

## Membrane and bending

The Nastran formulation does not use full 2x2 quadrature. Each block combines:

1. a uniform strain evaluated at the center;
2. two bilinear `hourglass` modes orthogonalized against linear fields;
3. distortion normalization;
4. a metric built from lines joining the midpoints of opposite sides.

In abbreviated form:

```text
K0  = A*s*B0.T*D*B0
h   = h0 - G0.T*X.T*h0,  h0 = [1,-1,1,-1]
chi = 1 / (1 - ||h-h0||²/12)
M   = A*(a⊗a/||a||⁴ + b⊗b/||b||⁴)
Khg = s*L*(E*chi*M/12)*L.T
K   = K0 + Khg
```

Here `s=t` for membrane and `s=(t³/12)*(12I/T³)` for bending.

## Transverse shear

The kinematic operator is derived from the four edges, following the public
NASTRAN-95 `UEV/ANGLEI/EDGSHR` construction. Two samples are taken on each
natural centerline:

```text
(0, -1/sqrt(3)), (0, +1/sqrt(3)),
(-1/sqrt(3), 0), (+1/sqrt(3), 0)
```

The stiffness is obtained from an augmented flexibility:

```text
Z = Z_residual_bending + Z_MID3
Kshear = S.T * solve(Z, S)
```

`MID3` is added in series as flexibility, so activating it makes the element
softer. The aspect correction uses `epsilon=0.025`.

For oblique affine quadrilaterals, the antisymmetric mode of each pair is
scaled by:

```text
m(c²) = 1 + 0.4100259495463487*c² / (1 - 0.9156556465380681*c²)
```

Here `c` is the cosine between the midlines. The law equals 1 for a rectangle,
remains positive over the complete domain, and was identified from moderate-
obliquity DMIG cases 81/89 and high-obliquity case 90. Loaded case 60 checks
an intermediate angle. This is a project identification, not a constant
published by MSC.

## K6ROT

The implementation extends the CTRIA3 formulation, which has already been
validated, to four nodes. At every vertex it penalizes:

```text
r_i = theta_z_i - 0.5*(dv/dx - du/dy)_i
K6 = sum_i K6ROT*1e-6*G*t*(4*detJ_i)*(r_i outer r_i)
```

This is not an isolated diagonal stiffness. It preserves rigid-body motion
and couples nodal drilling rotation to the membrane macro-rotation.

## Recovery

Membrane forces `N`, moments `M`, top/bottom fiber stresses, principal
stresses, and von Mises stresses are evaluated at the center. These
quantities, together with stiffness, displacements, and reactions, closely
match Nastran.

The shell-element result interface uses one name for each fiber component:
`stress_top_fiber_x`, `stress_top_fiber_y`, `tau_top_fiber_xy`,
`stress_bot_fiber_x`, `stress_bot_fiber_y`, and `tau_bot_fiber_xy`. `top`
means `z=+t/2` and `bot` means `z=-t/2` in element axes. No `Z1/Z2`-based
aliases are exposed because that label depends on the external solver's
output convention.

`Qx/Qy` are recovered from the four generalized shear forces. MSC's
proprietary centroidal convention is not reproduced exactly: in case 91 the
largest difference is 47.1% on a small component, 0.0216 in absolute value,
while `N`, `M`, and fiber stresses remain below 0.1%. This is a recovery
limitation, not a stiffness or equilibrium error.

## Validation

`tests/solver_tst/elements/test_E2E_CQUAD_Nastran_K_solution.py` compares
DMIG data using:

- global relative Frobenius error;
- termwise comparisons with relative and absolute tolerances;
- element-axis membrane, plate/shear, and drilling blocks;
- drilling--membrane coupling where applicable;
- optional export through `debug_export` and `output_k_file`.

`tests/solver_tst/test_validation_end_to_end_cquad4.py` stores displacement,
load, reaction, and element-result references directly; it does not parse the
F06 during a test.

Case 90 required changing only the geometric continuation of the torsional
mode. With the final law, its full DMIG has a Frobenius error below `1e-6`, its
plate/shear block is below `4e-6`, and its terms satisfy `rtol=1e-5`,
`atol=2e-3`. Case 91 gives vector errors of `2.07e-5` for displacements and
`6.30e-5` for reactions.

Cases 57 and 58 are not clean formulation controls: their MID3 modulus is of
order `1e-20`, with condition numbers `2.90e25` and `8.52e26`. Their comments
and tolerances document which small components are sensitive. Case 59 is
warped and is retained only as a scope `xfail`.

## Limitations

- Coplanar, non-degenerate quadrilaterals only.
- Isotropic `MAT1`; no laminates, `MAT8`, or composite coupling.
- No `ZOFFS`, nodal thicknesses, `TFLAG`, `MID4`, or face loads.
- No `SNORM` smoothing; use independent normals, equivalent to `SNORM,0`, in
  faceted assemblies.
- `Qx/Qy` recovery is approximate relative to MSC Nastran output.

The mathematical derivation and detailed code correspondence are in
[`cquad4_formulation.md`](cquad4_formulation.md). See also the general
[`architecture.md`](architecture.md), [`formulation.md`](formulation.md), and
[`validation.md`](validation.md) documents.

## References

- R. H. MacNeal, *A Simple Quadrilateral Shell Element*, Computers &
  Structures 8, 175--183, 1978.
- R. H. MacNeal, [*Derivation of element stiffness matrices by assumed strain
  distributions*](https://doi.org/10.1016/0029-5493(82)90262-X), 1982.
- NASA, [NASTRAN-95 `quad4d.f`](https://github.com/nasa/NASTRAN-95/blob/master/mis/quad4d.f).
- Hexagon, *MSC Nastran Reference Guide*.
