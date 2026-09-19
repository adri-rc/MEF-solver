# Mathematical formulation

## General assumptions

The linear solver covers small-deformation statics with isotropic materials:

```text
K u = F
```

Every node reserves six degrees of freedom
`[Tx, Ty, Tz, Rx, Ry, Rz]`. Each element first builds stiffness in its local
frame and transforms it as:

```text
u_local  = T u_global
K_global = T.T K_local T
```

The global matrix is assembled sparsely. There is no general `AUTOSPC`: a
mechanism or a free DOF without stiffness produces a singularity error.

## Isotropic material

`MAT1` accepts any two of `E`, `G`, and `nu` and completes the third
from the isotropic relation:

```text
G = E / (2*(1+nu))
```

Shell membrane and bending use plane stress:

```text
D = E/(1-nu²) * [[1, nu, 0],
                  [nu, 1, 0],
                  [0, 0, (1-nu)/2]]
```

## CROD

CROD is a straight two-node rod. In its local frame:

```text
K_axial   = E*A/L * [[1,-1],[-1,1]]
K_torsion = G*J/L * [[1,-1],[-1,1]]
```

Recovery returns axial strain/stress and force, twist per unit length, torque,
and torsional stress.

### Geometric nonlinear CROD

`CRODNonlinear` uses reference positions `X1`, `X2` and total
displacements to form the current chord:

```text
d       = (X2 + u2) - (X1 + u1)
l       = ||d||
n       = d/l
L0      = ||X2-X1||
epsilon = (l-L0)/L0
sigma   = E*epsilon
N       = A*sigma
```

The translational internal force in basic coordinates is
`[-N*n, +N*n]`. Its consistent derivative contains material and geometric
parts:

```text
k   = (E*A/L0) (n outer n) + (N/l) (I - n outer n)
K_e = [[ k, -k],
       [-k,  k]]
```

The second term becomes negative in compression and is required for stiffness
loss, large rotation, and snap-through. Area and modulus remain referenced to
the initial state; the implementation has neither finite material strain nor
an isochoric section update.

At every increment, full Newton solves:

```text
r(u, lambda) = lambda F_ext - F_int(u)
K_T(u) Delta_u = r(u, lambda)
```

Torsion remains uncoupled:
`T=G*J*(theta2-theta1).n0/L0`. Cases 116/117 confirm that axial extension
does not alter this reference torsional stiffness.

`LGDISP=-1` retains reference directions and stiffness; `LGDISP=1` or
`2` enables the current-chord formulation. With supported conservative loads
and `FOLLOWK=NO`, values 1 and 2 give the same response in cases 112/118.
See [nonlinear_crod.md](nonlinear_crod.md).

## CBAR

Linear CBAR includes axial extension, torsion, and bending in two planes. Its
reference triad is defined by the chord and orientation vector or `G0`.
`PBAR` supplies `A`, `I1`, `I2`, `I12`, `J`, and `K1/K2`.
Offsets and end releases are included.

### Geometric nonlinear CBAR

`CBARNonlinear` removes rigid motion with a corotated frame. Nodal rotations
use the exponential map. The frame follows the chord through the shortest
rotation and adds the mean section twist about it; offsets rotate with their
nodes.

Cases 130/131 identify logarithmic axial strain:

```text
epsilon_log = log(l/L0)
sigma       = E*epsilon_log
N           = A*sigma
```

`L0*epsilon_log` and relative end rotations enter the validated local PBAR
matrix. Local actions rotate back to the basic frame, and offsets add
`r x f`. The tangent is a central derivative of the complete spatial force
and may be nonsymmetric when Newton is parameterized by total rotation
vectors.

`LGDISP=-1` retains reference kinematics; `LGDISP=1/2` enables the
corotational formulation. `INT2=LC` uses linear axial and cubic transverse
Euler--Bernoulli interpolation, independent of `K1/K2`. `INT2=LS` uses
linear displacement/rotation, constant curvature, and midpoint shear with
transverse stiffness `K1*A*G` and `K2*A*G`. See
[nonlinear_cbar.md](nonlinear_cbar.md).

## CTRIA3

CTRIA3 combines a constant-strain membrane with a Reissner--Mindlin plate.
Its physical local DOFs are `[u,v]` for membrane and
`[w,theta_x,theta_y]` for the plate. When `K6ROT>0`, normal rotation is
stabilized by penalizing its difference from membrane macro-rotation:

```text
omega = 0.5*(dv/dx - du/dy)
```

Thickness multiplies the membrane constitutive matrix. The plate applies the
two PSHELL factors independently:

```text
D_b = (12I/T^3) * E2*t^3/[12*(1-nu2^2)] * C(nu2)
D_s = (TS/T) * G3*t * I
```

`TS/T` therefore replaces the default `5/6` shear-correction factor; it is
not the `12I/T^3` factor, which scales bending only. Fibre-stress recovery
uses the effective section modulus `(12I/T^3)*t^2/6`, matching MSC output
for non-default factors.

The three PSHELL material fields do not inherit one another:

- blank `MID1` disables membrane;
- blank `MID2` disables the complete plate block;
- `MID3` without `MID2` is rejected, matching MSC Nastran `FATAL 9994`;
- with active `MID2` and blank `MID3`, bending remains but the explicit
  physical Reissner--Mindlin shear term is zero.

This matches material selection and inactive-block semantics; it does not make
the local plate identical to MSC's native CTRIA3. MSC retains transverse
response with active `MID2` and blank `MID3`, whereas the explicit local
shear term is zero. The DMIG references validate block activation, not
complete equality of plate stiffness.

Centroidal three-dimensional CTRIA3 moments/shears are not considered
equivalent to MSC within the current scope. The physical T-frame correlates
global response, membrane, and fibre stresses but not `M/Q`. A diagnostic
with `E(MID3)=1e-20` reduces moment error to `7.106e-5` and makes shear
engineering zero in both solvers, while creating an unusable near-mechanism.
This isolates the limitation without validating the degenerate model or
justifying a stiffness recalibration. See
[validation.md](validation.md#complementary-diagnostic-with-negligible-mid3).

## CQUAD4

CQUAD4 has four nodes and 24 DOFs. The active formulation reconstructs
Nastran's QUAD4 behaviour for **coplanar** geometry. Its local frame uses:

- origin at the mean of the four nodes;
- normal from diagonal cross product `G1-G3 x G2-G4`;
- local x along the bisector of unit diagonals `G1-G3` and `G4-G2`;
- local y completing the right-handed triad.

Bilinear shape functions are defined on `[-1,1]²`. Membrane and bending
combine the uniform central strain with two normalized bilinear stabilization
modes:

```text
K0  = A*s*B0.T*D*B0
h   = h0 - G0.T*X.T*h0
chi = 1/(1 - ||h-h0||²/12)
M   = A*(a outer a/||a||⁴ + b outer b/||b||⁴)
Khg = s*L*(E*chi*M/12)*L.T
```

Here `s=t` for membrane and
`s=(t³/12)*(12I/T³)` for bending.

Transverse shear uses MacNeal edge strings sampled along both natural
midlines. It is expressed through augmented flexibility:

```text
Z = Z_residual_bending + Z_MID3
K_shear = S.T * solve(Z, S)
```

`MID3` contributes flexibility in series. Blank `MID3` does not inherit
`MID2`, and blank `MID2` does not inherit `MID1`. The aspect-ratio
correction uses `epsilon=0.025`.

For skew affine elements, the torsional mode uses:

```text
m(c²) = 1 + 0.4100259495463487*c² /
             (1 - 0.9156556465380681*c²)
```

The function was identified from cases 81/89/90 and confirmed by loaded case
60; case 91 exercises two highly skew panels.

CQUAD4 extends the CTRIA3 drilling penalty to all four corners:

```text
r_i = theta_z_i - omega_i
K6  = sum_i K6ROT*1e-6*G*t*(4*detJ_i)*(r_i outer r_i)
```

All nodes must satisfy coplanarity tolerance
`max(1e-10, 1e-8*L)`. Warping is rejected rather than projected away.
`SNORM` smoothing is not implemented; faceted assemblies must use
independent normals equivalent to `SNORM,0`.

Recovery at the centre returns membrane forces, moments, shears, top/bottom
fibre stresses, principal stresses, and von Mises values. Centroidal
`Qx/Qy` may differ from MSC's reporting convention; stiffness,
displacements, reactions, membrane forces, moments, and fibres have
independent validation.

The complete derivation is in
[cquad4_formulation.md](cquad4_formulation.md).

## Boundary conditions and reactions

Free DOFs `f` and prescribed DOFs `s` are partitioned as:

```text
K_ff u_f = F_f - K_fs u_s
```

This supports nonzero SPC values. Conflicting values on the same DOF raise an
explicit error. Reactions use the original assembled system:

```text
R = K u - F
```
