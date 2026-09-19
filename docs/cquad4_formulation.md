# Planar CQUAD4: formulation, implementation, and scope

## 1. Purpose of this document

This document provides a self-contained description of the project's
four-node quadrilateral element. The explanation is derived from the active
Python code and validation cases and does not depend on any artifact outside
the repository's Markdown documentation.

The global solver uses `CQUAD4NastranElementSolver`, a reconstruction of the
Nastran QUAD4 structure for **planar** elements. The classical
Reissner--Mindlin formulation is retained in a separate class only as an
academic reference; it is not assembled by default.

The principal files are:

- `src/solver/elements/cquad4_element.py`: geometry, transformation,
  constitutive matrices, `K6ROT`, and shared recovery;
- `src/solver/elements/cquad4_nastran_element.py`: active Nastran-type
  stiffness;
- `src/solver/elements/cquad4_reissner_mindlin_element.py`: independent
  classical isoparametric formulation;
- `tests/solver_tst/elements/test_cquad4_element.py`: unit properties and
  geometric validation;
- `tests/solver_tst/elements/test_E2E_CQUAD_Nastran_K_solution.py`: DMIG
  matrix comparisons;
- `tests/solver_tst/test_validation_end_to_end_cquad4.py`: displacements,
  loads, reactions, and element results against Nastran.

For project-wide flow and coordinate conventions, see
[`architecture.md`](architecture.md). The general linear equations are in
[`formulation.md`](formulation.md), and the complete evidence map is in
[`validation.md`](validation.md).

## 2. Physical and geometric assumptions

The formulation assumes:

- small strains and small displacements;
- isotropic, linear-elastic `MAT1` material;
- plane stress for membrane and bending;
- four nodes ordered consecutively around the boundary;
- a non-degenerate isoparametric transformation with positive Jacobian;
- four nodes coplanar within the numerical tolerance.

The element may be inclined and oriented arbitrarily in three-dimensional
space. A structure may contain panels in different planes, provided every
individual CQUAD4 is planar.

### 2.1. Coplanarity restriction

If `x_i` are the nodal coordinates, `x_bar` is their average, `e_z` is the
local normal, and `L` is the greatest side length, the code checks:

$$
\max_i[(\mathbf{x}_i-\bar{\mathbf{x}})\cdot\mathbf{e}_z]
-\min_i[(\mathbf{x}_i-\bar{\mathbf{x}})\cdot\mathbf{e}_z]
\leq \max(10^{-10},10^{-8}L).
$$

If this condition is not satisfied, the constructor raises `ValueError`. The
solver neither projects the nodes onto a mean plane nor applies a warping
correction. This limitation is deliberate and is part of the final scope.

The Jacobian is also checked at all four natural vertices when the element is
constructed, independently of `K6ROT` and of which `PSHELL` blocks are active.
Concave, inverted, or cross-numbered quadrilaterals are therefore rejected
before assembly.

## 3. Local system and degrees of freedom

The local origin is the average of the four positions:

$$
\bar{\mathbf{x}}=\frac{1}{4}\sum_{i=1}^{4}\mathbf{x}_i.
$$

Define the diagonals:

$$
\mathbf{d}_{13}=\mathbf{x}_3-\mathbf{x}_1,
\qquad
\mathbf{d}_{24}=\mathbf{x}_4-\mathbf{x}_2,
\qquad
\mathbf{d}_{42}=\mathbf{x}_2-\mathbf{x}_4.
$$

The orthonormal basis used by the code is:

$$
\mathbf{e}_z=
\frac{\mathbf{d}_{13}\times\mathbf{d}_{24}}
     {\|\mathbf{d}_{13}\times\mathbf{d}_{24}\|},
$$

$$
\mathbf{e}_x=
\frac{\widehat{\mathbf{d}}_{13}+\widehat{\mathbf{d}}_{42}}
     {\|\widehat{\mathbf{d}}_{13}+\widehat{\mathbf{d}}_{42}\|},
\qquad
\mathbf{e}_y=\mathbf{e}_z\times\mathbf{e}_x.
$$

Both diagonals must be normalized before they are added. Using their original
lengths would make the local orientation depend improperly on the longer
diagonal of an irregular quadrilateral.

The planar coordinates of node `i` are:

$$
x_i=(\mathbf{x}_i-\bar{\mathbf{x}})\cdot\mathbf{e}_x,
\qquad
y_i=(\mathbf{x}_i-\bar{\mathbf{x}})\cdot\mathbf{e}_y.
$$

Each node reserves six degrees of freedom. The complete local vector has 24
components:

$$
\boldsymbol{\delta}^{e}=
[u_1,v_1,w_1,\theta_{x1},\theta_{y1},\theta_{z1},\ldots,
u_4,v_4,w_4,\theta_{x4},\theta_{y4},\theta_{z4}]^T.
$$

- `u`, `v`: membrane displacements;
- `w`, `theta_x`, `theta_y`: plate displacement and rotations;
- `theta_z`: normal rotation, or *drilling DOF*.

Membrane and plate are uncoupled in the local system. Coupling to `theta_z`
comes only from `K6ROT` stabilization.

## 4. Isoparametric interpolation

The natural domain is `(xi, eta) in [-1,1]^2`, with nodes:

```text
G1 = (-1,-1)    G2 = ( 1,-1)
G4 = (-1, 1)    G3 = ( 1, 1)
```

The bilinear functions are:

$$
N_1=\frac14(1-\xi)(1-\eta),
\qquad
N_2=\frac14(1+\xi)(1-\eta),
$$

$$
N_3=\frac14(1+\xi)(1+\eta),
\qquad
N_4=\frac14(1-\xi)(1+\eta).
$$

The geometry and any bilinear nodal field are interpolated as:

$$
x(\xi,\eta)=\sum_{i=1}^{4}N_i x_i,
\qquad
y(\xi,\eta)=\sum_{i=1}^{4}N_i y_i.
$$

The Jacobian used in Python is:

$$
\mathbf{J}=
\begin{bmatrix}
\partial x/\partial\xi & \partial y/\partial\xi\\
\partial x/\partial\eta& \partial y/\partial\eta
\end{bmatrix}.
$$

Physical derivatives are evaluated through:

$$
\begin{bmatrix}N_{i,x}\\N_{i,y}\end{bmatrix}
=\mathbf{J}^{-1}
\begin{bmatrix}N_{i,\xi}\\N_{i,\eta}\end{bmatrix}.
$$

The formulation requires `det(J) > 0` at every point it uses.

## 5. Material, constitutive matrix, and PSHELL semantics

For an isotropic material, the plane-stress matrix is:

$$
\mathbf{D}(E,\nu)=\frac{E}{1-\nu^2}
\begin{bmatrix}
1&\nu&0\\
\nu&1&0\\
0&0&(1-\nu)/2
\end{bmatrix}.
$$

`PSHELL` separates the physical mechanisms:

- `MID1`: membrane material;
- `MID2`: bending material;
- `MID3`: material for physical transverse-shear flexibility;
- `12I/T^3`: bending-inertia ratio;
- `TS/T`: effective shear-thickness ratio.

The active CQUAD4 requires `MID1`. A bending-only property with blank `MID1`
is outside the supported subset and is rejected explicitly. The remaining
blank-field semantics are also significant:

- blank `MID2` disables bending and does **not** inherit `MID1`;
- blank `MID3` disables physical shear flexibility and does **not** inherit
  `MID2`;
- even when `MID3` is blank, the QUAD4 residual bending flexibility remains
  when `MID2` is present.

The section factors are:

$$
A_m=t,
\qquad
I_b=\frac{t^3}{12}r_I,
\qquad
t_s=t\,r_s,
$$

where:

$$
r_I=\frac{12I}{t^3},
\qquad
r_s=\frac{TS}{t}.
$$

The default values are `r_I=1` and `r_s=5/6`.

The `THETA` field orients the material axes. With isotropic `MAT1`, the
constitutive matrix is invariant under that rotation. CQUAD4 forces and
stresses are recovered in the element's own axes, matching the Nastran output
used for validation.

## 6. Strain--displacement matrices

### 6.1. Membrane

For the vector:

$$
\mathbf{u}_m=[u_1,v_1,u_2,v_2,u_3,v_3,u_4,v_4]^T,
$$

the membrane matrix is:

$$
\mathbf{B}_m=
\begin{bmatrix}
N_{1,x}&0&\cdots&N_{4,x}&0\\
0&N_{1,y}&\cdots&0&N_{4,y}\\
N_{1,y}&N_{1,x}&\cdots&N_{4,y}&N_{4,x}
\end{bmatrix}.
$$

Therefore:

$$
\boldsymbol{\varepsilon}_m=\mathbf{B}_m\mathbf{u}_m,
\qquad
\boldsymbol{\sigma}_m=\mathbf{D}_m\boldsymbol{\varepsilon}_m.
$$

### 6.2. Plate curvature

For:

$$
\mathbf{q}_p=[w_1,\theta_{x1},\theta_{y1},\ldots,
w_4,\theta_{x4},\theta_{y4}]^T,
$$

the convention is:

$$
\boldsymbol{\kappa}=
\begin{bmatrix}
\partial\theta_y/\partial x\\
-\partial\theta_x/\partial y\\
\partial\theta_y/\partial y-\partial\theta_x/\partial x
\end{bmatrix}
=\mathbf{B}_b\mathbf{q}_p.
$$

The bending stiffness is first built as an eight-degree-of-freedom planar
block and is then transformed to `[w, theta_x, theta_y]` order through the
rotation mapping used in the code.

## 7. Nastran-type membrane and bending stiffness

The active implementation does **not** obtain these blocks with full `2x2`
Gauss integration. Membrane and curvature each combine a central uniform
strain with two stabilized bilinear modes.

Let:

- `B0 = B(0,0)`;
- `J0 = det(J(0,0))`;
- `A = 4 J0`;
- `s=t` for membrane;
- `s=I_b` for bending.

The uniform contribution is:

$$
\mathbf{K}_0=A\,s\,\mathbf{B}_0^T\mathbf{D}\mathbf{B}_0.
$$

### 7.1. Bilinear mode

Start from:

$$
\mathbf{h}_0=[1,-1,1,-1]^T.
$$

If `G0` contains physical shape-function gradients at the center and `X`
contains local nodal coordinates, remove the part represented by linear
fields:

$$
\mathbf{h}=\mathbf{h}_0-\mathbf{G}_0^T\mathbf{X}^T\mathbf{h}_0.
$$

The modal operator `L` places `h` in the two in-plane directions. Mode energy
must be corrected when the diagonals no longer bisect each other. The
normalization is:

$$
\chi=\left(1-\frac{\|\mathbf{h}-\mathbf{h}_0\|^2}{12}\right)^{-1}.
$$

The denominator must be positive. Otherwise the geometry cannot be normalized
and is rejected.

### 7.2. Midline metric

Define the vectors between the midpoints of opposite sides:

$$
\mathbf{a}=\frac12(\mathbf{x}_2+\mathbf{x}_3-
                         \mathbf{x}_1-\mathbf{x}_4),
$$

$$
\mathbf{b}=\frac12(\mathbf{x}_3+\mathbf{x}_4-
                         \mathbf{x}_1-\mathbf{x}_2).
$$

The geometric metric is:

$$
\mathbf{M}=A\left(
\frac{\mathbf{a}\otimes\mathbf{a}}
     {(\mathbf{a}\cdot\mathbf{a})^2}
+
\frac{\mathbf{b}\otimes\mathbf{b}}
     {(\mathbf{b}\cdot\mathbf{b})^2}
\right).
$$

The modal stiffness is:

$$
\mathbf{K}_{hg}=s\,\mathbf{L}
\left(\frac{E\chi}{12}\mathbf{M}\right)\mathbf{L}^T.
$$

The final block is:

$$
\mathbf{K}_{mb}=\mathbf{K}_0+\mathbf{K}_{hg}.
$$

The same construction is used for membrane and bending with different
materials, section factors, and degree-of-freedom mappings. The bilinear modes
prevent *hourglass* mechanisms without penalizing uniform fields.

## 8. MacNeal transverse shear

The reference Reissner--Mindlin class uses the conventional strains:

$$
\gamma_{xz}=w_{,x}+\theta_y,
\qquad
\gamma_{yz}=w_{,y}-\theta_x.
$$

The active Nastran class does not interpolate these strains directly over an
irregular quadrilateral. It constructs a dual operator from the four edges,
following the `UEV/ANGLEI/EDGSHR` structure in the public NASTRAN-95 code.

### 8.1. Edge vectors

Let `e_i` be the unit vector of the edge starting at node `i`, let `l_i` be
its length, and let:

$$
a_i=-\mathbf{e}_i\cdot\mathbf{e}_{i-1}.
$$

The interpolated dual vector for edge `i` is:

$$
\mathbf{s}_i=
\frac{\mathbf{e}_i+a_i\mathbf{e}_{i-1}}{1-a_i^2}N_i
+
\frac{\mathbf{e}_i+a_{i+1}\mathbf{e}_{i+1}}{1-a_{i+1}^2}N_{i+1}.
$$

Differences in `w` and averages of rotations on the incident edges generate a
`4x12` operator `S`.

### 8.2. Sampling points

Two samples in each direction are taken on the natural centerlines:

```text
gamma_x: (0, -1/sqrt(3)), (0, +1/sqrt(3))
gamma_y: (-1/sqrt(3), 0), (+1/sqrt(3), 0)
```

### 8.3. Flexibility formulation

The stiffness is expressed as:

$$
\mathbf{K}_s=\mathbf{S}^T\mathbf{Z}^{-1}\mathbf{S},
\qquad
\mathbf{Z}=\mathbf{Z}_{r}+\mathbf{Z}_{MID3}.
$$

The code does not form `Z^-1` explicitly. It solves:

$$
\mathbf{Z}\mathbf{Y}=\mathbf{S},
\qquad
\mathbf{K}_s=\mathbf{S}^T\mathbf{Y}.
$$

This improves numerical stability and preserves symmetry.

### 8.4. Residual bending flexibility

If `Delta_x` and `Delta_y` are the midline projections and `D_b,11` and
`D_b,22` are the bending stiffnesses:

$$
f_x=\frac{\Delta x^2}{12A D_{b,11}},
\qquad
f_y=\frac{\Delta y^2}{12A D_{b,22}}.
$$

The aspect correction uses `epsilon=0.025`:

$$
\alpha=\frac{\varepsilon}
{\varepsilon+(1-\varepsilon)(\Delta x/\Delta y)^2},
$$

$$
\beta=\frac{\varepsilon}
{\varepsilon+(1-\varepsilon)(\Delta y/\Delta x)^2}.
$$

### 8.5. Torsional correction for oblique elements

Let `c` be the cosine between the two midlines:

$$
c=\frac{\mathbf{a}\cdot\mathbf{b}}
        {\|\mathbf{a}\|\,\|\mathbf{b}\|}.
$$

The factor is:

$$
m(c^2)=1+
\frac{0.4100259495463487\,c^2}
     {1-0.9156556465380681\,c^2}.
$$

Properties:

- `m(0)=1` for orthogonal midlines;
- the denominator is positive for `0 <= c^2 <= 1`;
- the function is positive, monotonic, and finite over the domain;
- it reproduces moderate-obliquity DMIG controls 81/89 and high-obliquity
  control 90;
- loaded case 60 checks an intermediate obliquity.

This law is a numerical identification specific to the project. It is not
presented as a constant published by MSC.

The residual blocks are:

$$
\mathbf{Z}_{r,x}=f_x
\begin{bmatrix}
1+\alpha m&1-\alpha m\\
1-\alpha m&1+\alpha m
\end{bmatrix},
$$

$$
\mathbf{Z}_{r,y}=f_y
\begin{bmatrix}
1+\beta m&1-\beta m\\
1-\beta m&1+\beta m
\end{bmatrix}.
$$

### 8.6. Physical MID3 flexibility

When `MID3` is present, for every sample `j`:

$$
f_{s,j}=\frac{1}{2J_j t_s G_3}.
$$

The pair mode uses:

$$
\mathbf{P}=\frac12
\begin{bmatrix}
1+m&1-m\\
1-m&1+m
\end{bmatrix}.
$$

Congruence with the roots of `f_s,j` preserves symmetry and positivity when
the Jacobian differs between the two samples. MID3 is added **in series** as
flexibility. Activating it therefore makes the element more flexible instead
of adding a parallel Reissner--Mindlin stiffness.

## 9. Drilling-rotation stabilization: K6ROT

The planar theory supplies no independent physical stiffness for `theta_z`.
The implementation extends to four nodes the penalty already validated for
CTRIA3.

At every vertex, the membrane macro-rotation is:

$$
\omega_i=\frac12(v_{,x}-u_{,y})_i.
$$

The nodal residual is:

$$
r_i=\theta_{zi}-\omega_i.
$$

The complete contribution is:

$$
\mathbf{K}_{6ROT}=\sum_{i=1}^{4}
\left(K6ROT\,10^{-6}G_1t\right)
\left(4J_i\right)\mathbf{r}_i\mathbf{r}_i^T.
$$

Derivatives and `J_i` are evaluated at the corresponding natural vertex. This
does not merely add four values to the diagonal: the residual couples
`theta_z` to the `u,v` displacements that form the macro-rotation.

This form:

- preserves rigid-body motion;
- follows the principle used for CTRIA3;
- stabilizes `theta_z` without creating an independent physical rotation;
- reproduces the drilling--membrane coupling in the validated DMIG matrices.

If `K6ROT=0`, this contribution is absent. The solver performs no general
`AUTOSPC` and does not remove inactive degrees of freedom automatically. Any
rotation without stiffness must be constrained explicitly in the model; a
free zero row or structural mechanism produces a singular-matrix error.

## 10. Local assembly and global transformation

The blocks are scattered into the `24x24` local matrix:

$$
\mathbf{K}_{loc}=
\mathbf{K}_m+\mathbf{K}_b+\mathbf{K}_s+\mathbf{K}_{6ROT}.
$$

The local indices are:

- membrane: `[T1,T2]` at every node;
- plate: `[T3,R1,R2]` at every node;
- drilling: `[R3]`, coupled to `[T1,T2]` by `K6ROT`.

The local rotation matrix has rows `[e_x, e_y, e_z]` and is repeated for the
translations and rotations of all four nodes. Transformation to basic
coordinates is:

$$
\mathbf{K}_{basic}=\mathbf{T}^T\mathbf{K}_{loc}\mathbf{T}.
$$

The result is symmetrized to remove rounding noise:

$$
\mathbf{K}\leftarrow\frac12(\mathbf{K}+\mathbf{K}^T).
$$

Global assembly additionally applies the transformation associated with the
`CD` field of every `GRID`. Nodal degrees of freedom may therefore be defined
in different systems without changing the element formulation.

## 11. Force and stress recovery

Recovery is performed at the center `(xi,eta)=(0,0)`.

### 11.1. Membrane

$$
\boldsymbol{\varepsilon}_m=\mathbf{B}_m(0,0)\mathbf{u}_m,
$$

$$
\boldsymbol{\sigma}_m=\mathbf{D}_m\boldsymbol{\varepsilon}_m,
\qquad
\mathbf{N}=t\boldsymbol{\sigma}_m.
$$

The method returns `Nx`, `Ny`, and `Nxy`.

### 11.2. Moments

$$
\boldsymbol{\kappa}=\mathbf{B}_b(0,0)\mathbf{q}_p,
$$

$$
\mathbf{M}=-I_b\mathbf{D}_b\boldsymbol{\kappa}.
$$

The negative sign adapts the internal curvature convention to the Nastran
moment convention used by the project.

### 11.3. Fiber stresses

Assuming a linear distribution through the thickness:

$$
\boldsymbol{\sigma}_{bot}=
\boldsymbol{\sigma}_m+\frac{6\mathbf{M}}{t^2},
$$

$$
\boldsymbol{\sigma}_{top}=
\boldsymbol{\sigma}_m-\frac{6\mathbf{M}}{t^2}.
$$

The result API uses one geometric convention: `top` is the fiber at
`z=+t/2`, and `bot` is the fiber at `z=-t/2`, always in element axes. The six
components are `stress_top_fiber_x`, `stress_top_fiber_y`,
`tau_top_fiber_xy`, `stress_bot_fiber_x`, `stress_bot_fiber_y`, and
`tau_bot_fiber_xy`. No `Z1/Z2` aliases or name permutations are retained;
they duplicated the same value, and `Z1/Z2` correspondence is not a portable
geometric convention across formulations and programs.

Principal stresses for each fiber are:

$$
\sigma_{1,2}=\frac{\sigma_x+\sigma_y}{2}
\pm\sqrt{\left(\frac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2},
$$

and plane-stress von Mises stress is:

$$
\sigma_{VM}=\sqrt{
\sigma_x^2-\sigma_x\sigma_y+\sigma_y^2+3\tau_{xy}^2}.
$$

### 11.4. Shear resultants Qx/Qy

The nodal forces produced by `K_s` and the plate displacement vector are
calculated first. The four generalized forces of operator `S` are then
reconstructed by least squares. The integration measure `2J_j` is removed,
and the two samples in each direction are averaged.

This recovery is consistent with the implemented stiffness but does not
exactly reproduce MSC Nastran's proprietary convention for its printed
centroidal value. Small `Qx/Qy` components may have large relative differences
without a comparable difference in the matrix, displacements, reactions, or
moments.

## 12. Validation strategy

The implementation is validated at several levels.

### 12.1. Unit tests

- `24x24` dimensions;
- Maxwell--Betti symmetry;
- six rigid-body modes;
- orthonormal local axes;
- constant membrane strain;
- membrane--plate decoupling;
- independent use of membrane and bending materials;
- exact form of the `K6ROT` penalty;
- control points and monotonicity of the oblique metric;
- rejection of collinearity, non-positive Jacobian, and warping.

### 12.2. DMIG matrix comparison

Relative Frobenius error is defined as:

$$
e_F=\frac{\|\mathbf{K}_{solver}-\mathbf{K}_{Nastran}\|_F}
          {\|\mathbf{K}_{Nastran}\|_F},
$$

where:

$$
\|\mathbf{A}\|_F=\sqrt{\sum_i\sum_j A_{ij}^2}.
$$

This measure summarizes global error, but a dominant block could conceal a
small erroneous term. Every case therefore also applies:

- termwise comparison with `rtol` and `atol`;
- comparison of the membrane block `[T1,T2]`;
- comparison of the plate/shear block `[T3,R1,R2]`;
- comparison of the drilling block `[R3]`;
- comparison of drilling--membrane coupling where applicable.

The blocks are first rotated into element axes. A much stiffer membrane block
therefore cannot conceal a shear or `K6ROT` error.

The test can export the calculated matrix to CSV through:

```python
"debug_export": True,
"output_k_file": "draft_matrix_name_python.csv",
```

### 12.3. Covered cases

- 51--55: membrane, bending, combination, three-dimensional orientation, and
  `K6ROT`;
- 65--75: residual and physical shear, `G3`, `TS/T`, absence of MID3, aspect
  correction, and skew;
- 81, 83, and 89: oblique affine quadrilaterals, active/blank MID3, and the
  full matrix;
- 90: high-obliquity DMIG control;
- 56--58, 60--61, 66, 68, and 76--91: loaded models, multiple properties,
  inclined panels, mixed structures, and prescribed displacements.

Nastran displacement, load, reaction, and result values are written directly
in the E2E test. The test does not parse the F06 file at runtime.

## 13. Closure-case results

### 13.1. Case 90

Case 90 has greater obliquity than the previous controls. It showed that the
former scalar law was too flexible. The investigation isolated the error in
the geometric continuation of the torsional mode; there was no need to alter:

- the edge operator;
- base residual flexibility;
- the MID3 contribution;
- bilinear membrane or bending;
- `K6ROT`.

With the final rational law:

- the full matrix has Frobenius error below `1e-6`;
- the plate/shear block is below `4e-6`;
- block terms satisfy `rtol=1e-5`, `atol=2e-3`.

### 13.2. Case 91

Case 91 contains two highly oblique, coplanar CQUAD4 elements with different
properties, forces, and moments. Its results are:

- displacement vector error: `2.07e-5`;
- reaction vector error: `6.30e-5`;
- membrane forces, moments, and fiber stresses: below `0.1%` per component;
- greatest `Qx/Qy` difference: `47.1%`, equal to `0.0216` in absolute value.

The last difference is classified as a centroidal-recovery limitation, not a
stiffness error.

### 13.3. Cases 57 and 58

These cases must not be used as primary controls for recalibrating stiffness.
Their `MID3` modulus is of order `1e-20`:

- case 57: `cond(K)=2.90e25`; vector errors `U=1.72e-3`, `R=4.07e-3`, and
  result error `9.28e-4`;
- case 58: `cond(K)=8.52e26`; translations are within approximately `0.3%`,
  but small rotations, secondary moments, and `Qy` are sensitive; result
  vector error is `2.94e-2`.

A large relative difference in an almost-zero quantity does not imply that
the dominant structural response is incorrect.

### 13.4. Out-of-scope cases

The six strict `xfail` cases are:

| Cases | Reason | Expected result |
|---|---|---|
| 54, 59, 62, 63, 64 | One or more warped CQUAD4 elements | The constructor rejects the geometry |
| 78 | Default `SNORM` normal smoothing | Does not match without `SNORM` support |

Equivalent faceted cases 82/88 with `SNORM,0` pass. These `xfail` cases are
not pending planar-CQUAD4 defects; they protect the scope boundary.

The complete suite observed on 19 September 2026 finished with:

```text
718 passed, 6 xfailed, 2 warnings
```

The warnings were local `PytestCacheWarning` messages caused by denied cache
write access, not numerical failures.

## 14. Explicit limitations

The implementation does not include:

- warped CQUAD4 elements or projection onto a mean plane;
- `SNORM` normal smoothing between elements;
- orthotropic materials, `MAT8`, or laminates;
- membrane--bending coupling for composites;
- `MID4`;
- `ZOFFS` or midsurface eccentricity;
- variable nodal thicknesses and `TFLAG`;
- distributed loads or face pressure;
- exact recovery of MSC's proprietary centroidal `Qx/Qy` convention;
- nonlinear analysis, large displacements, or material nonlinearity.

The formulation does not attempt to reproduce every option of the commercial
CQUAD4. Only fields covered by the reader, data model, and tests are supported.

## 15. Closure criterion

The planar CQUAD4 stiffness is considered complete because:

1. symmetry, rigid modes, axes, and uniform fields are verified;
2. active DMIG matrices are compared globally, blockwise, and term by term;
3. MID1, MID2, MID3, `12I/T^3`, `TS/T`, and `K6ROT` have been isolated;
4. rectangular, trapezoidal, skewed, and highly oblique controls exist;
5. multiple-property and mixed models reproduce displacements, reactions,
   forces, moments, and stresses;
6. unsupported geometries are rejected explicitly.

Non-blocking future improvements are:

- reproduce printed `Qx/Qy` more accurately;
- perform mesh-refinement convergence studies;
- optimize sparse assembly after measurement;
- implement `SNORM` only if future scope requires it.

Do not relax the tolerance to accept warping. A future non-coplanar
formulation should be treated as a new element with its own kinematics and
DMIG validation.

## 16. References

- R. H. MacNeal, *A Simple Quadrilateral Shell Element*, Computers &
  Structures, volume 8, pages 175--183, 1978.
- R. H. MacNeal, [*Derivation of Element Stiffness Matrices by Assumed Strain
  Distributions*](https://doi.org/10.1016/0029-5493(82)90262-X), Nuclear
  Engineering and Design 70(1), 3--12, 1982.
- NASA, [NASTRAN-95: `quad4d.f` routine](https://github.com/nasa/NASTRAN-95/blob/master/mis/quad4d.f).
- Hexagon Manufacturing Intelligence, *MSC Nastran Reference Guide*.
