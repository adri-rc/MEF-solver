# Final linear T case: three-dimensional aircraft frame

This `SOL 101` model is the three-dimensional linear integration case for
the master's thesis. A perforated web based on the historical planar model and
a transverse lower flange form a T section. The purpose is not to represent a
certifiable component or replace experimental correlation, but to exercise
CROD, CBAR, CTRIA3, and CQUAD4 in one 3D assembly with local-coordinate loads,
a nodal moment, and prescribed displacements.

The model uses coherent `mm - N - MPa` units, small deformation, and
isotropic linear elasticity. `PARAM,AUTOSPC,NO` is explicit: only declared
DOFs are constrained.

![Meshed T-frame geometry](tfm_final_airframe_frame_T_patran.png)

## Geometry and structural interpretation

The segment runs along global `X`. Its `1200 x 600 mm` web lies in the
`XY` plane at `Z=0`, with two `150 mm` radius cut-outs centred at
`(300,300,0)` and `(900,300,0)`. A `1200 x 300 mm` lower flange lies in
the `XZ` plane at `Y=0`; the two planes meet along `Y=Z=0`.

Perimeter CBAR members idealize the cap attached to a real frame web. They add
axial area and bending inertia, carry load around the cut-outs, and prevent the
sheet edge from being the only load path. Secondary CBAR members connect the
upper/lower cut-out boundaries with the cap. The idealization is centred and
has no cap/web offset.

Twelve new CBAR form six diagonal segments at `X=0, 600, 1200 mm`,
symmetrically on both web faces. Each diagonal is split at its midpoint: PID
22 lies next to the web and PID 23 next to the flange. PID 23 doubles `A`,
`I1`, `I2`, and `J` to idealize rib thickening near the frame base.

Two CROD replace the former upper triangular bracing. They lie in the web
plane and connect free load points A, GRID 376 and 377, to upper nodes 171 and
359. `CORD2R 101/102` align DOF 1 and force vectors with each tie.
Constraints `23456` at A remove only modes that an isolated CROD cannot
resist; they are neither frame supports nor rivets.

## Mesh inventory

The model contains 417 GRID and 506 elements:

| Family | Count | PID | Function |
|---|---:|---:|---|
| CQUAD4 | 256 | 10 | `t=2 mm` quadrilateral web panels |
| CTRIA3 | 128 | 11 | `t=2 mm` transitions around both cut-outs |
| CQUAD4 | 32 | 12 | `t=4 mm` T-section flange |
| CBAR | 56 | 20 | Perimeter cap, longitudinal edges, and central post |
| CBAR | 20 | 21 | Secondary cut-out reinforcement |
| CBAR | 6 | 22 | Rib halves next to the web |
| CBAR | 6 | 23 | Rib halves next to the flange |
| CROD | 2 | 30 | Load-introduction ties at A |

Totals by family are 288 CQUAD4, 128 CTRIA3, 88 CBAR, and 2 CROD. The twelve
PID 22/23 members are the added 3D ribs; the other 76 retain the starting
planar model's web reinforcement.

## Materials and properties

All elements use `MAT1 1` with `E=71700 MPa`, `nu=0.33`, and
`G=26954.887 MPa`.

| Card | Main parameters | Interpretation |
|---|---|---|
| PSHELL 10 | `t=2`, `MID1=MID2=MID3=1`, `TS/T=5/6` | Web CQUAD4 with physical membrane, bending, and shear |
| PSHELL 11 | `t=2`, `MID1=MID2=MID3=1`, `TS/T=5/6` | CTRIA3 rings with the same material |
| PSHELL 12 | `t=4`, `MID1=MID2=MID3=1`, `TS/T=5/6` | Thicker flange CQUAD4 |
| PBAR 20 | `A=235`, `I1=21389`, `I2=14854`, `J=1078` | Equivalent `40x4 + 25x3 mm` cap T section |
| PBAR 21 | `A=120`, `I1=6000`, `I2=3000`, `J=350` | Secondary reinforcement |
| PBAR 22 | `A=120`, `I1=6000`, `I2=3000`, `J=350` | Rib half next to the web |
| PBAR 23 | `A=240`, `I1=12000`, `I2=6000`, `J=700` | Doubled rib half next to the flange |
| PROD 30 | `A=50`, `J=20`, `C=2.5` | Load-introduction CROD |

This three-dimensional case uses physical shell bending and shear.
`PARAM,K6ROT,100` supplies the documented drilling stiffness and
`PARAM,SNORM,0` keeps independent normals across intersecting planes.

## Loads and boundary conditions

Four flange stations idealize discrete rivets or fixture attachment points.
The complete edge is not fixed, allowing the flange to deform between
connections.

| Flange side | GRID | Coordinates `(X,Y,Z)` mm | Condition |
|---|---|---|---|
| `Z=-150` | 405, 413 | `(300,0,-150)`, `(900,0,-150)` | Fixed `123456=0` |
| `Z=+150` | 422, 430 | `(300,0,150)`, `(900,0,150)` | `13456=0`, `T2=+0.2 mm` |

The equal displacement at both second-side rivets introduces a small
differential fixture motion. It belongs to SPC 120 and is combined with fixed
points and guides by SPCADD 130; it is not `AUTOSPC`.

The only subcase applies:

- `FORCE=6000 N` at each point A (GRID 376/377), directed outward along the
  corresponding CROD through CID 101/102;
- `MOMENT Mz=250000 Nmm` at upper midpoint B, GRID 165;
- prescribed `T2=+0.2 mm` at GRID 422/430.

![Supports, guides, and prescribed displacements](tfm_final_airframe_frame_T_patran_spc.png)

![Forces at A and moment at B](tfm_final_airframe_frame_T_patran_force_moment.png)

The following image overlays the undeformed structure in blue and an
automatically amplified deformation in white. The scale is visual only; the
analysis remains small-deformation linear statics.

![Original and amplified deformed geometry](tfm_final_airframe_frame_T_patran_deformed.png)

The shell contour shows `Stress Tensor / At Z2 / VONM`. It is a smoothed
nodal fringe for interpretation; acceptance values come from explicit F06
element stations, not graphical interpolation.

![Positive-fibre shell von Mises stress](tfm_final_airframe_frame_T_patran_von_mises_top.png)

## MSC Nastran reference

The unchanged [DAT](tfm_final_airframe_frame_T.dat) ran in MSC Nastran Student
Edition 2026.1 on 24 August 2026. The
[audited F06](tfm_final_airframe_frame_T.f06) reaches `END OF JOB`, contains
no numerical user/system warnings or `FATAL` messages, and was solved with
`AUTOSPC=NO`. Its literal `Warning` occurrence is legal header text.

MSC reports:

- stiffness-matrix condition number: `3.172218e9`;
- energy epsilon: `-3.0259522e-13`;
- largest translation magnitude: `3.985038 mm` at GRID 434;
- largest rotation magnitude: `0.0181810 rad` at GRID 415;
- largest force reaction: `2520.765 N` at GRID 413;
- largest moment reaction: `178885.39 Nmm` at GRID 413.

| File | SHA-256 |
|---|---|
| `tfm_final_airframe_frame_T.dat` | `1c3f4b0300815fc2106a28adfb5655494c459dc54424c2eb5ed6c804168b851a` |
| `tfm_final_airframe_frame_T.f06` | `edb202790970e5f9ec41a91d341112bdfcb9dfc55276cf36a466c46b699a99a9` |
| `tfm_final_airframe_frame_t.op2` (local audit) | `a3ccf11aa227ed9137dbdd2fc06e8e2222a5cf6e9f06757ae439e1950b944bac` |

## Automated validation

[`test_linear_final_t_validation_case.py`](../../test_linear_final_t_validation_case.py)
contains five checks for this physical model and four for the negligible-MID3
diagnostic. Pytest reads neither F06 nor OP2; selected MSC values are literal
test data.

The comparison includes 15 displacement GRID, all six reaction stations, both
CROD, and 15 elements from each remaining family:

- CBAR: `7001, 7025, 7026, 7030, 7056, 7066, 7072, 7076, 7101, 7102, 7105,
  7106, 7108, 7109, 7112`;
- CTRIA3: `5001, 5008, 5016, 5024, 5032, 5040, 5048, 5050, 5051, 5064, 5065,
  5096, 5111, 5112, 5128`;
- CQUAD4: `1001, 1064, 1096, 1106, 1128, 1151, 1152, 1155, 1192, 1230,
  1231, 1256, 2001, 2013, 2032`.

CBAR references include all eight F06 actions and normal stresses at four A/B
fibres. Accepted shell values include `Nx, Ny, Nxy` and positive-fibre
`sx, sy, txy, von Mises`, covering web, flange, perimeter, secondary
reinforcement, and both halves of the new ribs.

## Results and limitations

Relative error is
`||result_MEF - result_MSC||2 / ||result_MSC||2`. Full-field auditing used
binary OP2 precision; regression literals use six-digit F06 output.

| Full field | Relative L2 error | Note |
|---|---:|---|
| Translations | `8.206e-3` | Dominated by transverse response |
| Rotations | `2.547e-2` | Out-of-plane rotations are most sensitive |
| Force reactions | `1.705e-4` | Correct global closure |
| Moment reactions | `3.553e-4` | Maximum absolute difference `59.66 Nmm` |
| CROD force/stress | `1.871e-6` | Axial response nearly machine-precision identical |
| CROD torsion | `8.125e-3` | Secondary action of about `1.37 Nmm` |
| Complete CBAR force vector | `2.744e-3` | Dominant in-plane actions correlate more closely |
| CBAR normal stress | `3.930e-3` | Maximum absolute difference `0.303 MPa` |
| CQUAD4 membrane | `3.268e-5` | Web and flange |
| CQUAD4 moments / shears | `6.177e-3 / 5.883e-3` | Centroidal 3D recovery remains close |
| CQUAD4 fibre stress | `6.582e-3` | Selected core `4.790e-3` |
| CTRIA3 membrane | `7.218e-6` | Very close agreement |
| CTRIA3 fibre stress | `1.020e-2` | Selected core `6.775e-3` |
| CTRIA3 moments / shears | `6.820e-1 / 1.730e1` | Not accepted as correlation |

Maximum displacement differences are `0.000179 mm` in T1, `0.010416 mm`
in T2, and `0.050705 mm` in T3; rotation differences remain below
`1.93e-4 rad`. GRID 434 has `T2=3.985038 mm` in MSC and
`3.993304 mm` locally. At GRID 413, MSC `Fy=-2509.189 N` and
`Mx=178884.3 Nmm` become `-2509.618 N` and `178944.0 Nmm`. Both CROD
recover `6000 N` and `120 MPa` axial response in both solvers.

The nonuniform discrepancy is evidence, not a failure to hide with one broad
tolerance. Dominant CBAR `MA2`, `MB2`, `V2`, and axial values have L2
errors below `5e-4`; small orthogonal `MA1`, `MB1`, `V1`, and torsion
remain diagnostics with errors up to `0.163`. CTRIA3 membrane and fibre
stress are accepted, but centroidal moments/shears are not validated by this
deck. Their difference belongs to result recovery and does not justify
recalibrating element stiffness.

The A/B control
[case_linear_test_final_T_no_CTRIA_MID3](../case_linear_test_final_T_no_CTRIA_MID3/README.md)
changes only CTRIA3 `MID3` to an `E=1e-20` material. Moment L2 error falls
from `0.682` to `7.106e-5`, linking the discrepancy to shear
formulation/recovery. The same control has a `T3` near-mechanism and
condition number `5.660208e28`; it is diagnostic evidence, not a physically
valid replacement.

After transforming reactions from each `GRID.CD` to the basic frame, local
equilibrium closes below `1e-6 N` and `2e-3 Nmm`.

## Local execution

```powershell
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_linear_test_final_T\tfm_final_airframe_frame_T.dat
.\.venv\Scripts\python.exe -m pytest -q tests\solver_tst\test_linear_final_t_validation_case.py
```
