# Linear T-frame diagnostic with negligible CTRIA3 MID3 stiffness

This model is the A/B control for
[case_linear_test_final_T](../case_linear_test_final_T/README.md). Geometry,
connectivity, structural properties, loads, and boundary conditions are
unchanged. The **only** modification replaces transverse-shear `MID3` in the
128 CTRIA3 with an auxiliary `E=1e-20 MPa` material. Its purpose is to
isolate the discrepancy in centroidal CTRIA3 moment/shear recovery; it is not
a physically usable frame variant.

The name `no_CTRIA_MID3` abbreviates the numerical experiment. PSHELL
`MID3` is not blank: it points to `MAT1 99`, whose stiffness is negligible.
This preserves the remainder of the model and lets both codes recover shell
terms under the same nearly zero physical-shear assumption.

## Exact change from the physical model

Original T-frame triangle property:

```text
PSHELL,11,1,2.0,1,,1,0.8333333333
```

Diagnostic property:

```text
PSHELL,11,1,2.0,1,,99,0.8333333333
MAT1,99,1.0E-20,,0.33
```

Thus `MID1=1` and `MID2=1` retain aluminium membrane and bending; only
`MID3` changes from `MAT1 1` to `MAT1 99`. PID 10/12 CQUAD4 continue to
use `MID3=1`. No GRID, EID, PID, thickness, orientation, coordinate system,
load, SPC, or solution parameter changes.

## Preserved model

The coherent `mm - N - MPa`, small-deformation T frame still contains:

| Family | Count | PID | Function |
|---|---:|---:|---|
| CQUAD4 | 256 | 10 | `t=2 mm` perforated web panels |
| CTRIA3 | 128 | 11 | `t=2 mm` cut-out transitions |
| CQUAD4 | 32 | 12 | `t=4 mm` T flange |
| CBAR | 56 | 20 | Perimeter cap and central post |
| CBAR | 20 | 21 | Secondary reinforcement |
| CBAR | 6 | 22 | Rib halves next to the web |
| CBAR | 6 | 23 | Doubled rib halves next to the flange |
| CROD | 2 | 30 | Load-introduction ties |

Totals remain 417 GRID and 506 elements: 288 CQUAD4, 128 CTRIA3, 88 CBAR,
and 2 CROD. Both `6000 N` forces at A, the `250000 Nmm` moment at B, and
the `T2=+0.2 mm` differential motion remain. Two flange points are fixed,
two are guided prescribed-displacement rivets, and CROD free-end guides are
identical to the original model.

`PARAM,AUTOSPC,NO` remains explicit. Neither MSC Nastran nor the local solver
should hide the deliberately introduced transverse stiffness loss through
automatic stabilization.

## MSC Nastran reference

[tfm_final_airframe_frame_T_no_CTRIA_MID3.dat](tfm_final_airframe_frame_T_no_CTRIA_MID3.dat)
was run with MSC Nastran Student Edition 2026.1 on 24 August 2026. The
[F06](tfm_final_airframe_frame_t_no_ctria_mid3.f06) reaches `END OF JOB`,
contains no `FATAL`, and exits with `NSEXIT: EXIT(0)`.

MSC reports:

- stiffness-matrix condition number: `5.660208e28`;
- energy epsilon: `1.2989111e-13`;
- largest translation magnitude: `3.982957 mm` at GRID 434;
- largest `|T3|`: `0.466976792 mm` at GRID 364.

The condition number is almost twenty orders of magnitude above the physical
model's `3.172218e9`. Normal MSC termination does not make this a
well-conditioned structure.

| File | SHA-256 |
|---|---|
| `tfm_final_airframe_frame_T_no_CTRIA_MID3.dat` | `9088444b921ee252f53310b76aa2e2bc9f9bfeab951c4dd63e216dfc803a789e` |
| `tfm_final_airframe_frame_t_no_ctria_mid3.f06` | `5d82a82c9a2bf9951e8b4f1211674f9e503726e8dab244eb9d06eadf600baec0` |
| `tfm_final_airframe_frame_t_no_ctria_mid3.op2` | `f52d9e98fe108fc60882c82d1efdd3501dbdeeda042a27c1d439d2476e8f4236` |

## What the comparison establishes

The table uses full OP2 fields and relative error
`||result_MEF - result_MSC||2 / ||result_MSC||2`. CTRIA3 `core 7` combines
`Nx, Ny, Nxy` with the four positive-fibre stress quantities used by
validation. Selected printed results are embedded in
[`test_linear_final_t_validation_case.py`](../../test_linear_final_t_validation_case.py);
pytest reads neither F06 nor OP2.

| Quantity | Physical MID3 T frame | `E(MID3)=1e-20` | Interpretation |
|---|---:|---:|---|
| Complete translation vector | `8.206e-3` | `1.745e9` | Second value is dominated by the local `T3` near-mechanism |
| Rotations | `2.546e-2` | `5.079e-6` | Very close A/B control correlation |
| Force reactions | `1.705e-4` | `2.298e-7` | Close global closure |
| Moment reactions | `3.553e-4` | `4.897e-8` | Close global closure |
| CTRIA3 membrane `N` | `7.218e-6` | `7.218e-6` | MID3 does not alter membrane |
| CTRIA3 moments `M` | `6.820e-1` | `7.106e-5` | About `9.6e3` improvement |
| Relative CTRIA3 shears `Q` | `1.730e1` | `4.612e12` | Second ratio is meaningless because MSC is nearly zero |
| CTRIA3 `core 7` | `5.820e-3` | `6.880e-6` | Membrane and positive fibre correlate |

Controlled displacement components also correlate: L2 errors are
`1.795e-7` for `T1`, `3.233e-8` for `T2`, and `5.079e-6` for
rotations. Maximum absolute differences are `2.806e-7 mm` in T1 (GRID 164),
`3.831e-7 mm` in T2 (GRID 364), and below `1.40e-7 rad` in rotations. At
six reaction stations, maxima are `4.692e-4 N` in F3 (GRID 430) and
`1.213e-2 Nmm` in M1 (GRID 422).

Both CROD have aggregate error `1.248e-10`; maxima are `7.618e-7 Nmm` in
torque and `9.522e-8 MPa` in torsional stress. For 15 selected CBAR, force
and stress L2 errors are `2.237e-6 / 2.140e-6`; maxima are `0.20584 Nmm`
in EID 7056 `MA2` and `3.608e-4 MPa` in its `MINA`. The 15 CQUAD4 have
`2.160e-5` complete-resultant and `5.269e-6` positive-fibre errors; maxima
are `1.579e-3 N/mm` in `Qy` and `1.418e-4 MPa` in von Mises, both EID
1151.

For selected CTRIA3, maximum membrane difference is `2.247e-4 N/mm` in EID
5065 `Ny`. Component moment L2 errors are
`4.11e-5 / 4.64e-5 / 7.36e-5` for `Mx/My/Mxy`, with maximum
`1.455e-5 Nmm/mm` in EID 5112 `My`. Positive-fibre error is `3.675e-6`,
with maximum `1.122e-4 MPa` in EID 5111 `smin`.

For `Q`, the complete local-field norm is about `4.002e-11 N/mm` and MSC's
is `8.677e-24 N/mm`; maximum absolute difference is only
`1.867e-11 N/mm` in EID 5084 `Qy`. Dividing by a nearly zero reference
creates the apparent `4.612e12` ratio. Normalized by physical-case `Q`
scale, the difference is `6.817e-10`. The correct conclusion is that both
recover engineering-zero shear, not that a twelve-order relative discrepancy
is physically significant.

## Transverse near-mechanism and evidence boundary

The local solver obtains `|T3|=7.3586443229e10 mm` at GRID 203, while MSC
keeps `T3=-0.342126399 mm` there and a finite transverse field. `T3`
relative L2 error is `2.163e10`. This must not be hidden behind a broad
tolerance. With physical MID3 effectively removed, the local formulation
leaves triangle transverse motion nearly uncontrolled. MSC's finite field is
consistent with different residual flexibility or stabilization, but DAT/F06
cannot uniquely identify that internal treatment. The condition number
confirms a numerical limit, not a robust structure.

This case therefore gives causal evidence about **element-result recovery**:
CTRIA3 moment discrepancy drops from about 68% to `7.1e-5`, while membrane,
fibre, and other selected results remain correlated. It does not validate
`T3`, replace the physical-MID3 model, or expand solver scope. Linear
validation is closed only within these explicit limits.

## Local execution

```powershell
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_linear_test_final_T_no_CTRIA_MID3\tfm_final_airframe_frame_T_no_CTRIA_MID3.dat
.\.venv\Scripts\python.exe -m pytest -q tests\solver_tst\test_linear_final_t_validation_case.py
```

Geometry/support/load/deformation/stress images for the physical model remain
in the original case. They are not presented as validation images for this
variant because the local `T3` deformation is dominated by the deliberately
introduced near-mechanism.
