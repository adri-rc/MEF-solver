# Final linear case: aircraft-frame coupon

This `SOL 101` integration model was prepared for the master's thesis
defence. It represents a planar `1200 x 600 mm` aircraft frame/web coupon
with two `150 mm` radius cut-outs, six idealized lower rivets, caps,
stiffeners, and two upper load-introduction clevises.

The model is neither experimental correlation nor general validation of
Nastran. The unchanged DAT completed normally in MSC Nastran Student Edition
2026.1 on 21 August 2026 with no `FATAL` messages, closing numerical
cross-validation of this particular model. The
[audited F06](tfm_final_airframe_frame.f06) is preserved beside the deck, and
selected references are transcribed into the E2E test. The run used build
`2026.1-2a19f742`; the F06 SHA-256 is
`80aea144d3a084c886353e56d3e88e4290d625b6a894ec9203ffaad6ba721a5e`.

## Dimensions and contents

The model uses coherent `mm - N - MPa` units and small deformation:

| Part | Idealization | Property |
|---|---|---|
| Web | 256 CQUAD4 | PSHELL 10, aluminium, `t=2.0 mm` |
| Two cut-out boundaries | 128 CTRIA3 | PSHELL 11, `t=2.0 mm`, nearly zero MID3 |
| Perimeter cap and central post | 56 CBAR | PBAR 20, equivalent `40x4 + 25x3 mm` T section |
| Upper/lower cut-out stiffeners | 20 CBAR | PBAR 21, secondary equivalent section |
| Pull clevises | 4 CROD | PROD 30, `A=50 mm2`, `J=20 mm4` |

There are 377 GRID and 464 elements. The web is coplanar. Its nodal
`T3/R1/R2` DOFs are guided by the fixture so that the test exercises
in-plane membrane response. This avoids assigning physical accuracy to
CTRIA3 transverse-shear flexibility, which is a declared solver limitation.
`R3` remains active for shell drilling and in-plane CBAR bending.

PBAR 20 represents a centred T section with a `40 x 4 mm` flange and
`25 x 3 mm` web: `A=235 mm2`, `I1=21389 mm4`, `I2=14854 mm4`, and
`J=1078 mm4`. No offsets are applied because cap/web eccentricity would
require an additional modelling decision.

## Loads and boundary conditions

- SPC 100 provides planar fixture guidance and constrains inactive clevis DOFs.
- SPC 110 fixes GRID `181, 184, 187, 366, 369, 372` on the lower row; with
  SPC 100 this is equivalent to `123456`.
- SPC 120 prescribes `T2=+0.8 mm` at the right clevis.
- SPCADD 130 is selected by the subcase.
- LOAD 201 applies `Fy=+12000 N` at the left clevis and `Fx=+2500 N` at
  the upper central-post joint.
- LOAD 202 applies `Mz=+250000 Nmm` at the same joint.
- LOAD 200 combines 201 and 202.

Force, moment, and prescribed displacement share one subcase. This is useful
as a final integration demonstration. A more diagnostic study would duplicate
the deck into three independent load cases without changing mesh or
properties.

## Generation and local execution

From the repository root:

```powershell
.\.venv\Scripts\python.exe tests\solver_tst\validation_cases\case_linear_test_final\section_frame_model.py
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_linear_test_final\tfm_final_airframe_frame.dat
```

The report is written to `output/tfm_final_airframe_frame.rpt`. The
deterministic generator is the editable source of the mesh.

## Correlation with MSC Nastran 2026.1

Representative points are:

| Quantity | MEF-solver | MSC Nastran F06 |
|---|---:|---:|
| Left clevis, GRID 376, `T2` | `0.744011933 mm` | `0.7440119 mm` |
| Right clevis, GRID 377, `T1` | `0.055691117 mm` | `0.05569107 mm` |
| Actuator reaction at GRID 377, `F2` | `12993.822458 N` | `12993.82 N` |
| Maximum CROD axial stress, EID 8003/8004 | `152.976475 MPa` | `152.9765 MPa` |
| Maximum CBAR normal stress, EID 7056 | `131.623096 MPa` | `131.6232 MPa` |
| Maximum shell von Mises, EID 5030 | `43.074543 MPa` | `43.07457 MPa` |

Full OP2 arrays were compared before selecting representative points:

| Full field | Relative L2 error |
|---|---:|
| Displacements | `1.052e-6` |
| SPC reactions | `1.450e-5` |
| CROD forces / stresses | `4.945e-8 / 5.701e-8` |
| CBAR forces / stresses | `5.138e-6 / 5.361e-6` |
| CTRIA3 resultants / stresses | `6.798e-6 / 5.701e-6` |
| CQUAD4 resultants / stresses | `1.073e-5 / 8.485e-6` |

The automated test closes global force equilibrium below `1e-6 N` and
moment equilibrium below `1e-3 Nmm`. Values for several GRID/EID from every
family are embedded explicitly in `test_linear_final_validation_case.py`;
pytest does not read F06 or OP2 at runtime. Printed F06 tables are rounded,
whereas the vector errors above used OP2 precision.

## Patran import and execution

The workflow follows the **Read Input File** process in the Patran/MSC Nastran
interface guide:

1. Create a database with `File > New`, open
   `Preferences > Analysis`, and select `MSC Nastran`. Patran assigns no
   units; retain `mm-N-MPa`.
2. In `Analysis`, choose `Action = Read Input File`, select
   `tfm_final_airframe_frame.dat`, and keep nodes, elements, properties,
   materials, loads, constraints, and load cases enabled.
3. Review `Summary Data`, `Reject Cards`, and `View Log File`. GRID,
   CQUAD4, CTRIA3, CBAR, CROD, PSHELL, PBAR, PROD, MAT1, FORCE, MOMENT,
   SPC/SPC1/SPCADD, LOAD, and the subcase should all be present.
4. To solve from the imported database, select linear static analysis, the
   complete model, and a full run. Preserve Patran's regenerated BDF so its
   intent can be compared with the source deck.
5. The most auditable route is to run the original DAT directly in MSC
   Nastran and then use `Analysis > Access Results > Read Output2`.
   `PARAM,POST,-1` and `PLOT` requests are already present.

Reference: [Patran Interface to MSC Nastran Preference Guide, Volume 1,
Structural Analysis](https://documentation-be.hexagon.com/bundle/Patran_2024.1_Interface_to_MSC_NASTRAN_Preference_Guide_Volume_1_Structural_Analysis/raw/resource/enus/Patran_2024.1_Interface_to_MSC_NASTRAN_Preference_Guide_Volume_1_Structural_Analysis.pdf?save_local=true),
chapters 3--5.

For presentation, group entities by PID/family, display PSHELL thickness,
render beam sections, and overlay FORCE, MOMENT, and SPC symbols. A second
view with amplified deformation and top/bottom von Mises contours makes the
directly comparable quantities clear.
