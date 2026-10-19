# Final nonlinear validation case: three-span spatial arch bridge

## Purpose and scope

This case validates the mixed nonlinear `CBAR/CROD` assembly on a comparatively
large structure. The model represents an idealized fragment of a through-arch
railway bridge with three spans, a spatial deck, hangers, bracing, and explicit
abutment and pier bases. The supplied photographs were used only as conceptual
references for the sequence of spans, hangers, and foundations; they provide no
dimensions, properties, or loads for the analysis.

The `_updated` suffix is retained solely for traceability with the delivered
thesis and its archived validation discussion. This directory contains the
authoritative public validation case; the suffix does not indicate a dependency
on another repository case.

The model uses `mm-N-MPa` units and remains within the subset actually
implemented by the project:

- `SOL 400`, `ANALYSIS=NLSTAT`, and geometric nonlinearity;
- isotropic elastic `MAT1` materials;
- `PBAR/PBARN1/CBAR` and `PROD/CROD`;
- conservative nodal `FORCE` and `MOMENT` loads;
- explicit `SPC/SPC1/SPCADD` constraints;
- full Newton `PFNT` with fixed increments.

It does not use solids, shells, contact, plasticity, automatic gravity,
distributed loads, follower loads, or `AUTOSPC`.

## Geometry

The longitudinal axis is `X`, the bridge width lies along `Y`, and `Z` is
vertical. The deck is `30 000 mm` long and `3 600 mm` wide. It has 31 stations
spaced at `1 000 mm` and four longitudinal lines at
`Y=-1800,-600,+600,+1800 mm`.

The supports are located at `X=0, 8000, 20000, 30000 mm`, giving span lengths
of `8000`, `12000`, and `10000 mm`. Each outer plane contains three consecutive
parabolic segments with rises of `1900`, `3200`, and `2400 mm`, respectively.
The arch springing nodes coincide with the outer deck girders.

The four support lines include two columns beneath the deck edges, four inclined
legs, and four ground nodes. The ground elevations are staggered so that the
abutments and the two piers have different heights:

| Support X [mm] | Column-bottom Z [mm] | Ground elevation Z [mm] |
|---:|---:|---:|
| 0 | -900 | -1800 |
| 8000 | -1800 | -3200 |
| 20000 | -2200 | -3600 |
| 30000 | -1100 | -2100 |

## Entity inventory

| Entity | Count |
|---|---:|
| GRID | 202 |
| CBAR | 144 |
| CROD | 428 |
| Total elements | 572 |
| MAT1 | 2 |
| PBAR + PBARN1 | 4 + 4 |
| PROD | 7 |
| FORCE | 44 |
| MOMENT | 3 |
| LOAD | 1 |

GRID numbering is grouped as follows:

| Group | GRID | Count |
|---|---|---:|
| Four deck lines | `1001-1031`, `1101-1131`, `1201-1231`, `1301-1331` | 124 |
| Outer arch at `Y=-1800` | block `2002-2030`, except the reused springing nodes | 27 |
| Outer arch at `Y=+1800` | block `2102-2130`, except the reused springing nodes | 27 |
| Lower pier caps | `3001/3002`, `3011/3012`, `3021/3022`, `3031/3032` | 8 |
| Ground bases | `4001-4004`, `4011-4014`, `4021-4024`, `4031-4034` | 16 |

## Elements and associated properties

### CBAR

`CBAR` elements are reserved for members whose bending and torsional stiffness
is relevant. Every `PBAR` has an associated `PBARN1` with `SECT=S`, `BEH2=BAR`,
and `INT2=LC`.

| PID | Main EIDs | Count | Representation | A [mm²] | I1 [mm⁴] | I2 [mm⁴] | J [mm⁴] |
|---:|---|---:|---|---:|---:|---:|---:|
| 101 | `10001-10030`, `10101-10130` | 60 | two outer girders/ties | 18000 | 1.20e9 | 4.50e8 | 9.00e7 |
| 103 | `11001-11030`, `11101-11130` | 60 | two three-span arch ribs | 12000 | 7.00e8 | 3.20e8 | 5.00e7 |
| 105 | `12001-12008` | 8 | pier and abutment columns | 24000 | 1.60e9 | 1.00e9 | 1.50e8 |
| 106 | `12101-12116` | 16 | inclined foundation legs | 16000 | 8.00e8 | 5.00e8 | 9.00e7 |

### CROD

`CROD` elements form the pin-jointed part of the truss. They are bilateral
axial elements: a hanger can mathematically enter compression, and neither
slackening nor `tension-only` cable behavior is modelled.

| PID | Main EIDs | Count | Representation | A [mm²] | MID |
|---:|---|---:|---|---:|---:|
| 201 | `20401-20454` | 54 | outer vertical hangers | 800 | 2 |
| 202 | `20001-20030`, `20101-20130` | 60 | axial inner girders | 2500 | 1 |
| 203 | `20201-20293` | 93 | transverse deck links | 1400 | 1 |
| 204 | `20701-20730`, `21001-21060` | 90 | outer and inner plan bracing | 900 | 2 |
| 205 | `20501-20527`, `20601-20630` | 57 | ties and diagonals between arch planes | 650 | 2 |
| 206 | `20801-20812` | 12 | axial ground-base beams | 3000 | 1 |
| 207 | `20901-20962` | 62 | network diagonal hangers and cap links | 700 | 2 |

The two materials are:

| MID | Use | E [MPa] | G [MPa] | ν |
|---:|---|---:|---:|---:|
| 1 | primary members and deck/base truss | 210000 | 80769.2308 | 0.30 |
| 2 | hangers and bracing | 195000 | 75000 | 0.30 |

## Structural idealization

The outer girders act as arch-system ties and provide flexural continuity to
the deck. The upper ribs transfer load through two hanger families: vertical
hangers in each outer plane and diagonal hangers connecting the inner lines to
the arch. The second family provides a direct vertical load path and prevents
the inner lines from forming an out-of-plane mechanism.

Plan bracing connects the inner longitudinal girders to the outer girders;
without it, each inner `CROD` line could translate rigidly in `X`. Transverse
ties and wind diagonals alternate between the two arch planes. Piers and
abutments are represented by `CBAR` columns and flared legs anchored at four
ground GRID points. This is a 1D idealization of foundation load
transfer, not a concrete or soil-structure interaction model.

## Boundary conditions

`PARAM,AUTOSPC,NO` is written explicitly in the DAT. Every degree of freedom
without physical stiffness is handled explicitly:

- the 16 ground GRID points restrain `T1`, `T2`, `R1`, `R2`, and `R3`;
- the four GRID points at the left support also have `T3=0`;
- the other three bases receive uniform imposed settlements of `-2`, `-6`,
  and `-10 mm` in `T3`, ramped with the load factor;
- the 62 GRID points on the two inner girders, connected only by `CROD`, have
  their `R1-R3` rotations restrained because those rotations do not belong to
  the axial truss kinematics.

These constraints do not silently stabilize a mechanism: they are declared
parts of the pin-jointed idealization and ground anchorage. The first mesh
topology, before the network hangers and inner bracing were added, correctly
produced a singular matrix with `AUTOSPC=NO`; the preserved DAT contains the
corrected stable topology.

## Loads

`LOAD 500` combines three families:

| SID | Actions | Count | Factor in LOAD 500 |
|---:|---|---:|---:|
| 501 | asymmetric vertical traffic on the two inner lines | 24 FORCE | 1.00 |
| 502 | lateral wind at eight stations on both arches | 16 FORCE | 0.75 |
| 503 | longitudinal braking and nodal couples | 4 FORCE + 3 MOMENT | 1.00 |

The applied resultants printed by MSC in the basic system are `Fx=350000 N`,
`Fy=363525 N`, and `Fz=-2592000 N`. The SID 503 nodal moments are two couples
about `X` (`+1.10e8` and `-8.00e7 Nmm`) and one about `Y` (`-1.40e8 Nmm`). All
directions are fixed in the basic system; `PARAM,FOLLOWK,NO` excludes follower
load stiffness.

## Nonlinear control and execution

The step uses `NLSTEP 500`, twenty fixed increments, up to 80 iterations and
24 permitted bisections per increment, residual criterion `PV=3e-7`, and full
Newton `PFNT`:

```text
NLSTEP,500,1.0,SEVERELY
,GENERAL,80,1,24
,FIXED,20
,MECH,PV,3.E-7,3.E-7,3.E-7,PFNT
```

To regenerate and solve locally:

```powershell
.\.venv\Scripts\python.exe tests\solver_tst\validation_cases\case_nonlinear_test_final_updated\bridge_model_updated.py
.\.venv\Scripts\python.exe main.py --file tests\solver_tst\validation_cases\case_nonlinear_test_final_updated\tfm_final_nonlinear_bridge_updated.dat
```

For MSC Nastran 2026.1 Student Edition, from the case directory:

```powershell
& 'C:\Program Files\MSC.Software\NaPa_SE\20261\Nastran\bin\nastran.exe' `
  'tfm_final_nonlinear_bridge_updated.dat' 'scr=yes' 'old=no'
```

MSC reached `LOAD STEP=1.0`, completed all 20 increments without bisections,
printed `JOB CONVERGES FOR THE CURRENT STEP`, and terminated with `END OF JOB`,
without `USER FATAL` or `SYSTEM FATAL` messages. The final-state residual ratio
is `6.369e-8`, below the requested criterion. The quantitative comparison and
the values transcribed into the end-to-end test are detailed in
[`RESULTS_ANALYSIS.md`](RESULTS_ANALYSIS.md).

## Patran views

The following images were captured in Patran 2026.1 after importing this
case's DAT and OP2. Values next to the glyphs were hidden to avoid overlap;
the exact magnitudes and resultants remain documented above and in the F06.

Spatial mesh and staggered bases:

![Bridge geometry in Patran](tfm_final_nonlinear_bridge_updated_patran_geometry_isometric.png)

Physical supports and imposed settlements. For readability, the image shows
the `SPC1 10/11` anchors and `SPC 20` settlements but omits `SPC1 30`, which
only blocks the kinematically inactive rotations of the 62 `CROD`-only GRID
points:

![Supports and settlements in Patran](tfm_final_nonlinear_bridge_updated_patran_supports_isometric.png)

Vertical traffic loads in an oblique view:

![Traffic loads in Patran](tfm_final_nonlinear_bridge_updated_patran_loads_traffic_oblique.png)

Lateral wind loads, shown in perspective to make their direction visible:

![Wind loads in Patran](tfm_final_nonlinear_bridge_updated_patran_loads_wind.png)

The longitudinal braking forces and nodal moments are documented numerically
in the load section above. No separate verified Patran
image is retained for this load family.

Combined view of all three load families, supports, and settlements:

![Combined loads and boundary conditions in Patran](tfm_final_nonlinear_bridge_updated_patran_loads_supports.png)

Final deformed configuration, with the original geometry in blue and the
deformed geometry in white. Patran magnifies displacements graphically to make
the response visible; the plotted deformation is not to scale:

![Final deformed configuration in Patran](tfm_final_nonlinear_bridge_updated_patran_deformed_final.png)

Translational displacement magnitude in the final converged state:

![Final displacement magnitude in Patran](tfm_final_nonlinear_bridge_updated_patran_displacement_magnitude.png)

No stress contour is included. The results database that Patran built from the
consolidated OP2 exposes only applied loads, displacements, and constraint
forces for the 20 states. That OP2 contains no `force/stress/strain` tables for
`CROD/CBAR`, although element results are printed in the F06 and are used for
the correlation. Producing a contour from that source would have attributed to
Patran a field that is not actually available.

## Evidence files

- `tfm_final_nonlinear_bridge_updated.dat`: exact model solved by both codes;
- `tfm_final_nonlinear_bridge_updated.f06`: MSC Nastran listing for all 20
  incremental states;
- `tfm_final_nonlinear_bridge_updated.pch`: structured displacements, loads,
  and reactions;
- `bridge_model_updated.py`: deterministic DAT generator;
- `RESULTS_ANALYSIS.md`: correlation metrics and scope observations;
- the ten `*_patran_*.png` files: verified Patran views of geometry, boundary
  conditions, loads, and the converged response.

The three numerical reference artifacts are identified by SHA-256:

| File | SHA-256 |
|---|---|
| `tfm_final_nonlinear_bridge_updated.dat` | `7c5f374b77a0abf83aab1e4c31320425faac5bd16335d32c9c0e5f9beac05439` |
| `tfm_final_nonlinear_bridge_updated.f06` | `4180efc0678eeec374d8c14d9f7f712d57985c300564e1d1da25010a660a91be` |
| `tfm_final_nonlinear_bridge_updated.pch` | `232dce9c28494fbee56cadfbbed0d03afbf8a80172d00051b82109407311c75e` |

The OP2 is retained locally to build the Patran views but is not treated as a
versioned source of truth. Test references are literals extracted from the
final F06/PCH block, and the test does not read external result files while it
runs.

## Interpretation limits

This case validates the academic solver within its scope; it does not design a
real bridge. The steel remains elastic even if a recovered stress might exceed
an allowable value. Buckling, fatigue, yielding, connections, cable slackening,
unilateral supports, and soil behavior are not assessed. Bases, cross-members,
and hangers are lines without volume. A Patran plot for this mesh must present
axial quantities or `CBAR/CROD` fibre stresses separately; there is no
continuous von Mises surface field comparable to a shell model.
