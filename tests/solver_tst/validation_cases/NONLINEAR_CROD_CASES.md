# MSC Nastran SOL 400 campaign for CROD

The initial campaign contains 20 cases, from `case_92` through `case_111`.
After their results were audited, eleven focused controls, `case_112` through
`case_122`, were added for objectivity, a central-difference tangent, torsion,
behavior near the limit point, and a linear limit sensitive to geometry
changes. Cases 123--127 form a structural-scale campaign with 10--30 rods,
mixed properties, multiple loads, and prescribed displacements; their MSC
references have been audited and incorporated into the end-to-end tests.

## Purpose of each case

| Case | Model | Main check |
|---:|---|---|
| 92 | Axial CROD, `LGDISP=-1` | SOL 400 reference without geometric nonlinearity |
| 93 | Same CROD, small load, `LGDISP=1` | Small-strain limit against 92 |
| 94 | CROD, large tension, 20 increments | Axial curve and final state with `FIXED` |
| 95 | Replica of 94 with `ADAPT` | Stepping independence, iterations, and residuals |
| 96 | CROD, finite compression | Signs, shortening, force, and stress |
| 97 | Inclined 3D CROD with `GRID.CD`/`FORCE.CID` | Transformations before iteration |
| 98 | Two identical CRODs in series | Assembly and internal-force continuity |
| 99 | Series with two `MAT1` and two `PROD` entries | Element-specific states and compatibility |
| 100 | Two coincident CRODs in parallel | Nonlinear load sharing |
| 101 | Axial CROD with `LOAD` | Combination reduction before stepping |
| 102 | Symmetric two-CROD V | First large rotation and symmetry |
| 103 | Asymmetric V with combined loading | Rotation, different sections, and 2D equilibrium |
| 104 | Shallow V, `ADAPT`, advanced diagnostics | Convergence near stiffness loss |
| 105 | Closed triangular truss | Roller support, compatibility, and horizontal reaction |
| 106 | X-braced square | Redundant 2D assembly and sharing among six rods |
| 107 | Tetrahedral tripod | Equilibrium and 3D geometric rotation |
| 108 | Square-base pyramid | Four CRODs, two sections, and 3D loading |
| 109 | Shallow V with 5 increments | Step sensitivity, paired with 110 |
| 110 | Shallow V with 40 increments | Convergence under increment refinement |
| 111 | Spatial mast, 9 CRODs, 2 materials | Final end-to-end case, `LOAD`, energy, and adaptive stepping |
| 112 | Symmetric V stable at moderate load | Common reference for cases 113--115 and 118 |
| 113 | Copy of 112 rigidly rotated in 3D | Objectivity and `GRID.CD`/`FORCE.CID` transformations |
| 114 | Case 112 with a `+1000` horizontal perturbation | Positive branch of a central-difference tangent probe |
| 115 | Case 112 with a `-1000` horizontal perturbation | Negative branch, symmetry, and deformed tangent |
| 116 | One CROD in pure torsion with `J>0` | Rotational DOF, internal torque, and torsional stress |
| 117 | One CROD under combined tension and torsion | Axial-torsional coupling/independence at finite length |
| 118 | Copy of 112 with `LGDISP=2` | Equivalence with non-follower forces |
| 119 | Symmetric V at `3.8e6`, 60 increments | Refined stable state immediately before the jump in 102 |
| 120 | V from 102 with aggressive `ADAPT` | Step changes, cutbacks, and jump; may fail informatively |
| 121 | Symmetric V with a small load and `LGDISP=-1` | Geometrically sensitive linear limit |
| 122 | Copy of 121 with `LGDISP=1` | Convergence to the linear limit in a truss that rotates |
| 123 | 2D arched roof, 10 CRODs | Redundant assembly, three PROD entries, and seven nodal forces |
| 124 | 2D truss bridge, 20 CRODs | Five PROD entries, distributed nodal loads, and a support with prescribed `T1=-4` |
| 125 | Space dome, 30 CRODs | Large 3D displacement, four supports, and two materials |
| 126 | Space canopy, 24 CRODs | Combined 3D forces and a prescribed support displacement |
| 127 | 3D braced tower, 30 CRODs | Six PROD entries, two materials, and `ADAPT` stepping |

## Complete audit

Results observed with **MSC Nastran Student Edition 2026.1**, using the
`NLPERF` engine, on 20 August 2026:

- the 30 converged jobs 92--119/121--122 finished with `exit number 0`,
  reached load factor `1.0`, and contain no fatal errors or non-finite values;
- 92/93 agree at `u=1.0` and reaction `-2.1e5`; with `LGDISP=-1`, case 92
  does not print the `RODNL` table, so its check uses displacement/reaction;
- 94/95 agree exactly in the final state (`u=20`, `sigma=42000`, reaction
  `-4.2e6`); the rerun of 95 uses nine nonuniform adaptive increments;
- 104 and 110 agree over the complete printed curve. Case 109, with five
  steps, reaches the same displacement, stress, and reaction, but its
  accumulated energy/work differs by `0.65 %`, which is precisely the intended
  sensitivity;
- 102 jumps in the final increment from `uy=-42.830` to `uy=-223.31`. This is
  not an erroneous termination: it converged in 14 iterations with a residual
  ratio of `4.66e-16`, but belongs to the post-*snap-through* branch and should
  not be used as the solver's first debugging case;
- 105 has another converged jump in increment 17 at load factor `0.48571`:
  `uy` changes from `-33.975` to `-218.50` in 14 iterations. This is also a
  *snap-through* check, not a job failure;
- the 95/104/109/110/111 reruns no longer contain warning 602;
  95/104/111 use 9/14/15 actual adaptive increments;
- 112/113/118 agree when the rotated model is expressed in its output axes;
  114/115 give a central tangent compliance of `6.127186e-6`;
- 116/117 recover `R1=0.9904762`, torque `40000`, and torsional stress `2000`;
- 119 converges before the limit with `uy=-40.32496`; 120 terminates
  abnormally at `lambda=0.98387` after reaching the minimum step, a result
  retained as a negative reference;
- 121/122 give `uy=-0.269374/-0.2699212`, a geometric difference of `0.2031 %`;
- 123--127 finished with `exit number 0`, reached load factor `1.0`, and
  reproduced `25/30/30/40/14` increments. The exhaustive comparison of the
  printed tables obtained maximum relative errors of `3.77e-7` in
  displacements, `3.76e-7` in reactions, and `3.72e-7` in stress/strain;
- the complete adaptive sequence for 127 agrees from `0.025` through `1.0`,
  including the growing steps and the final five increments of `0.1`;
- 124 and 126 exactly confirm the final prescribed displacements `T1=-4.0`
  and `T1=+5.0`, together with their conjugate reactions;
- 125 exhibits *snap-through* at GRID 5 between `lambda=0.533333` and
  `0.566667`, where `T3` jumps from `-13.5215` to `-69.40416`. The final
  mechanical state agrees, but MSC's `TotStrEn=4.7083e7` accumulates work
  through the jump and is not compared with the state elastic energy
  `1.3756e7`.

No job from 92--127 remains to be rerun.

## Execution and files to preserve

Run each `.dat` from its own directory so that the outputs remain next to the
model. Preserve at least the `.dat` and `.f06` in the repository. The F06
provides the printed results, `NLSTEP` summary, iterations, and MSC version.
Also preserve the `.pch` only when a state extraction will be automated.
`.op2`, `.f04`, `.itr`, `.log`, and `.sts` are working artifacts; for a failed
case, retain `.log` and `.sts` temporarily until it is diagnosed, but they are
not permanent references by default.

Do not rename the job root: the tests associate references through the
`ejNNN_...` prefix.

For each F06, record whether MSC used the classic `CNLS` engine or the newer
`NLPERF` engine, together with the exact version. If a case does not converge,
do not change tolerances or loading yet: also return the failure's `.f06`,
`.sts`, and `.log` so that the deck or physical phenomenon can be assessed
before choosing another strategy.

## Priority intermediate data

- Cases 92/93: final displacement/force and their difference; the `RODNL`
  stress table appears only in 93.
- Cases 94/95: every load--displacement and axial-force--increment pair.
- Cases 104/109/110: the F06 `%` table, worst residual, bisections, and the
  load--displacement curve for node 3.
- Case 111: `GPFORCE`, energy, forces in all nine CRODs, and displacements of
  nodes 5/6 at every output state.
- Cases 112/113/118: equality of results when expressed in the local axes of
  the base model.
- Cases 114/115: calculate `(u_plus-u_minus)/2000` and check the mean against
  112 to obtain a reference for the deformed tangent compliance.
- Cases 116/117: `R1` rotation, torque, and torsional stresses and strains.
- Cases 119/120: load--displacement curve, number of cycles, cutbacks, and the
  final residual before/after the limit point.
- Cases 121/122: relative difference in displacement and reactions for the
  small-rotation limit of a transversely loaded truss.
- Cases 123--127: displacements of every free node, reactions at every support,
  and stresses in at least one element for each PROD; for 124/126, also verify
  that the final prescribed displacement is exactly `-4.0/+5.0` and record its
  conjugate reaction; for 127, record every adaptive load factor.

The campaign does not request a tangent DMIG for each iteration: MSC Nastran
does not provide a standard entry equivalent to the linear export used for
shells. If the preceding balance is insufficient to identify the formulation,
add a dedicated DMAP deck without modifying these campaign controls.
