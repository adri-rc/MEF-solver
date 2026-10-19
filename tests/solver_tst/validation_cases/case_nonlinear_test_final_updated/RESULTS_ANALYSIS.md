# Results analysis and correlation with MSC Nastran

## Reference and method

The independent reference is the final `LOAD STEP = 1.00000E+00` block from
`tfm_final_nonlinear_bridge_updated.f06`, generated with MSC Nastran Student
Edition 2026.1 on 25 August 2026. Displacements and reactions were also checked
against the final PCH block. MSC and MEF-solver solve the same DAT without an
intermediate Patran regeneration step.

For a complete vector or result family, the relative L2 error is

```text
epsilon_L2 = ||x_local - x_MSC||_2 / ||x_MSC||_2
```

Componentwise relative errors are interpreted only when the reference has
sufficient magnitude. The absolute difference is also reported near zero,
where a percentage ratio can be misleading.

As in the rest of the mixed-element campaign, centroidal `CBAR` axial stresses
are reconstructed by averaging the C/D/E/F stresses at both ends. Fibre
stresses are compared directly point by point. Current length is computed from
the original coordinates and MSC's final `DISPLACEMENT VECTOR`. Test references
are literals written in the code: the end-to-end test does not open the F06 or
PCH at runtime.

## Convergence

| Quantity | MSC Nastran | MEF-solver |
|---|---:|---:|
| Final load factor | 1.0 | 1.0 |
| Converged increments | 20 | 20 |
| Newton iterations | 40 | 40 |
| Cutbacks/bisections | 0 | 0 |
| MSC final residual ratio | 6.369e-8 | — |
| Strain energy [Nmm] | 1.092467e7 | 1.0924822116e7 |

MSC prints two iterations for every increment, `JOB CONVERGES FOR THE CURRENT
STEP`, and `END OF JOB`. There are no `USER FATAL`, `SYSTEM FATAL`,
`USER WARNING`, or `SYSTEM WARNING` messages. MSC's final absolute residual is
`2.919e-2 N`; the normalized ratio remains below the `3e-7` criterion. Energy
differs by a relative `1.392e-5`. Residual values retain the precision visible
in the final F06.

## Global response

The complete comparison, rather than only the points included in the
end-to-end test, gives:

| Complete family | Relative L2 error |
|---|---:|
| Translations, 202 GRID | 3.280e-5 |
| Rotations, 202 GRID | 2.416e-5 |
| Force reactions | 2.898e-5 |
| Moment reactions | 4.250e-5 |
| Strain energy | 1.392e-5 |

The largest translation norm occurs at `GRID 1127`:

| Code | T1 [mm] | T2 [mm] | T3 [mm] | Norm [mm] |
|---|---:|---:|---:|---:|
| MSC | 1.36855459 | 1.83992660 | -15.23815823 | 15.40972867 |
| MEF-solver | 1.36860658 | 1.84006768 | -15.23830206 | 15.40989236 |

MSC's largest rotation is `R1=0.04409650 rad` at `GRID 1015`. The
`-2/-6/-10 mm` settlements are recovered exactly at the GRID points of the
three affected bases. No drift of the imposed values is detected.

## CROD results

The correlation across the axial stresses of all 428 `CROD` elements gives
`epsilon_L2=2.125e-5`. The largest absolute difference is `0.006625 MPa` at
`EID 20610`. The largest raw percentage, `0.780 %`, occurs at `EID 20122`,
whose reference is nearly zero (`-0.01918 MPa`); it therefore does not govern
the physical assessment. When the comparison is restricted to
`|sigma_MSC|>1 MPa`, the maximum relative error is `0.0537 %`.

| Complete CROD result | Relative L2 error |
|---|---:|
| Axial stress | 2.125056e-5 |
| Axial strain | 2.125684e-5 |
| Axial force | 2.123120e-5 |
| Current length | 8.367439e-9 |

At `EID 20610`, which governs the absolute differences, the reconstructed
force is `-67029.495 N` in MSC and `-67033.8012 N` locally; the lengths are
`3818.973969` and `3818.973840 mm`.

The largest absolute axial stress is `+220.4665 MPa` in network hanger
`EID 20929`, `PID 207`. Its sign denotes tension. The model material is elastic
without a yield limit: this value compares codes rather than accepting a steel
grade or sizing the hanger. Other hangers can appear in compression because
`CROD` is bilateral and does not implement slackening.

## CBAR results

The centroidal axial stress of all 144 `CBAR` elements gives
`epsilon_L2=1.379e-5`; the largest absolute difference is `0.000661 MPa` at
`EID 12111`. The axial response and global force transfer therefore correlate
at the same level as the displacements.

| Complete CBAR result | Relative L2 error |
|---|---:|
| Centroidal axial stress | 1.379175e-5 |
| Axial strain | 1.379135e-5 |
| Axial force | 1.627389e-5 |
| Current length | 1.047300e-9 |

At `EID 12111`, axial force is `-193312.564 N` in MSC and `-193323.145 N`
locally; the length difference is only `-4.885e-6 mm`.

C/D/E/F recovery is more sensitive. Comparing all 1152 fibres at both ends:

| CBAR group | Fibre relative L2 error |
|---|---:|
| All CBAR elements | 1.673e-2 |
| Deck girders, PID 101 | 4.212e-2 |
| Arches, PID 103 | 8.779e-4 |
| Columns, PID 105 | 2.294e-4 |
| Foundation legs, PID 106 | 2.171e-4 |

The largest absolute difference is `2.529 MPa` at `EID 10014`, end B, fibre
E/P2. Among fibres whose reference exceeds `10 MPa`, the worst percentage is
`11.57 %` at `EID 10015`, end B, fibre D/P1: MSC prints `19.48201 MPa` and
MEF-solver recovers `21.73605 MPa`.

The discrepancy does not appear in centroidal axial results and is concentrated
in the continuous deck girders, where end moments and the corotational fibre
recovery convention have greater influence. Neither stiffness nor formulation
has been recalibrated to force this fit. The end-to-end test uses a separate
tolerance for the `stress_points_*` paths, so the fibre envelope does not weaken
the displacement, reaction, energy, axial-force, or `CROD` checks.

The maximum MSC normal fibre stress is `121.1869 MPa` at `EID 12110`, at the
end associated with `GRID 3021`, fibre F. It is below the maximum axial stress
in the hangers, although both values belong to an elastic code comparison and
not to a resistance check.

## Geometric-nonlinearity signal

An otherwise identical local control was solved with only the evaluation
changed to `LGDISP=-1`; the preserved DAT was not modified. At `GRID 1015/T2`:

| Evaluation | T2 [mm] |
|---|---:|
| `LGDISP=1` | 2.20017237795 |
| `LGDISP=-1` | 2.15731057757 |

The change is `+1.9868164 %`. The maximum norm changes from `15.40989236 mm`
to `15.45020968 mm`, and the L2 error between all free translations in the two
controls is `0.3092466 %`. The updated geometry therefore exercises a
measurable nonlinear signal without approaching a limit point or requiring
cutbacks.

## Regression selection

The entry in `test_nonlinear_mixed_validation_end_to_end.py` retains a sample
distributed across all three spans, both sides, piers, and foundations:

- deck, arch, inner-line, and base GRID points, including displacements and
  reactions;
- `CBAR` elements from girders, arches, columns, and foundation legs;
- `CROD` elements from inner girders, cross-members, hangers, bracing, base
  beams, and network hangers;
- final energy and all 20 increments;
- the `EID 10014/10015` fibres that govern the recovery differences.

The test contains 230 numerical reference values: 54 displacement components,
24 reaction components, 76 results from 15 CBAR elements (including 16 fibre
values), 75 results from 15 CROD elements, and the energy. It also checks 30
element types, the final factor, and all 20 increments. None of these values is
obtained by reading external files during the test.

The common tolerance remains separate from the fibre tolerance. This prevents
the `11.57 %` envelope, which belongs to only one recovered quantity, from
hiding a global regression.

The entry's envelopes are separated by result family: `0.1 %` for
displacements and energy, `0.6 %` for reactions (governed by a vertical
reaction of only `236.6 N`), `1 %` for centroidal element states (governed by
the nearly unloaded `CROD 20122`), and `12 %` exclusively for the governing
`stress_points_*` fibres. `abs_tol=1e-5` does not mask any of those governing
components.

## Visual postprocessing evidence

Patran 2026.1 imported the DAT and all 20 OP2 states. Ten real captures are
preserved for geometry, supports and settlements, loads, displacement views,
and the final deformed configuration. The final view overlays the original blue
geometry with the amplified white response:

![Final deformed configuration in Patran](tfm_final_nonlinear_bridge_updated_patran_deformed_final.png)

The complete set and the explanation of which SPC entries are visible in each
view are in the case README's [Patran views](README.md#patran-views) section.

No stress contour was obtained because the results database that Patran builds
from the OP2 offers only loads, displacements, and constraint forces. The
consolidated OP2 contains no `CROD/CBAR` force, stress, or strain tables,
although those tables are printed in the F06 and support the comparison above.
This export limitation is documented instead of presenting a contour that does
not originate from the actual element results.

## Observations and interpretation limits

- The four bases begin at different elevations; they represent different
  abutment/pier heights and staggered ground, not one horizontal soil plane.
- Arch, deck, and pier share a node at the springings; there are no bearings,
  releases, or unilateral contact.
- Inner lines, hangers, and `CROD` bracing provide no bending stiffness and do
  not buckle. Their inactive rotations are restrained explicitly.
- Automatic self-weight, distributed loads, plasticity, fatigue, buckling,
  cable slackening, contact, and soil-structure interaction are omitted.
- Patran stress views would have to separate `CBAR` and `CROD`; an entirely
  1D mesh has no continuous von Mises surface field.

Subject to these qualifications, the case closes the validation of mixed
assembly, geometrically nonlinear kinematics, and axial recovery at a scale of
approximately 200 nodes. Deck fibre recovery remains explicitly quantified and
monitored rather than presented as an exact match.
