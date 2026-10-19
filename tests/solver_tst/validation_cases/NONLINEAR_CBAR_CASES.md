# MSC Nastran SOL 400 campaign for CBAR

## Purpose

Decks 128--166 provide the external references for the corotational CBAR.
The campaign contains 39 cases: the 21 initial cases 128--148, five larger
models 149--153, seven interpolation diagnostics 154--160, and six fibre
recovery diagnostics 161--166.

Each folder contains one `.dat` file. Every Bulk Data card is preceded by
comments naming its fields so the deck can be audited without consulting
another file. The models use consistent N-mm-MPa units.

## Case matrix

| Case | Model or variation | Main check |
|---:|---|---|
| 128 | small axial CBAR, `LGDISP=-1` | SOL 400 reference in the linear limit |
| 129 | small bending, `LGDISP=-1` | linear bending and sign reference |
| 130 | finite axial tension | extension, axial force, and material stiffness |
| 131 | subcritical compression | geometric stiffness in compression |
| 132 | large force in local Y | bending about `I1`, rotation, and P-delta effect |
| 133 | large force in local Z | bending about `I2` and sign convention |
| 134 | pure moment about local Z | constant curvature in the first plane |
| 135 | pure moment about local Y | constant curvature in the second plane |
| 136 | pure torsion | `GJ`, relative twist, and torque |
| 137 | base combined load | finite axial/biaxial/torsional state |
| 138 | base plus a `+100 N` Y probe | positive tangent-compliance difference |
| 139 | base plus a `-100 N` Y probe | central tangent-compliance difference |
| 140 | inclined 3D CBAR, `G0`, loads in `CID` | objectivity, orientation, and load transformation |
| 141 | asymmetric section, `I12 != 0` | coupling between the two bending planes |
| 142 | short beam, `K1/K2`, `INT2=LS` | one-point linear shear interpolation |
| 143 | prescribed transverse displacement | nonzero SPC and reaction without external load |
| 144 | prescribed displacement/rotation and FORCE | multiple nonzero SPCs combined with a load |
| 145 | two elements with different MAT1/PBAR | assembly and interface between properties |
| 146 | cantilever refined to four elements | mesh convergence relative to case 132 |
| 147 | rectangular frame with `PB=6` | stable condensation of an end release |
| 148 | spatial frame with offsets | integration of materials, loads, SPC, and `LROFF` |
| 149 | two-story planar frame, 15 CBAR | three-property assembly and global P-delta effect |
| 150 | two-level spatial frame, 20 CBAR | rotated loads, 3D bracing, and adaptive stepping |
| 151 | semicircular arch, 20 CBAR | curved geometry, four properties, and imposed settlement |
| 152 | braced spatial building, 20 CBAR | `LROFF` offsets, combined load, and prescribed displacement |
| 153 | 3D zigzag boom, 15 CBAR | five sections, three materials, `I12`, shear, and imposed rotation |

## Closing `INT2=LC/LS`: cases 154--160

The initial implementation showed that case 142 could not be reproduced by
replacing MSC's `LS` integration with the closed Timoshenko matrix of `PBAR`.
Seven minimal decks isolate the correct law:

| Case | INT2 | K1/K2 | Load | Comparison |
|---:|---|---|---|---|
| 154 | LC | 0.75/0.50 | Y + 0.5 Z | direct counterpart to 142 |
| 155 | LS | 0.75/0.50 | Y | shear/K1 plane |
| 156 | LC | 0.75/0.50 | Y | counterpart to 155 |
| 157 | LS | 0.75/0.50 | Z | shear/K2 plane |
| 158 | LC | 0.75/0.50 | Z | counterpart to 157 |
| 159 | LS | 0/0 | Y + 0.5 Z | LS effect without PBAR factors |
| 160 | LC | 0/0 | Y + 0.5 Z | counterpart to 159 |

All use the same 30 increments and `OUTCTRL=(SOLUTION,INTERM)`. Their F06
files end without `FATAL` and reach load factor 1.0. Cases 154/160 demonstrate
that `LC` ignores `K1/K2`; 155/157 identify the linear shear law; and 159
identifies the unit factor used by `LS` when `K1/K2=0`.

Final references are transcribed in
`test_nonlinear_cbar_validation_end_to_end.py`. States 0.2, 0.5, and 0.8 from
155--159 are also retained in
`test_nonlinear_cbar_intermediate_states.py`. No test reads Nastran output
at runtime.

## Closing CBAR fibre recovery: cases 161--166

These six decks progressively isolate operations that can alter the sign
convention or the axes at which MSC evaluates recovery points C--F:

| Case | Model | Main check |
|---:|---|---|
| 161 | one X-aligned CBAR with symmetric C--F points | combined axial/biaxial bending/torsion signs |
| 162 | one CBAR inclined in 3D | transport of fibres with the corotated triad |
| 163 | asymmetric section and `I12 != 0` | biaxial coupling without inferring the centroid from C--F |
| 164 | `G0` orientation and loads in `CORD2R` | orientation, load transformation, and finite rotation |
| 165 | two-member spatial elbow with two MAT1/PBAR pairs | A/B end convention and continuity between distinct triads |
| 166 | three CBAR, `LROFF`, and prescribed rotation | fibres with offsets, multiple elements, and a nonzero SPC |

Their F06/PCH references end without `FATAL`, reach load factor 1.0, and are
transcribed into `test_nonlinear_cbar_validation_end_to_end.py`. The end-to-end
checks include C--F values at both ends as well as nodal response, centroidal
axial state, and selected end actions. Nearly neutral fibres and the moment
distribution in the second segment of case 165 and third segment of case 166
remain explicit diagnostics rather than broadly tolerated or expected
failures. The wider 161--176 closure campaign is described in
[`NONLINEAR_FINAL_CLOSURE_CASES.md`](NONLINEAR_FINAL_CLOSURE_CASES.md).

## Audit of the first run (2026-08-20)

The original F06 files for cases 128--148 ended normally: they contain no
`USER/SYSTEM FATAL` messages, complete step 1, and reach load factor 1.0.
However, that run is not a valid nonlinear CBAR reference. The original DAT
files contained only `PBAR`; SOL 400 requires an associated `PBARN1` to enable
the advanced CBAR formulation. The response of perturbed cases 137--139 was
exactly linear, and the F06/OP2 files contained no CBAR element results,
revealing the problem.

DAT files 128--148 were corrected with `PBARN1` (`BEH2=BAR`, `INT2=LC`, except
for `INT2=LS` in case 142), `NLSTRESS`/`BOUTPUT`, and explicit `STRAIN` and
`ELFORCE` requests. The corrected run closed 23 of the 26 cases at load factor
1.0. The remaining three were audited as follows:

- Cases 140 and 150 solved the first increment, but intermediate export ended
  with `USER FATAL 2025`: SOL 400 incremental files did not retain the
  `CORD2R` used by `GRID.CD`. Their DAT files now output results in the basic
  system while retaining `FORCE/MOMENT.CID` to validate the input
  transformation.
- Case 148 converged at every increment, but adaptive control exhausted
  `MAXINC=400` at load factor 0.775568. Its DAT now uses 100 fixed increments.

Cases 140, 148, and 150 were subsequently rerun. Their current F06/PCH files
reach load factor 1.0 and replace the failed versions described above.

Cases 137--139 must be run without modifying any other input. At the final
state, the Y-direction tangent compliance at the free end can be approximated
as:

```text
d(uy)/d(Fy) ~= (uy_138 - uy_139) / 200 N
```

Case 137 checks that the mean of the two perturbed solutions remains close to
the base state.

## Intermediate references for implementing the tangent

All decks request `NLOPRM OUTCTRL=(SOLUTION,INTERM) NLDBG=ADVDBG`. For each
converged increment, Nastran writes displacements, loads, reactions,
work/energy, nonlinear CBAR stresses and strains, and the iteration and
residual history. Cases 132 and 133 isolate the two bending planes; 136
isolates torsion; 145--146 isolate assembly and refinement; and 149--153 check
larger assemblies.

Cases 137--139 are especially useful because they share exactly 40 fixed
increments. Because Nastran also multiplies the probe by load factor `lambda`,
the following can be evaluated at any common nonzero increment:

```text
d(uy)/d(Fy) | lambda ~= (uy_138(lambda) - uy_139(lambda)) / (200*lambda N)
```

This provides an external tangent-compliance reference along the whole path,
not only at the final state. The F06 does not print the complete global
tangent matrix or, with `PBARN1`, the classic `FORCES IN BAR ELEMENTS` table.
Debugging therefore combines this probe with nodal reactions and the
`NLSTRESS/BOUTPUT` tables. The implementation uses these probes in the
intermediate test. Diagnostics 155--159 also provide interpolation-specific
intermediate references that detect an `LC/LS` regression without relying
only on the final state.

## Evidence retained from the runs

At least the `.f06` and `.pch` files should be returned with each `.dat` for
an initial audit. The `.pch` simplifies vector transcription; the `.f06`
retains convergence, warnings, and element tables. Any additional large
artifact must provide information that cannot be reconstructed from these
sources before it is added to the repository.

Check the following in every F06:

- `FATAL MESSAGE`, `USER FATAL MESSAGE`, `WARNING MESSAGE`, and job
  termination;
- `NON-LINEAR ITERATION MODULE OUTPUT` and the final load percentage/factor;
- `D I S P L A C E M E N T   V E C T O R`;
- `F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T`;
- nonlinear CBAR tables with `STRESS`, `TOTAL STRAIN`, and C--F stations;
- strain energy and any intermediate table requested by
  `OUTCTRL=(SOLUTION,INTERM)`.

SOL 400 CBAR output uses corotated local axes. When transcribing it, record
whether a column is an internal section resultant or an element nodal action;
its sign need not match the code's `forces_n1/forces_n2` directly.

## Recommended run order

Run 128--136 first and verify that `BOUTPUT` or `NLSTRESS` results appear for
CBAR. If they end normally, continue with 137--146, then 147--148, 149--153,
154--160, and finally 161--166. An isolated failure in 141, 142, 147, 148,
151, or 152 points respectively to `I12`, `LS` integration, releases, offsets,
an imposed settlement, or offsets in a larger mesh, without invalidating the
basic cases.

Do not rename jobs or their folders: the case number matches the `SPC`, `LOAD`,
and `NLSTEP` SIDs, except for auxiliary sets documented in the deck itself.
