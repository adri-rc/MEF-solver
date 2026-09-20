# SOL 400 closure: CBAR fibres and mixed CBAR/CROD models

## Objective and status

Cases 161--176 form the final geometric-nonlinear validation campaign. The
DAT decks pass preprocessing, rank checks, and the local solver. Their MSC
Nastran F06/PCH runs finish without a `FATAL` message, reach load factor 1.0,
and their final references are transcribed into the CBAR and mixed end-to-end
tests.

References for cases 161--166 are in
`test_nonlinear_cbar_validation_end_to_end.py`. Cases 167--176 are kept
separately in `test_nonlinear_mixed_validation_end_to_end.py`; each compares
at least one CBAR result and one CROD result in addition to the nodal response.

All models use consistent N-mm-MPa units, `SOL 400`, `ANALYSIS=NLSTAT`,
`PARAM,LGDISP,1`, conservative nodal loads, and fixed load increments. Each
Bulk Data card is preceded by a comment identifying its fields. The decks
request `OUTCTRL=(SOLUTION,INTERM)` and the `NLSTRESS`, `BOUTPUT`, `ELFORCE`,
`STRAIN`, displacement, reaction, and energy tables.

## CBAR fibre-stress diagnostics: cases 161--166

This first family isolates the operations that can change the sign convention
or the axis on which MSC evaluates the C--F recovery points.

| Case | Model | Identified behaviour |
|---:|---|---|
| 161 | one CBAR along X, symmetric C--F points | sign reference under combined axial load, biaxial bending, and torsion |
| 162 | one CBAR inclined in 3D | transport of the fibres with the corotated frame |
| 163 | asymmetric section and `I12 != 0` | biaxial coupling without inferring the centroid by averaging C--F |
| 164 | `G0` orientation and loads in `CORD2R` | composition of orientation, load transformation, and finite rotation |
| 165 | two-CBAR spatial elbow, two MAT1/PBAR definitions | A/B end convention and continuity between different frames |
| 166 | three CBAR elements, `LROFF`, and imposed rotation | fibres with offsets, multiple elements, and a nonzero rotational SPC |

Closing the recovery requires more than the maximum and minimum values. The
last state and at least two intermediate states must retain the following for
every CBAR:

- C, D, E, and F normal stresses at ends A and B;
- axial force, local moments 1/2, shear forces, and torque;
- axial strain and the available curvatures or section strains;
- global nodal displacements and rotations;
- an unambiguous load-factor identifier for every result block.

Cases 161 and 162 establish the signs first. Case 163 determines how `I12` and
the asymmetric points operate. Cases 164--166 must be interpreted only after
those results, so a transformation error cannot accidentally compensate for a
sign error.

## Mixed CBAR/CROD models: cases 167--176

| Case | Nodes | CBAR | CROD | Main configuration |
|---:|---:|---:|---:|---|
| 167 | 3 | 1 | 1 | cantilever with tie, FORCE, and MOMENT |
| 168 | 4 | 3 | 1 | braced frame with two materials |
| 169 | 6 | 2 | 3 | spatial mast with three guy rods and distinct CROD properties |
| 170 | 6 | 5 | 4 | two-bay cross-braced frame with several nodal loads |
| 171 | 8 | 8 | 6 | spatial platform with columns, ring members, and diagonals |
| 172 | 9 | 10 | 8 | two-storey frame with imposed vertical settlement |
| 173 | 9 | 12 | 10 | pyramidal roof with several sections and `LC/LS` integration |
| 174 | 12 | 16 | 8 | two-level space frame with loads in `CORD2R` and imposed displacement |
| 175 | 10 | 13 | 15 | spatial building with offsets, three materials, and imposed displacement |
| 176 | 12 | 20 | 10 | three-level braced tower with several FORCE/MOMENT loads and imposed settlement |

The sequence deliberately increases the model size to 30 elements and
combines `MAT1`, `PBAR/PBARN1`, `PROD`, `CBAR`, `CROD`, `FORCE`, `MOMENT`,
`LOAD`, `SPC1`, `SPC`, `SPCADD`, load coordinate systems, offsets, and
different properties. It does not use `K6ROT`: that parameter stabilizes shell
rotations and is not part of the CROD/CBAR formulation.

## Final integration model

The numbered campaign culminates in the unnumbered structural model retained
in the same mixed end-to-end test:

| Case | GRID | CBAR | CROD | Action and purpose |
|---|---:|---:|---:|---|
| `case_nonlinear_test_final_updated` | 202 | 144 | 428 | three-span bridge, explicit bases, 44 FORCE cards, 3 MOMENT cards, and 20 increments |

The final model correlates the complete result families with L2 errors of
approximately `10^-5` for displacements, reactions, and axial results. CBAR
fibres in the continuous deck require separate interpretation: `1.673 %` L2
error and `11.57 %` maximum significant error. The test therefore retains
separate tolerances for the global response, reactions, centroidal states,
and `stress_points_*` paths; the fibre envelope does not apply to the other
results.

## Execution order and audit

Run cases 161--163 first. If `NLSTRESS/BOUTPUT` results are present for CBAR
and every case reaches load factor 1.0, continue with 164--166. Then run
167--171 and finally 172--176. Do not change the job name: the directory,
case number, and main SID agree.

Retain at least the `.f06` and `.pch` files alongside each DAT deck. Also keep
the `.itr` and `.sts` files during the initial audit because they can clarify
an incremental failure; decide afterwards whether they warrant versioning.

Check the following in every F06 file:

1. no `USER FATAL MESSAGE` or `SYSTEM FATAL MESSAGE`;
2. normal termination and a final load factor of 1.0;
3. the actual number of converged increments and any bisections;
4. `DISPLACEMENT VECTOR` and `FORCES OF SINGLE-POINT CONSTRAINT`;
5. nonlinear CBAR and CROD tables, including C--F values in cases 161--166;
6. equilibrium among external loads, reactions, and element actions;
7. intermediate results requested by `OUTCTRL=(SOLUTION,INTERM)`.

An F06 file that terminates normally but contains no CBAR results does not
validate the model: confirm that MSC recognizes `PBARN1`. If a complex case
fails, do not rerun the entire campaign; retain the failed F06 file to
determine whether the cause is convergence, card formatting, or an unsupported
combination.
