# CTRIA3 PSHELL material semantics

These seven decks isolate the activation of `MID1`, `MID2`, and `MID3`
in a single CTRIA3. They share geometry, thickness, and three distinct
materials; `K6ROT=0` prevents drilling stabilization from contaminating the
comparison.

The models were run with MSC Nastran Student Edition 2026.1 on 3 September
2026. Cases 01, 02, 03, 04, and 07 completed without `FATAL` messages and
wrote the `KAAX` matrix to PCH. Cases 05 and 06 deliberately end with
`USER FATAL MESSAGE 9994 (BULAEM)` and
`MID3 = 0 when MID2 = 0`, so they do not produce PCH output.

| Case | Membrane | Plate | Result |
|---|---|---|---|
| 01: `MID1` only | active | inactive | valid |
| 02: `MID1+MID2` | active | active | valid |
| 03: `MID1+MID2+MID3` | active | active | valid |
| 04: `MID2` only | inactive | active | valid |
| 05: `MID3` only | -- | -- | `FATAL 9994` |
| 06: `MID1+MID3` | -- | -- | `FATAL 9994` |
| 07: `MID2+MID3` | inactive | active | valid |

The references establish independent `MID1` and `MID2` assignment and show
that `MID3` is valid only when `MID2` exists. They are not used to
recalibrate the solver's Reissner--Mindlin formulation: these tests validate
block activation semantics, not term-by-term plate-block equivalence.
