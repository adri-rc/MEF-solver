# CTRIA3 PSHELL factors: stiffness matrices

These four one-element models isolate the two scalar `PSHELL` factors. K01
and K02 vary `TS/T` while keeping `12I/T^3=1`; K03 and K04 vary
`12I/T^3` while keeping `TS/T=0.833333`.

Every deck sets `AUTOSPC=NO`, `K6ROT=0`, and requests the `KAAX` matrix
through `EXTOUT=DMIGPCH`. Preserved PCH files and their derived CSV matrices
provide the MSC references.

MSC Nastran Student Edition 2026.1 completed all four decks on 3 September
2026 with `END OF JOB` and no `FATAL`. Solver/MSC relative Frobenius errors
are `9.96e-5` and `4.98e-4` for K01/K02 and below `1.9e-10` for
K03/K04. Representative K01/K02 terms remain below `1.23e-3` relative
error.

K01/K02 deliberately use a flexible `MID3` so its contribution can be
isolated without turning this check into a plate-formulation recalibration.
K03/K04 use a nearly zero `MID3` to isolate bending scaling. The campaign
validates application of factors read from the deck, not general equivalence
between the local CTRIA3 plate and MSC's native formulation.
