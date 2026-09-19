# CTRIA3 PSHELL factors: result validation

Each model is a four-CTRIA3, six-node cantilever. Geometry, materials, loads,
and supports remain fixed; only the PSHELL factors change:

| Case | `12I/T^3` | `TS/T` |
|---|---:|---:|
| R01 | 1.00 | 0.25 |
| R02 | 1.00 | 1.25 |
| R03 | 0.40 | 0.50 |

Every deck sets `AUTOSPC=NO`, constrains the left edge explicitly, and uses
`K6ROT=100` for drilling rotation. Combined in-plane tension/shear, normal
forces, and a nodal moment exercise displacement, reaction, membrane, moment,
shear, and outer-fibre recovery.

References come from MSC Nastran Student Edition 2026.1. Recovered-shear
tolerance remains separate from other quantities so correct factor parsing is
not conflated with the known CTRIA3 centroidal-recovery difference.

All three decks completed on 3 September 2026 with `END OF JOB` and no
`FATAL`. Selected-value maxima are 1.10% for displacement, 0.61% for
reaction, 4.02% for moment, and 4.41% for fibre stress; membrane error remains
below `2.7e-8` relative. Significant shear differs by at most 0.87%, while
near-zero components differ by less than `8.3e-3 N/mm` absolute.
