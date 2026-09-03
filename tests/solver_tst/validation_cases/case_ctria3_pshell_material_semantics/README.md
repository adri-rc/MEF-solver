# Semántica de materiales PSHELL en CTRIA3

Estos siete decks aíslan la activación de `MID1`, `MID2` y `MID3` en un
único CTRIA3. Todos comparten geometría, espesor y tres materiales distintos;
`K6ROT=0` evita que la estabilización de giro normal contamine la comparación.

MSC Nastran Student Edition 2026.1 se ejecutó el 3 de septiembre de 2026. Los
casos 01, 02, 03, 04 y 07 terminaron sin mensajes `FATAL` y escribieron la
matriz `KAAX` en sus PCH. Los casos 05 y 06 terminaron deliberadamente con
`USER FATAL MESSAGE 9994 (BULAEM)` y la indicación
`MID3 = 0 when MID2 = 0`; por eso no producen PCH.

Resultados estructurales de las matrices MSC:

| Caso | Membrana | Placa | Resultado |
|---|---|---|---|
| 01: sólo `MID1` | activa | inactiva | válido |
| 02: `MID1+MID2` | activa | activa | válido |
| 03: `MID1+MID2+MID3` | activa | activa | válido |
| 04: sólo `MID2` | inactiva | activa | válido |
| 05: sólo `MID3` | -- | -- | `FATAL 9994` |
| 06: `MID1+MID3` | -- | -- | `FATAL 9994` |
| 07: `MID2+MID3` | inactiva | activa | válido |

La referencia demuestra la asignación independiente de `MID1` y `MID2`, y
que `MID3` sólo es admisible si existe `MID2`. No se emplea para recalibrar la
formulación Reissner--Mindlin del solver: los tests de este caso comprueban la
semántica de activación, no equivalencia término a término del bloque de placa.
