# Factores PSHELL del CTRIA3: validación de resultados

Los tres modelos forman un voladizo de cuatro CTRIA3 y seis nodos. Comparten
geometría, materiales, cargas y empotramiento; sólo cambian los factores
`TS/T` y `12I/T^3` de PSHELL:

| Caso | `12I/T^3` | `TS/T` |
|---|---:|---:|
| R01 | 1.00 | 0.25 |
| R02 | 1.00 | 1.25 |
| R03 | 0.40 | 0.50 |

Todos declaran `AUTOSPC=NO`, fijan explícitamente el borde izquierdo y usan
`K6ROT=100` para el giro de perforación. La carga combina tracción en el plano,
cortante en el plano, fuerzas normales y un momento nodal. De este modo se
ejercitan desplazamientos, reacciones, membrana, momentos, cortantes y
tensiones de las fibras exteriores.

Las referencias proceden de MSC Nastran Student Edition 2026.1. La tolerancia
del cortante recuperado se mantiene separada de las restantes magnitudes para
no confundir la lectura correcta de los factores con la discrepancia conocida
de recuperación centroidal del CTRIA3.

Los tres decks se ejecutaron el 3 de septiembre de 2026 sin mensajes `FATAL`
y terminaron con `END OF JOB`. En los valores seleccionados, los máximos son
`1.10 %` para desplazamientos, `0.61 %` para reacciones, `4.02 %` para
momentos y `4.41 %` para fibras; membrana queda por debajo de `2.7e-8`
relativo. El cortante significativo difiere como máximo un `0.87 %` y las
componentes casi nulas menos de `8.3e-3 N/mm` en valor absoluto.
