# Estado y validación

## Línea base reproducida

Fecha: **13 de agosto de 2026**. Entorno: `.venv` del repositorio.

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

Resultado:

```text
316 passed, 5 xfailed in 2.27s
```

También se ejecutó `python -m compileall -q src tests main.py` sin errores. Los
cinco `xfail` son discrepancias CQUAD4 3D/multielemento conocidas, declaradas con
`strict=True`; no hay fallos inesperados.

## Estrategia

La validación distingue geometría, matrices constitutivas, rigidez elemental,
transformación global, ensamblaje, condiciones de contorno, solución,
reacciones y recuperación. Se usan pruebas analíticas, invariantes, *patch
tests* y comparaciones con F06/PCH/CSV de Nastran.

## Cobertura resumida

- preproceso: `GRID`, materiales, propiedades, elementos, `INCLUDE`, Case
  Control, combinaciones `LOAD`, `SPCADD` y sistemas de coordenadas;
- núcleo: mapeo de GDL, cuadratura, transformaciones, ensamblaje disperso,
  singularidades, desplazamientos impuestos y reacciones;
- CROD/CBAR: rigidez, transformaciones 3D, fuerzas y tensiones;
- CTRIA3: membrana, flexión, cortante, orientaciones 3D y comparación de
  matrices globales;
- CQUAD4: dos formulaciones separadas, postproceso, `K6ROT`, modos rígidos,
  comparación DMIG de los casos 51-55, 65, 67 y 69-75, y respuesta E2E de los casos
  56-64, 66 y 68 con referencias F06 escritas directamente en los tests.

## CQUAD4 frente a Nastran

| Caso | Comparación | Error de Frobenius | Criterio |
|---|---|---:|---:|
| 51 | bloque de membrana | `0.152382 %` | `< 1 %` |
| 52 | bloque de flexión | `0.152382 %` | `< 1 %` |
| 53 | bloque activo combinado | `0.372560 %` | `< 1 %` |
| 54 | matriz global completa `24 x 24` | `2.113059 %` | `< 3 %` |
| 55 | matriz `24 x 24`, 3D y `K6ROT=100` | `1.077026 %` | `< 3 %` |
| 65 | cortante aislado, `K6ROT=0` | `5.74e-9 %` | `< 1e-6 %` |
| 67 | MID1/MID2/MID3 y `K6ROT=100` | `0.000298 %` | `< 0.001 %` |
| 69-71, 73 | MID3/G3/TST rectangulares | `< 7.3e-9 %` | `< 1e-6 %` |
| 72 | MID3 vacío, `K6ROT=100` | `0.000121 %` | `< 0.001 %` |
| 74-75 | geometría sesgada, con/sin MID3 | `0.72145 %` | `< 1 %` |

El caso 54 es inclinado y alabeado: todos los GDL globales reciben
contribuciones. La comparación no elimina filas ni columnas y verifica la
corrección de fuerzas y momentos entre el plano medio y los nodos reales. El
resultado cumple el umbral del proyecto, pero debe describirse como
coincidencia aproximada, no como identidad numérica.

Los casos 65/69 y 67/72 forman pares con y sin `MID3`; 70, 71 y 73 cambian
solo `G3` o `TS/T`. Todos ellos verifican tanto el error global como los
términos individuales. Los F06 66 y 68 confirman que desplazamientos, cargas,
reacciones y resultantes recuperados coinciden dentro de sus tolerancias. El
error casi idéntico de 74 y 75 muestra que el `0.72 %` residual pertenece a la
aproximación geométrica irregular y no a la flexibilidad física de `MID3`.

El detalle de formulación, resultados y fuentes está en [cquad4.md](cquad4.md).

## Criterios para nuevas capacidades

- simetría de `K` próxima a precisión de máquina;
- invariancia ante traslación y rotación rígida;
- número y forma de los modos de cuerpo rígido esperados;
- energía no negativa;
- *patch test* y solución analítica en geometrías regulares;
- sensibilidad controlada a distorsión y espesor;
- comparación independiente con tolerancias justificadas;
- convergencia con refinamiento de malla;
- equilibrio entre cargas y reacciones.

Las entradas próximas a cero deben evaluarse con tolerancia absoluta. Para una
matriz completa se informa además del máximo absoluto y del peor término. El
error relativo de Frobenius sigue siendo un resumen global, pero no sustituye
la comparación por bloques ni el chequeo término a término: una submatriz grande
puede ocultar errores importantes en otra más pequeña.

## Trabajo pendiente de validación

1. Reducir el error geométrico residual de los CQUAD4 sesgados y alabeados.
2. Ampliar los estudios de convergencia y *patch tests* multielemento.
3. Documentar la versión exacta de Nastran que genera cada referencia.
4. Consolidar los casos CTRIA3 con diferencias pequeñas aún conocidas.
5. Incorporar CI para publicar automáticamente el resultado de la suite.
