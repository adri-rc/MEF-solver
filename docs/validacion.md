# Estado y validación

## Línea base reproducida

Fecha: **13 de agosto de 2026**. Entorno: `.venv` del repositorio.

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

Resultado:

```text
333 passed, 4 xfailed in 3.33s
```

También se ejecutó `python -m compileall -q src tests main.py` sin errores. Los
cuatro `xfail` son discrepancias CQUAD4 3D/multielemento conocidas, declaradas con
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
  comparación DMIG de los casos 51-55, 65, 67, 69-75, 81, 83 y 89, y
  respuesta E2E de los casos 56-64, 66, 68, 76-80, 82 y 84-88 con referencias F06 escritas
  directamente en los tests.

## CQUAD4 frente a Nastran

| Caso | Comparación | Error de Frobenius | Criterio |
|---|---|---:|---:|
| 51 | matriz completa; `MID2` vacío deja placa nula | `0.039670 %` | `< 1 %` |
| 52 | bloque de flexión | `0.039670 %` | `< 1 %` |
| 53 | bloque activo combinado | `0.023385 %` | `< 1 %` |
| 54 | matriz global completa `24 x 24` | `2.079993 %` | `< 3 %` |
| 55 | matriz `24 x 24`, 3D y `K6ROT=100` | `0.028215 %` | `< 3 %` |
| 65 | cortante aislado, `K6ROT=0` | `5.74e-9 %` | `< 1e-6 %` |
| 67 | MID1/MID2/MID3 y `K6ROT=100` | `0.000298 %` | `< 0.001 %` |
| 69-71, 73 | MID3/G3/TST rectangulares | `< 7.3e-9 %` | `< 1e-6 %` |
| 72 | MID3 vacío, `K6ROT=100` | `0.000121 %` | `< 0.001 %` |
| 74-75 | geometría sesgada, con/sin MID3 | `0.05436 %` | `< 1 %` |
| 81 | paralelogramo oblicuo, propiedades completas | `0.0000103 %` | `< 0.00002 %` |
| 83 | paralelogramo oblicuo, MID3 aislado | `0.000314 %` | `< 0.0004 %` |
| 89 | paralelogramo oblicuo, MID3 vacío | `0.0000889 %` | `< 0.0002 %` |

El caso 54 es inclinado y alabeado: todos los GDL globales reciben
contribuciones. La comparación no elimina filas ni columnas y verifica la
corrección de fuerzas y momentos entre el plano medio y los nodos reales. El
resultado cumple el umbral del proyecto, pero debe describirse como
coincidencia aproximada, no como identidad numérica.

Los casos 65/69 y 67/72 forman pares con y sin `MID3`; 70, 71 y 73 cambian
solo `G3` o `TS/T`. Todos ellos verifican tanto el error global como los
términos individuales. Los F06 66 y 68 confirman que desplazamientos, cargas,
reacciones y resultantes recuperados coinciden dentro de sus tolerancias. El
error casi idéntico de 74 y 75 muestra que el `0.054 %` residual pertenece a la
aproximación geométrica irregular y no a la flexibilidad física de `MID3`.

Los casos 81/83/89 separan el paralelogramo oblicuo en propiedades completas,
MID3 aislado y MID3 vacío. Tras añadir el acoplamiento de las dos muestras de
cortante, el bloque `[T3,R1,R2]` baja a `3.14e-6`, `3.14e-6` y `8.89e-7`. El
test sigue comprobando membrana, placa/cortante y *drilling* por separado y
añade comparación término a término con tolerancia absoluta; por tanto, la
norma global no puede ocultar una submatriz pequeña incorrecta.

Los cuatro `xfail` E2E restantes tienen un alcance diferente:

- 62 contiene un CQUAD4 intencionadamente alabeado (`4.15 %` de la longitud
  característica fuera del plano), por lo que prueba la corrección aproximada
  entre el plano medio y los nodos físicos;
- 63 y 64 ensamblan CQUAD4 alabeados con propiedades discontinuas. Sus
  desviaciones fuera del plano están entre `2.3 %` y `7.7 %`, acumulando
  proyección, corrección de alabeo y recuperación multielemento.

Son casos de interés para cerrar la generalización 3D, pero no contradicen la
validación plana. Los casos coplanares 60, 61 y 76 se han cerrado tras aplicar
la interpolación *edge-string shear*: sus errores de desplazamiento en norma
son `0.055 %`, `0.031 %` y `0.0116 %`. Los F06 de 76-80 ya están
transcritos explícitamente en la suite. Los casos 76, 77 (bracket plegado), 79
(viga cajón) y 80 (modelo mixto CBAR/CTRIA3/CQUAD4) validan desplazamientos,
cargas, reacciones, fuerzas y tensiones. Solo el 78 sigue como
`xfail(strict=True)` entre los nuevos casos planos: su discrepancia se concentra
en rotaciones de una cubierta facetada sometida a un momento con componentes
normales a las facetas. El caso 82, geométricamente idéntico pero con
`PARAM,SNORM,0.0`, ya confirma la hipótesis: frente a su F06 el solver presenta
`0.0227 %` de error de desplazamientos y `0.00695 %` de reacciones, mientras
que los desplazamientos de los propios F06 78 y 82 difieren un `8.25 %`.
Por tanto, el `xfail` 78 documenta la ausencia de promediado de normales shell
en GRID compartidos, no un error del CQUAD4 plano evaluado con normal propia.

Los F06 84-88 amplían la cobertura a un panel inclinado multipropiedad, dos
modelos mixtos CQUAD4/CTRIA3/CBAR/CROD, desplazamientos impuestos y una malla
facetada con `SNORM=0`. Los errores vectoriales de desplazamientos son
`1.22e-5`, `4.82e-6`, `1.69e-6`, `6.50e-6` y `1.69e-8`; los de reacciones son
`2.88e-5`, `1.66e-6`, `4.72e-6`, `8.35e-6` y `3.96e-8`. En 84 los resultantes
`N/M` y tensiones quedan por debajo de `0.003 %`, pero un cortante centroidal
pequeño difiere `0.00977` (`72 %` relativo). El test lo documenta con una
tolerancia exclusiva de `QX/QY`, sin relajar el resto de magnitudes.

La incorporación de estos casos también corrigió la recuperación con `THETA`:
MSC Nastran informa las componentes shell en los ejes del elemento. `THETA`
orienta el material y no debe rotar la salida de un `MAT1` isótropo hacia los
ejes materiales.

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

1. Ejecutar los controles predictivos 90/91 de alta oblicuidad e incorporar su
   DMIG/F06 sin reajustar primero la formulación.
2. Reproducir la recuperación centroidal `QX/QY` comercial.
3. Reducir el error geométrico residual de los CQUAD4 alabeados.
4. Ampliar los estudios de convergencia y *patch tests* multielemento.
5. Documentar la versión exacta de Nastran que genera cada referencia.
6. Consolidar los casos CTRIA3 con diferencias pequeñas aún conocidas.
7. Incorporar CI para publicar automáticamente el resultado de la suite.
