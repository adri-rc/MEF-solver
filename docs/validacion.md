# Estado y validación

## Línea base reproducida

Fecha: **5 de agosto de 2026**. Entorno: `.venv` del repositorio.

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

Resultado:

```text
298 passed in 1.88s
```

También se ejecutó `python -m compileall -q src tests main.py` sin errores. La
suite completa queda limpia; ya no es necesario excluir pruebas históricas.

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
- CQUAD4: dos formulaciones separadas, postproceso, `K6ROT`, modos rígidos y
  comparación automática de los casos 51 a 54.

## CQUAD4 frente a Nastran

| Caso | Comparación | Error de Frobenius | Criterio |
|---|---|---:|---:|
| 51 | bloque de membrana | `0.152382 %` | `< 1 %` |
| 52 | bloque de flexión | `0.152382 %` | `< 1 %` |
| 53 | bloque activo combinado | `0.372560 %` | `< 1 %` |
| 54 | matriz global completa `24 x 24` | `2.150270 %` | `< 3 %` |

El caso 54 es inclinado y alabeado: todos los GDL globales reciben
contribuciones. La comparación no elimina filas ni columnas y verifica la
corrección de fuerzas y momentos entre el plano medio y los nodos reales. El
resultado cumple el umbral del proyecto, pero debe describirse como
coincidencia aproximada, no como identidad numérica.

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
matriz completa se informa además del máximo absoluto, pero el criterio global
principal es el error relativo de Frobenius.

## Trabajo pendiente de validación

1. Añadir más DMIG CQUAD4 alabeadas y con cortante transversal no nulo.
2. Ampliar los estudios de convergencia y *patch tests* multielemento.
3. Documentar la versión exacta de Nastran que genera cada referencia.
4. Consolidar los casos CTRIA3 con diferencias pequeñas aún conocidas.
5. Incorporar CI para publicar automáticamente el resultado de la suite.
