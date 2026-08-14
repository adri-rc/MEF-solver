# Documentación de MEF-solver

MEF-solver es un solver académico de elementos finitos orientado al análisis estructural estático lineal y a la comparación con modelos de MSC Nastran. El proyecto dispone de un flujo `BDF/DAT -> modelo interno -> ensamblaje -> solución -> resultados`, pero no es un sustituto general de Nastran: implementa un subconjunto deliberadamente pequeño de tarjetas. La rigidez CQUAD4 se considera cerrada para cuadriláteros coplanares.

Esta documentación describe el estado observado en el repositorio el **14 de agosto de 2026**. La colección completa termina con `335 passed, 6 xfailed`.

## Índice

- [Documento de traspaso](TRASPASO.md): contexto autocontenido para continuar el proyecto en otra sesión u ordenador.
- [Uso del programa](uso.md): entorno, ejecución, formato de entrada y resultados.
- [Arquitectura](arquitectura.md): módulos, modelo de datos y flujo del análisis.
- [Formulación matemática](formulacion.md): ecuaciones, convenciones y formulación de los elementos.
- [Estado y validación](validacion.md): cobertura real, resultados de pruebas y casos Nastran.
- [Desarrollo y validación de CQUAD4](cquad4.md): estado detallado del elemento cuadrilateral.
- [Formulación completa de CQUAD4](cquad4_tesis.md): desarrollo matemático, parámetros, correspondencia con el código, evidencia y limitaciones.
- [Limitaciones y hoja de ruta](roadmap.md): deuda técnica y evolución propuesta hacia elementos 3D.

## Estado funcional resumido

| Capacidad | Lectura | Matriz de rigidez | Solución global | Recuperación de resultados | Validación observada |
|---|---:|---:|---:|---:|---|
| CROD + PROD + MAT1 | Sí | Sí | Sí | Axil y torsión | Unitaria y casos Nastran de extremo a extremo |
| CBAR + PBAR + MAT1 | Sí | Sí | Sí | Fuerzas de extremo y tensiones simplificadas | Unitaria y casos Nastran de extremo a extremo |
| CTRIA3 + PSHELL + MAT1 | Sí | Sí | Sí | Membrana, flexión y cortante | Amplia comparación de extremo a extremo y de `K_global` |
| CQUAD4 + PSHELL + MAT1 | Sí | Sí | Sí | Membrana, flexión, cortante y fibras | DMIG hasta alta oblicuidad y E2E; sólo coplanar |
| FORCE / MOMENT / LOAD | Sí | Sí | Sí | Sistemas de entrada transformados | Cubierto por pruebas |
| SPC / SPC1 / SPCADD | Sí | Sí | Sí | Reacciones y valores impuestos | Cubierto por pruebas |
| Sólidos 3D | No | No | No | No | Trabajo futuro |

CQUAD4 dispone de una reconstrucción Nastran y una formulación clásica Reissner-Mindlin en módulos separados. La primera incorpora estabilización *hourglass*, cortante de MacNeal y `K6ROT`; rechaza explícitamente el alabeo. Véanse [cquad4.md](cquad4.md) y el [documento de traspaso](TRASPASO.md).

## Alcance físico actual

- Análisis estático, lineal, de pequeñas deformaciones.
- Material elástico lineal isótropo (`MAT1`).
- Seis grados de libertad por nodo: `[Tx, Ty, Tz, Rx, Ry, Rz]`.
- Matriz global dispersa y resolución directa mediante SciPy.
- Cargas nodales y restricciones de punto único.
- Sin análisis dinámico, no linealidad, contacto, plasticidad, cargas distribuidas, MPC/RBE ni sólidos.

## Lectura recomendada

Para ejecutar el programa, comenzar por [uso.md](uso.md). Para continuar el
proyecto desde el punto actual, leer [cquad4.md](cquad4.md), la formulación
autocontenida [cquad4_tesis.md](cquad4_tesis.md) y después
[validacion.md](validacion.md). La [hoja de ruta](roadmap.md) propone un orden
de trabajo que prioriza la corrección y la trazabilidad numérica.
