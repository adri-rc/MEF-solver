# Documentación de MEF-solver

MEF-solver es un solver académico de elementos finitos orientado al análisis estructural estático y a la comparación con modelos de MSC Nastran. Incluye estática lineal para CROD, CBAR, CTRIA3 y CQUAD4 coplanar, y estática no lineal geométrica `SOL 400` para MAT1 con CROD o CBAR. El proyecto no es un sustituto general de Nastran: implementa un subconjunto deliberadamente pequeño de tarjetas. La rigidez CQUAD4 se considera cerrada para cuadriláteros coplanares.

Esta documentación describe el estado verificado en el repositorio el **21 de agosto de 2026**. La colección completa termina con `570 passed, 7 xfailed`.

## Índice

- [Documento de traspaso](TRASPASO.md): contexto autocontenido para continuar el proyecto en otra sesión u ordenador.
- [Uso del programa](uso.md): entorno, ejecución, formato de entrada y resultados.
- [Guía del informe RPT](result_report.md): significado de cada tabla, correspondencia F06 e historial de Newton.
- [Arquitectura](arquitectura.md): módulos, modelo de datos y flujo del análisis.
- [Formulación matemática](formulacion.md): ecuaciones, convenciones y formulación de los elementos.
- [Estado y validación](validacion.md): cobertura real, resultados de pruebas y casos Nastran.
- [Desarrollo y validación de CQUAD4](cquad4.md): estado detallado del elemento cuadrilateral.
- [Formulación completa de CQUAD4](cquad4_tesis.md): desarrollo matemático, parámetros, correspondencia con el código, evidencia y limitaciones.
- [Análisis no lineal de CROD](nonlinear_crod.md): formulación, tarjetas SOL 400, solver incremental y referencias 92--127.
- [Análisis no lineal de CBAR](nonlinear_cbar.md): formulación corrotacional, límites y campaña 128--153.
- [Limitaciones y hoja de ruta](roadmap.md): deuda técnica y evolución propuesta hacia elementos 3D.

## Estado funcional resumido

| Capacidad | Lectura | Matriz de rigidez | Solución global | Recuperación de resultados | Validación observada |
|---|---:|---:|---:|---:|---|
| CROD + PROD + MAT1 | Sí | Sí | Sí | Axil y torsión | Unitaria y casos Nastran de extremo a extremo |
| CBAR + PBAR + MAT1 | Sí | Sí | Sí | Fuerzas de extremo y tensiones simplificadas | Unitaria y casos Nastran de extremo a extremo |
| CTRIA3 + PSHELL + MAT1 | Sí | Sí | Sí | Membrana, flexión y cortante | Amplia comparación de extremo a extremo y de `K_global` |
| CQUAD4 + PSHELL + MAT1 | Sí | Sí | Sí | Membrana, flexión, cortante y fibras | DMIG hasta alta oblicuidad y E2E; sólo coplanar |
| SOL 400 + NLSTEP/NLPARM | Sí | Sí, tangente CROD/CBAR | Sí | Historial incremental y resultados CROD/CBAR | Campañas 92--153 de MSC |
| FORCE / MOMENT / LOAD | Sí | Sí | Sí | Sistemas de entrada transformados | Cubierto por pruebas |
| SPC / SPC1 / SPCADD | Sí | Sí | Sí | Reacciones y valores impuestos | Cubierto por pruebas |
| Sólidos 3D | No | No | No | No | Trabajo futuro |

CQUAD4 dispone de una reconstrucción Nastran y una formulación clásica Reissner-Mindlin en módulos separados. La primera incorpora estabilización *hourglass*, cortante de MacNeal y `K6ROT`; rechaza explícitamente el alabeo. Véanse [cquad4.md](cquad4.md) y el [documento de traspaso](TRASPASO.md).

## Alcance físico actual

- Solución estática lineal de pequeñas deformaciones.
- Estática no lineal geométrica CROD/CBAR con `SOL 400`, `LGDISP`, Newton completo y controles `FIXED`/`ADAPT`.
- Material elástico lineal isótropo (`MAT1`).
- Seis grados de libertad por nodo: `[Tx, Ty, Tz, Rx, Ry, Rz]`.
- Matriz global dispersa y resolución directa mediante SciPy.
- Cargas nodales y restricciones de punto único.
- Sin dinámica, no linealidad material, contacto, cargas seguidoras, cargas distribuidas, MPC/RBE ni sólidos; `PBARN1,INT2=LS` no se reproduce exactamente.

## Lectura recomendada

Para ejecutar el programa, comenzar por [uso.md](uso.md). Para continuar el
proyecto desde el punto actual, empezar por [TRASPASO.md](TRASPASO.md), seguir
el contexto técnico que allí se recomienda y respetar las reglas generales de
[`AGENTS.md`](../AGENTS.md). La [hoja de ruta](roadmap.md) recoge el trabajo
pendiente sin reabrir la rigidez CQUAD4 coplanar ya cerrada.
