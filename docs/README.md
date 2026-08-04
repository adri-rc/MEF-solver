# Documentación de MEF-solver

MEF-solver es un solver académico de elementos finitos orientado al análisis estructural estático lineal y a la comparación con modelos de MSC Nastran. El proyecto ya dispone de un flujo `BDF/DAT -> modelo interno -> ensamblaje -> solución -> resultados`, pero no debe considerarse todavía un sustituto general de Nastran: implementa un subconjunto deliberadamente pequeño de tarjetas y varias capacidades, especialmente CQUAD4, siguen en validación.

Esta documentación describe el estado observado en el repositorio el **4 de agosto de 2026**. Distingue entre una tarjeta que puede leerse, un elemento cuya rigidez puede ensamblarse y una capacidad validada de extremo a extremo.

## Índice

- [Documento de traspaso](TRASPASO.md): contexto autocontenido para continuar el proyecto en otra sesión u ordenador.
- [Uso del programa](uso.md): entorno, ejecución, formato de entrada y resultados.
- [Arquitectura](arquitectura.md): módulos, modelo de datos y flujo del análisis.
- [Formulación matemática](formulacion.md): ecuaciones, convenciones y formulación de los elementos.
- [Estado y validación](validacion.md): cobertura real, resultados de pruebas y casos Nastran.
- [Desarrollo y validación de CQUAD4](cquad4.md): estado detallado del elemento cuadrilateral.
- [Limitaciones y hoja de ruta](roadmap.md): deuda técnica y evolución propuesta hacia elementos 3D.

## Estado funcional resumido

| Capacidad | Lectura | Matriz de rigidez | Solución global | Recuperación de resultados | Validación observada |
|---|---:|---:|---:|---:|---|
| CROD + PROD + MAT1 | Sí | Sí | Sí | Axil y torsión | Unitaria y casos Nastran de extremo a extremo |
| CBAR + PBAR + MAT1 | Sí | Sí | Sí | Fuerzas de extremo y tensiones simplificadas | Unitaria y casos Nastran de extremo a extremo |
| CTRIA3 + PSHELL + MAT1 | Sí | Sí | Sí | Membrana, flexión y cortante | Amplia comparación de extremo a extremo y de `K_global` |
| CQUAD4 + PSHELL + MAT1 | Sí | Sí | Sí | Membrana, flexión y cortante | Integración global y comparación de los bloques activos de los casos Nastran 51/52 con error normalizado inferior al 3 % |
| FORCE / MOMENT | Sí | Sí | Sí | En coordenadas globales | Cubierto por pruebas |
| SPC / SPC1 | Sí | Sí | Sí | Reacciones | Cubierto por pruebas |
| Sólidos 3D | No | No | No | No | Trabajo futuro |

La formulación CQUAD4 sigue en validación: incorpora integración reducida selectiva, estabilización *hourglass* y *drilling*, y ya participa en el postproceso global. La comparación completa depende de la semántica deliberada de sustitución automática de los MID vacíos de `PSHELL`; véanse [cquad4.md](cquad4.md) y el [documento de traspaso](TRASPASO.md).

## Alcance físico actual

- Análisis estático, lineal, de pequeñas deformaciones.
- Material elástico lineal isótropo (`MAT1`).
- Seis grados de libertad por nodo: `[Tx, Ty, Tz, Rx, Ry, Rz]`.
- Matriz global dispersa y resolución directa mediante SciPy.
- Cargas nodales y restricciones de punto único.
- Sin análisis dinámico, no linealidad, contacto, plasticidad, cargas distribuidas, MPC/RBE ni sólidos.

## Lectura recomendada

Para ejecutar el programa, comenzar por [uso.md](uso.md). Para continuar el TFM desde el punto actual, leer primero [cquad4.md](cquad4.md) y después [validacion.md](validacion.md). La [hoja de ruta](roadmap.md) propone un orden de trabajo que prioriza la corrección y la trazabilidad numérica antes de incorporar elementos 3D.
