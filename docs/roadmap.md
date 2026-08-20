# Hoja de ruta

## Estado de base

El flujo DAT/BDF → modelo → ensamblaje disperso → solución → resultados está
operativo para CROD, CBAR, CTRIA3 y CQUAD4 plano. La colección completa termina
con `418 passed, 6 xfailed`; los xfail documentan geometría CQUAD4 alabeada y
suavizado `SNORM`, ambos fuera de alcance.

La rigidez CQUAD4 plana se considera cerrada tras los casos 90/91. No se
planifica extenderla a elementos alabeados en la siguiente fase.

Desde ese cierre se separó la escritura del informe en
`src/postprocessing/report_writer.py`, se añadieron sus pruebas y se verificó
el caso mixto 80 con CBAR, CTRIA3 y CQUAD4. `requirements.txt` declara
dependencias mínimas, aunque todavía falta un bloqueo reproducible.

El preproceso distingue estática lineal y `SOL 400`, normaliza `NLSTEP` o
`NLPARM`, `LGDISP` y `FOLLOWK`, y dispone de 31 casos CROD (92--122) ejecutados
con MSC Nastran 2026.1/NLPERF. `NonlinearStaticSolver` resuelve el subconjunto
MAT1/PROD/CROD mediante Newton completo, tangente consistente, stepping fijo o
adaptativo y recortes. Los 30 trabajos convergidos coinciden con referencias
explícitas y el 120 reproduce el límite bajo control de fuerza.

## Prioridad 0 completada: no linealidad geométrica CROD

1. Completada la campaña MSC y conservadas sus salidas auditables.
2. Transcritos estados intermedios, tangente por sonda y resultados finales.
3. Implementados `CRODNonlinear` y `NonlinearStaticSolver` sin modificar la
   formulación lineal.
4. Integrados CLI, historial incremental, recuperación e informe.
5. Documentado el punto límite 120 y el alcance sin arco de longitud.

La siguiente ampliación física será CBAR con grandes desplazamientos, como una
fase independiente y precedida por sus propios casos de referencia.

## Prioridad 1: consolidación

1. Añadir integración continua para ejecutar pytest en cada cambio.
2. Añadir `pyproject.toml` o un fichero de bloqueo con versiones reproducibles.
3. Añadir pruebas de convergencia sobre mallas refinadas, además de DMIG/E2E.
4. Definir una política de ejecución por lotes: detenerse o continuar y
   acumular errores por modelo.
5. Revisar y aislar los módulos y scripts antiguos o exploratorios antes de
   ampliar el alcance físico.

## Prioridad 2: rendimiento

1. Ensamblar mediante tripletas COO vectorizadas.
2. Cachear transformaciones y objetos de elemento entre ensamblaje y
   recuperación.
3. Medir lectura, rigidez elemental, ensamblaje, solución y postproceso sobre
   familias de mallas crecientes.
4. Estudiar reordenamientos y solvers iterativos cuando el tamaño lo justifique.

## Prioridad 3: nuevas capacidades

Cada ampliación debe incorporar un caso mínimo, una referencia externa y una
prueba E2E antes de anunciarse como soportada. Candidatos:

- presión de lámina, gravedad y otras cargas distribuidas;
- MPC y elementos rígidos;
- propiedades avanzadas de vigas;
- materiales ortótropos y laminados;
- elementos de mayor orden y sólidos;
- exportación a formatos de visualización.

## Mejoras opcionales de CQUAD4

- Investigar la recuperación centroidal propietaria de `Qx/Qy` sin tocar la
  rigidez validada.
- Añadir resultados en más estaciones o vértices.
- Implementar `SNORM` sólo si aparecen modelos facetados que lo requieran.
- Mantener el rechazo de nodos no coplanares. Una futura formulación alabeada
  sería un elemento nuevo, no una relajación de la tolerancia actual.
