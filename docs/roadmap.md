# Hoja de ruta

## Estado de base

El flujo DAT/BDF → modelo → ensamblaje disperso → solución → resultados está
operativo para CROD, CBAR, CTRIA3 y CQUAD4 plano. La colección completa termina
con `570 passed, 7 xfailed`; seis límites corresponden a CQUAD4 y el séptimo a
`PBARN1,INT2=LS`.

La rigidez CQUAD4 plana se considera cerrada tras los casos 90/91. No se
planifica extenderla a elementos alabeados en la siguiente fase.

Desde ese cierre se separó la escritura del informe en
`src/postprocessing/report_writer.py`, se añadieron sus pruebas y se verificó
el caso mixto 80 con CBAR, CTRIA3 y CQUAD4. `requirements.txt` declara
dependencias mínimas, aunque todavía falta un bloqueo reproducible.

El preproceso distingue estática lineal y `SOL 400`, normaliza `NLSTEP` o
`NLPARM`, `LGDISP`, `FOLLOWK` y los campos relevantes de `PBARN1`.
`NonlinearStaticSolver` resuelve MAT1 con CROD o CBAR mediante Newton completo,
stepping fijo/adaptativo y recortes. CROD dispone de 36 casos (92--127) y CBAR
de 26 (128--153), todos con referencias MSC auditables.

## Prioridad 0 completada: no linealidad geométrica CROD

1. Completada la campaña MSC y conservadas sus salidas auditables.
2. Transcritos estados intermedios, tangente por sonda y resultados finales.
3. Implementados `CRODNonlinear` y `NonlinearStaticSolver` sin modificar la
   formulación lineal.
4. Integrados CLI, recuperación e informe con configuración, historial de
   incrementos, correcciones de Newton y estado final diferenciado.
5. Documentado el punto límite 120 y el alcance sin arco de longitud.

## Prioridad 0 completada: no linealidad geométrica CBAR

1. Implementado `CBARNonlinear` corrotacional sin modificar la ruta lineal.
2. Integrados `I12`, `K1/K2`, offsets, releases, orientación `G0`, cargas y
   SPC impuestos.
3. Identificada la deformación axial logarítmica con 130/131 y validada la
   tangente mediante las sondas 137--139.
4. Transcritos 26 E2E explícitos: 25 pasan y 142 es `xfail` estricto.
5. Limitación abierta: reproducir la interpolación propietaria
   `PBARN1,INT2=LS` mediante los diagnósticos 154--160 y cerrar tensiones de
   fibra bajo rotación 3D combinada.

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
