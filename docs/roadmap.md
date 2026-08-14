# Hoja de ruta

## Estado de base

El flujo DAT/BDF → modelo → ensamblaje disperso → solución → resultados está
operativo para CROD, CBAR, CTRIA3 y CQUAD4 plano. La colección completa termina
con `337 passed, 6 xfailed`; los xfail documentan geometría CQUAD4 alabeada y
suavizado `SNORM`, ambos fuera de alcance.

La rigidez CQUAD4 plana se considera cerrada tras los casos 90/91. No se
planifica extenderla a elementos alabeados en la siguiente fase.

## Prioridad 1: consolidación

1. Añadir integración continua para ejecutar pytest en cada cambio.
2. Mantener un entorno reproducible y fijar versiones compatibles.
3. Separar el postproceso de `LinearStaticSolver`.
4. Completar el informe RPT con todas las magnitudes ya calculadas.
5. Añadir pruebas de convergencia sobre mallas refinadas, además de DMIG/E2E.

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
