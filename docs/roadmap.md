# Limitaciones y hoja de ruta

## Estado de partida

El flujo estático lineal y la suite de 337 pruebas están operativos: 333 pasan
y cuatro discrepancias conocidas están marcadas `xfail(strict=True)`. CROD,
CBAR, CTRIA3 y CQUAD4 se leen, ensamblan y recuperan resultados. CQUAD4 dispone
de formulaciones Nastran y Reissner-Mindlin separadas; los casos DMIG 51-55
cumplen sus tolerancias documentadas. Los casos 65-75 validan la flexibilidad
transversal `[T3,R1,R2]`, incluida la semántica de `MID3`, `G3` y `TS/T`.

El proyecto sigue siendo un solver académico de alcance controlado, no un
sustituto general de Nastran.

## Limitaciones principales

- solo estática lineal, pequeñas deformaciones y material `MAT1` isótropo;
- sin cargas distribuidas, gravedad, temperatura, MPC/RBE ni contacto;
- una ejecución resuelve un único subcaso seleccionado;
- campos shell como `MID4`, `Z1/Z2`, `ZOFFS`, `TFLAG` y espesores nodales no
  tienen todavía tratamiento físico completo;
- CQUAD4 alabeado coincide con la DMIG del caso 54 al `2.08 %`, no término a
  término, y necesita más referencias independientes;
- los CQUAD4 planos sesgados no afines 74-75 conservan un error de matriz del
  `0.046 %`; los afines oblicuos 81/83/89 quedan por debajo de `3.2e-6` de
  error relativo tras introducir el acoplamiento torsional de muestras;
- la recuperación centroidal `QX/QY` conserva diferencias localizadas aun
  cuando rigidez, desplazamientos, reacciones, `N/M` y tensiones coinciden;
- el solver no reproduce el promediado de normales shell de `SNORM` en GRID
  compartidos; con `SNORM=0`, el control 82 sí coincide con Nastran;
- la geometría alabeada necesita más referencias independientes;
- no hay sólidos 3D, dinámica, plasticidad ni no linealidad geométrica;
- el informe de texto aún puede especializarse mejor por tipo de elemento;
- `assemble()` no es idempotente y no existe CI publicada.

## Prioridades

### P0 - Consolidar la validación 2D

1. Refinar la cinemática geométrica de cuadriláteros sesgados/alabeados sin
   degradar las DMIG rectangulares 65-73 ya cerradas término a término.
2. Ejecutar *patch tests* multielemento y convergencia de placas.
3. Cerrar las discrepancias CTRIA3 conocidas con métricas absolutas y
   relativas separadas.
4. Registrar versión, deck y procedencia de cada DMIG/F06.
5. Ejecutar 90/91 para validar, sin recalibración previa, la continuación
   angular acotada en una oblicuidad mayor que 75/89.
6. Derivar la recuperación centroidal comercial de `QX/QY`; el caso 84 aísla
   una diferencia absoluta máxima de `0.00977` sin error apreciable en `N/M`.

### P1 - Completar semántica shell

1. Implementar `MID4`, `Z1/Z2`, `ZOFFS`, `TFLAG` y `Ti`; `12I/T^3` y `TS/T`
   ya se leen y participan en la rigidez/recuperación CQUAD4.
2. Validar materiales ortótropos y laminados antes de generalizar la
   estabilización isotrópica actual.
3. Separar recuperación y escritores en `src/postprocessing/`.

### P2 - Robustez y automatización

1. Añadir CI con pytest, `compileall` y resumen numérico de validación.
2. Hacer `assemble()` idempotente o impedir explícitamente la segunda llamada.
3. Mejorar diagnósticos de condicionamiento, mecanismos y componentes
   desconectadas.
4. Perfilar el ensamblaje antes de optimizarlo.

### P3 - Nuevas capacidades

1. Cargas distribuidas, gravedad, temperatura y restricciones multipunto.
2. Arquitectura para resolver varios subcasos en una ejecución.
3. `CTETRA4` como primer sólido, seguido de `CHEXA8` con *patch tests* 3D.
4. Solo después, dinámica y no linealidad si entran en el alcance del TFM.

## Entregables recomendados para la tesis

- alcance y subconjunto Nastran soportado a nivel de campo;
- derivación y comparación de las dos formulaciones CQUAD4;
- matriz de trazabilidad requisito-prueba-caso;
- resultados de DMIG, *patch tests*, convergencia y sensibilidad geométrica;
- explicación explícita de las diferencias con Nastran;
- caso demostrador reproducible desde el deck hasta el informe.
